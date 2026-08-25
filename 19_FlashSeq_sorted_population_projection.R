#!/usr/bin/env Rscript

#' FLASH-seq Sorted Population Projection
#'
#' Pilot validation of the CD27/CD28 surface marker strategy. CD27-high/CD28-high and
#' CD27-low/CD28-low malignant T cells were sorted from one patient into a 96-well plate and
#' profiled with the MERCURIUS FLASH-seq full-length protocol (Alithea Genomics). This script
#' builds the Seurat object from featureCounts output, applies QC, and projects the sorted cells
#' onto the CTCL and healthy-donor blood reference by Seurat label transfer.
#'
#' Corresponds to Suppl. Figure 5F-H.
#'
#' Note on the data: FLASH-seq is a full-length protocol without UMIs, so gene-level counts are
#' reads assigned to genes, not unique molecules. All depth-sensitive steps account for this.
#'
#' Dependencies: Seurat, qs2, dplyr, readr

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(readr)
})

param <- list(seed = 38, nfeatures = 3000, npcs_ref = 50, dims = 1:30)
set.seed(param$seed)
options(future.globals.maxSize = 32 * 1024^3)

# ---------------------------------------------------------------------------------------------
# Plate layout
# ---------------------------------------------------------------------------------------------

#' Build the well index to sorted population map
#'
#' 8 rows (A-H) x 12 columns, column-major well numbering: well 1 = A1, 2 = B1, ... 8 = H1,
#' 9 = A2, ... 96 = H12. Wells A1 (index 1) and A7 (index 49) were left empty.
#'
#' @param first_half Population sorted into plate columns 1-6. For the run analysed here
#'   (o41017) columns 1-6 are mTRM and columns 7-12 are mTCM, matching plate_layout.xlsx.
#' @return data.frame with well_num, col and population
#' @details The column-to-population assignment is NOT constant across runs. On a later run the
#'   sorted populations were dispensed into the opposite column halves relative to the layout
#'   file, which was only caught from the marker phenotype. This is therefore an explicit
#'   argument, and verify_orientation() below re-derives it from the data. Do not treat either
#'   as optional.
build_plate_map <- function(first_half = c("mTRM", "mTCM")) {
  first_half <- match.arg(first_half)
  other <- setdiff(c("mTRM", "mTCM"), first_half)

  # 8 rows x 12 columns, column-major: well 1 = A1, 8 = H1, 9 = A2, ... 96 = H12.
  map <- expand.grid(row = LETTERS[1:8], col = 1:12, stringsAsFactors = FALSE) |>
    dplyr::arrange(col, row) |>
    dplyr::mutate(
      well_num   = dplyr::row_number(),
      population = ifelse(col <= 6, first_half, other)
    )

  # Wells A1 (index 1) and A7 (index 49) were left empty; they are dropped downstream.
  map$population[map$well_num %in% c(1, 49)] <- NA_character_
  map
}

#' Verify the plate orientation from the marker phenotype
#'
#' The sort was on CD27 and CD28, so the wells labelled mTCM must show higher central-memory
#' marker expression than those labelled mTRM. If they do not, the column halves are swapped
#' and every downstream result is mislabelled. This stops rather than warns, because the failure
#' is silent and produces a plausible-looking but inverted answer.
#'
#' @param obj Seurat object with a population column and a log-normalised RNA assay
#' @param markers Central-memory markers expected to be higher in mTCM
verify_orientation <- function(obj, markers = c("CD27", "CCR7", "SELL", "TCF7", "LEF1")) {
  DefaultAssay(obj) <- "RNA"
  present <- intersect(markers, rownames(obj))
  if (length(present) < 2) {
    stop("Cannot verify plate orientation: fewer than two central-memory markers detected.")
  }

  score <- Matrix::colMeans(GetAssayData(obj, layer = "data")[present, , drop = FALSE])
  m_tcm <- median(score[obj$population == "mTCM"])
  m_trm <- median(score[obj$population == "mTRM"])

  cat(sprintf("Orientation check on %s: median score mTCM = %.3f, mTRM = %.3f\n",
              paste(present, collapse = "/"), m_tcm, m_trm))
  if (!(m_tcm > m_trm)) {
    stop("Plate orientation looks inverted: the wells labelled mTCM do not show higher ",
         "central-memory marker expression. Re-check the column-to-population assignment ",
         "with build_plate_map(first_half = ...) before going further.")
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------------------------
# Build the Seurat object from featureCounts
# ---------------------------------------------------------------------------------------------

#' Assemble a per-cell count matrix from featureCounts output
#'
#' Each per-cell file is a two-column TSV (Ensembl gene id, count). Cells are joined by gene
#' and missing entries filled with zero. Ensembl ids are converted to gene symbols because the
#' reference is symbol-based; the Ensembl id is kept where no symbol exists, and collisions are
#' resolved with make.unique.
#'
#' @param fc_dir featureCounts output directory for the order
#' @param annot_tsv Gene annotation (GENCODE GRCh38.p14 Release 48) with gene_id and gene_name
#' @return sparse count matrix, genes x cells
build_count_matrix <- function(fc_dir, annot_tsv) {
  files <- list.files(fc_dir, pattern = "\\.txt$", full.names = TRUE)
  files <- files[!grepl("-stats\\.txt$", files)]
  stopifnot(length(files) > 0)

  per_cell <- lapply(files, function(f) {
    d <- readr::read_tsv(f, show_col_types = FALSE)
    setNames(d[[2]], d[[1]])
  })
  names(per_cell) <- sub("\\.txt$", "", basename(files))

  genes <- sort(unique(unlist(lapply(per_cell, names))))
  mat <- vapply(per_cell, function(v) {
    out <- setNames(numeric(length(genes)), genes)
    out[names(v)] <- v
    out
  }, numeric(length(genes)))

  annot <- readr::read_tsv(annot_tsv, show_col_types = FALSE)
  sym <- annot$gene_name[match(rownames(mat), annot$gene_id)]
  sym[is.na(sym) | sym == ""] <- rownames(mat)[is.na(sym) | sym == ""]
  rownames(mat) <- make.unique(sym)

  Matrix::Matrix(mat, sparse = TRUE)
}

#' Build and QC-filter the FLASH-seq Seurat object
#'
#' @param count_mat Output of build_count_matrix
#' @param plate_map Output of build_plate_map
#' @return QC-filtered, SCTransform-normalised Seurat object
build_flashseq_object <- function(count_mat, plate_map) {
  # Well index is the integer immediately after the first underscore in the sample name.
  well <- as.integer(sub("^[^_]*_(\\d+).*$", "\\1", colnames(count_mat)))
  stopifnot(!any(is.na(well)))

  obj <- CreateSeuratObject(counts = count_mat, project = "FlashSeq")
  obj$well_num   <- well
  obj$population <- plate_map$population[match(well, plate_map$well_num)]

  # Drop the deliberately empty wells.
  obj <- obj[, !is.na(obj$population)]

  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  # QC. The count threshold is what removes wells that received no cell.
  obj <- subset(obj, nFeature_RNA >= 200 & nCount_RNA >= 1000 & percent_mt < 25)
  cat("Cells passing QC:", ncol(obj), "\n")
  print(table(obj$population))

  # SCTransform regressing out total counts. Without UMIs, read depth is the dominant
  # technical axis in a full-length library, so it must be regressed rather than ignored.
  obj <- SCTransform(obj, vars.to.regress = "nCount_RNA", verbose = FALSE)

  # A log-normalised RNA assay is required in parallel: label transfer onto the 10x-based
  # reference is cross-protocol and uses LogNormalize, not SCT.
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = param$nfeatures, verbose = FALSE)

  # Hard stop if the column halves are swapped. See verify_orientation().
  verify_orientation(obj)
  obj
}

# ---------------------------------------------------------------------------------------------
# Projection onto the blood reference
# ---------------------------------------------------------------------------------------------

#' Prepare the reference for anchoring and projection
#'
#' The UMAP model must be fitted on the SAME reduction the anchors use, with
#' return.model = TRUE, or MapQuery silently produces no projection.
#'
#' @param ref_path Path to CTCL_HD_blood_global.qs2
#' @return Seurat reference with pca and a projection-ready umap
prepare_reference <- function(ref_path) {
  ref <- qs2::qs_read(ref_path, nthreads = 1)
  DefaultAssay(ref) <- "RNA"

  # Build the composite label. This is load-bearing, not cosmetic: Malignant_T is NA for every
  # non-malignant reference cell, so transferring it alone forces every query cell onto a
  # malignant state and the sorted CD27-low/CD28-low cells cannot come back as a normal T-cell
  # subset. The composite label is the malignant state where a cell has one, and the detailed
  # cell type otherwise.
  detail_col <- NA_character_
  for (dc in c("predicted.celltype.l2", "celltype", "Global_celltype",
               "celltype.l2", "seurat_clusters")) {
    if (dc %in% colnames(ref@meta.data)) { detail_col <- dc; break }
  }
  if (is.na(detail_col) || !"Malignant_T" %in% colnames(ref@meta.data)) {
    stop("Reference lacks Malignant_T and/or a detailed cell-type column; ",
         "the composite label cannot be built and label transfer would be degenerate.")
  }
  cat("Detail label column:", detail_col, "\n")

  mal    <- as.character(ref$Malignant_T)
  detail <- as.character(ref@meta.data[[detail_col]])
  ref$composite_label <- dplyr::case_when(
    !is.na(mal) & mal != "healthyT" ~ mal,
    TRUE                            ~ detail
  )
  ref$detail_label <- detail

  n_na <- sum(is.na(ref$composite_label))
  cat("Composite label:", length(unique(na.omit(ref$composite_label))), "categories,",
      n_na, "NA\n")
  # A composite label that is mostly one value means the fallback did not fire.
  if (max(table(ref$composite_label)) > 0.9 * ncol(ref)) {
    stop("Composite label is >90% a single value; the malignant/detail merge is wrong.")
  }
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, nfeatures = param$nfeatures, verbose = FALSE)
  ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
  ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = param$npcs_ref, verbose = FALSE)
  ref <- RunUMAP(ref, reduction = "pca", dims = param$dims,
                 return.model = TRUE, verbose = FALSE)
  ref
}

#' Transfer reference labels onto the sorted cells and project into the reference UMAP
#'
#' Three label sets are transferred: the malignant state, the level-2 cell type, and a
#' composite label that uses the malignant state for malignant cells and the level-2 cell type
#' for everything else. The composite label is the one reported, because a sorted cell that is
#' not malignant must be allowed to receive a non-malignant identity.
#'
#' @param query QC-filtered FLASH-seq object
#' @param ref Output of prepare_reference
#' @return query with predicted labels and a ref.umap embedding
project_onto_reference <- function(query, ref) {
  anchors <- FindTransferAnchors(
    reference           = ref,
    query               = query,
    normalization.method = "LogNormalize",
    reference.reduction  = "pca",
    dims                 = param$dims,
    verbose              = FALSE
  )

  # composite_label is the reported one and must be present; prepare_reference() guarantees it.
  stopifnot("composite_label" %in% colnames(ref@meta.data))
  refdata <- list(
    composite_label = ref$composite_label,
    Malignant_T     = ref$Malignant_T,
    detail_label    = ref$detail_label
  )
  refdata <- lapply(refdata, function(v) setNames(as.character(v), colnames(ref)))

  query <- TransferData(anchorset = anchors, refdata = refdata,
                        query = query, dims = param$dims, verbose = FALSE)

  query <- MapQuery(anchorset = anchors, reference = ref, query = query,
                    refdata = refdata, reference.reduction = "pca",
                    reduction.model = "umap", verbose = FALSE)
  stopifnot(length(query@reductions) > 0)
  query
}

# ---------------------------------------------------------------------------------------------
# Differential expression between the sorted populations
# ---------------------------------------------------------------------------------------------

#' Differential expression between the two sorted populations
#'
#' @param obj Projected FLASH-seq object
#' @return FindMarkers data.frame, positive log2FC = higher in mTRM
run_sorted_de <- function(obj) {
  DefaultAssay(obj) <- "SCT"
  obj <- PrepSCTFindMarkers(obj)
  FindMarkers(obj, group.by = "population", ident.1 = "mTRM", ident.2 = "mTCM",
              test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0.1)
}

# ---------------------------------------------------------------------------------------------
# Workflow
# ---------------------------------------------------------------------------------------------

#' @param fc_dir featureCounts directory for the FLASH-seq order
#' @param annot_tsv GENCODE GRCh38.p14 Release 48 gene annotation
#' @param first_half Population in plate columns 1-6 ("mTRM" for o41017)
#' @param ref_path CTCL_HD_blood_global.qs2
#' @param out_dir Output directory
main <- function(fc_dir, annot_tsv, first_half, ref_path, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  plate_map <- build_plate_map(first_half)
  count_mat <- build_count_matrix(fc_dir, annot_tsv)
  obj       <- build_flashseq_object(count_mat, plate_map)

  ref <- prepare_reference(ref_path)
  obj <- project_onto_reference(obj, ref)
  rm(ref); gc()

  # Sorted population versus assigned label. This cross-tabulation IS the result of
  # Suppl. Figure 5F-G, so print it rather than only plotting it.
  pred_col <- grep("composite_label$", colnames(obj@meta.data), value = TRUE)
  pred_col <- pred_col[!grepl("score$", pred_col)][1]
  if (is.na(pred_col)) stop("No transferred composite label found on the query object.")
  ct <- table(obj@meta.data[[pred_col]], obj$population)
  cat("\nSorted population versus transferred label:\n"); print(ct)
  # If every sorted cell received the same label, the transfer collapsed and the result is
  # not interpretable. Fail loudly rather than writing a plausible-looking table.
  if (nrow(ct) < 2) {
    stop("All query cells received a single label; label transfer collapsed.")
  }

  degs <- run_sorted_de(obj)

  qs2::qs_save(obj, file.path(out_dir, "obj_flashseq_mapped.qs2"), nthreads = 8)
  write.csv(degs, file.path(out_dir, "DEGs_sorted_mTRM_vs_mTCM.csv"))
  cat("\nDone.\n")
}

if (sys.nframe() == 0 && !interactive()) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 5) {
    stop("Usage: 19_FlashSeq_sorted_population_projection.R ",
         "<fc_dir> <annot_tsv> <first_half:mTRM|mTCM> <ref_qs2> <out_dir>")
  }
  main(a[1], a[2], a[3], a[4], a[5])
}
