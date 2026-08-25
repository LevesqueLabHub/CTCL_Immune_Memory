#!/usr/bin/env Rscript

#' TCR Recovery and Library Complexity by Malignant State
#'
#' Quantifies how often a productive TCR clonotype is recovered from each malignant state, and
#' asks whether the difference is explained by library complexity rather than by biology.
#' The motivating observation is that mTRM cells yield productive clonotypes less often than
#' mTCM cells, which could trivially reflect shallower or less complex libraries.
#'
#' Corresponds to Suppl. Figure 2D-G:
#'   D  productive clonotype recovery rate per malignant state
#'   E  library complexity per state (detected genes)
#'   F  library complexity per state (counts per cell)
#'   G  TRAC and TRBC1/TRBC2 expression per state
#'
#' The comparison only works on an object built with combineTCR(removeNA = FALSE), because that
#' is what keeps the cells with no recovered contig. On an object where those cells were dropped,
#' the recovery rate is 100% by construction and this script will say so and stop.
#'
#' Dependencies: Seurat, qs2, dplyr, ggplot2, rstatix

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(ggplot2)
})

param <- list(seed = 38, state_col = "Malignant_T",
              tcr_genes = c("TRAC", "TRBC1", "TRBC2"))
set.seed(param$seed)

STATE_COLORS <- c(mTCM = "#4DAF4A", mTRM = "#984EA3", mCyto = "#FF7F00")

# ---------------------------------------------------------------------------------------------
# Recovery
# ---------------------------------------------------------------------------------------------

#' Per-cell productive-clonotype recovery flag
#'
#' A cell counts as recovered if scRepertoire assigned it a strict clonotype. Cells with no
#' contig carry NA in CTstrict, which is exactly why removeNA must have been FALSE upstream.
#'
#' @param obj Seurat object with CTstrict and a malignant state column
#' @return obj with a logical tcr_recovered column
flag_tcr_recovery <- function(obj) {
  md <- obj@meta.data
  stopifnot("CTstrict" %in% colnames(md), param$state_col %in% colnames(md))

  obj$tcr_recovered <- !is.na(md$CTstrict)

  # Positive control on the input, not on the result: if nothing is NA, the cells without a
  # contig were already discarded and any recovery rate computed here is meaningless.
  if (all(obj$tcr_recovered)) {
    stop("Every cell has a clonotype. This object was built with removeNA = TRUE, so ",
         "recovery rate cannot be measured from it. Rebuild with combineTCR(removeNA = FALSE).")
  }
  cat(sprintf("Cells with a productive clonotype: %d / %d (%.1f%%)\n",
              sum(obj$tcr_recovered), ncol(obj), 100 * mean(obj$tcr_recovered)))
  obj
}

#' Recovery rate per malignant state, tested per patient
#'
#' Tested as a patient-level paired comparison rather than a cell-level chi-square: cells within
#' a patient are not independent, and with a clonal malignancy a single patient can otherwise
#' drive the whole result.
#'
#' @param obj Object from flag_tcr_recovery
#' @return list(per_patient, per_state, test)
recovery_by_state <- function(obj) {
  md <- obj@meta.data
  stopifnot("patient" %in% colnames(md))

  per_patient <- md |>
    dplyr::filter(!is.na(.data[[param$state_col]])) |>
    dplyr::group_by(patient, state = .data[[param$state_col]]) |>
    dplyr::summarise(n = dplyr::n(), recovered = sum(tcr_recovered),
                     rate = recovered / n, .groups = "drop") |>
    # A rate over a handful of cells is noise; require a minimum denominator and say how many
    # patient-state combinations that drops.
    dplyr::filter(n >= 10)

  per_state <- per_patient |>
    dplyr::group_by(state) |>
    dplyr::summarise(n_patients = dplyr::n(), median_rate = median(rate),
                     min_rate = min(rate), max_rate = max(rate), .groups = "drop")

  cat("\nRecovery rate per state (patient-level):\n"); print(as.data.frame(per_state))

  test <- NULL
  if (all(c("mTCM", "mTRM") %in% per_patient$state)) {
    w <- per_patient |> dplyr::filter(state %in% c("mTCM", "mTRM")) |>
      dplyr::select(patient, state, rate) |>
      tidyr::pivot_wider(names_from = state, values_from = rate) |>
      dplyr::filter(!is.na(mTCM), !is.na(mTRM))
    if (nrow(w) >= 3) {
      test <- wilcox.test(w$mTCM, w$mTRM, paired = TRUE)
      cat(sprintf("Paired Wilcoxon, mTCM vs mTRM recovery, n = %d patients: p = %.4g\n",
                  nrow(w), test$p.value))
    } else {
      cat("Fewer than 3 patients contribute both states; no paired test attempted.\n")
    }
  }
  list(per_patient = per_patient, per_state = per_state, test = test)
}

# ---------------------------------------------------------------------------------------------
# Complexity and TCR transcript abundance
# ---------------------------------------------------------------------------------------------

#' Library complexity per malignant state
#'
#' If mTRM libraries are simply shallower, the recovery difference is technical. Reporting both
#' detected genes and counts per cell is the point: they can diverge, and only the pair
#' distinguishes "fewer molecules" from "fewer distinct molecules".
complexity_by_state <- function(obj) {
  md <- obj@meta.data |> dplyr::filter(!is.na(.data[[param$state_col]]))
  out <- md |>
    dplyr::group_by(state = .data[[param$state_col]]) |>
    dplyr::summarise(cells = dplyr::n(),
                     median_genes  = median(nFeature_RNA),
                     median_counts = median(nCount_RNA), .groups = "drop")
  cat("\nLibrary complexity per state:\n"); print(as.data.frame(out))
  out
}

#' TRAC and TRBC expression per malignant state
#'
#' Reads from the log-normalised RNA assay. If the TCR constant genes are themselves lower in
#' mTRM, the recovery gap has a transcriptional component and is not purely a contig-assembly
#' artefact.
tcr_gene_expression <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  genes <- intersect(param$tcr_genes, rownames(obj))
  if (length(genes) == 0) stop("None of ", paste(param$tcr_genes, collapse = "/"), " present.")
  if (length(genes) < length(param$tcr_genes)) {
    cat("Note: missing from the object and skipped:",
        paste(setdiff(param$tcr_genes, genes), collapse = ", "), "\n")
  }

  expr <- GetAssayData(obj, layer = "data")[genes, , drop = FALSE]
  df <- data.frame(state = obj@meta.data[[param$state_col]], t(as.matrix(expr)),
                   check.names = FALSE) |>
    dplyr::filter(!is.na(state)) |>
    tidyr::pivot_longer(-state, names_to = "gene", values_to = "expr")

  out <- df |> dplyr::group_by(gene, state) |>
    dplyr::summarise(median_expr = median(expr), pct_detected = 100 * mean(expr > 0),
                     .groups = "drop")
  cat("\nTCR constant-gene expression per state:\n"); print(as.data.frame(out))
  list(long = df, summary = out)
}

# ---------------------------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------------------------

#' @param states Character vector of state levels actually present
#' @return colour vector, interpolated if the object carries states we have no colour for
state_palette <- function(states) {
  if (all(states %in% names(STATE_COLORS))) return(STATE_COLORS[states])
  setNames(colorRampPalette(unname(STATE_COLORS))(length(states)), states)
}

make_plots <- function(rec, cx, tcr, out_dir) {
  pal <- state_palette(sort(unique(rec$per_patient$state)))

  # D: recovery rate, one point per patient, lines connecting the same patient
  pD <- ggplot(rec$per_patient, aes(state, 100 * rate)) +
    geom_boxplot(aes(fill = state), outlier.shape = NA, alpha = 0.6) +
    geom_line(aes(group = patient), colour = "grey60", linewidth = 0.3) +
    geom_point(size = 1.6) +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = "Cells with productive clonotype (%)",
         title = "TCR recovery per malignant state") +
    theme_bw() + theme(legend.position = "none")

  pE <- ggplot(cx, aes(state, median_genes, fill = state)) +
    geom_col() + scale_fill_manual(values = pal) +
    labs(x = NULL, y = "Median genes detected") + theme_bw() +
    theme(legend.position = "none")

  pF <- ggplot(cx, aes(state, median_counts, fill = state)) +
    geom_col() + scale_fill_manual(values = pal) +
    labs(x = NULL, y = "Median counts per cell") + theme_bw() +
    theme(legend.position = "none")

  pG <- ggplot(tcr$long, aes(state, expr, fill = state)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ gene, nrow = 1) +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = "Log-normalised expression") + theme_bw() +
    theme(legend.position = "none")

  for (nm in c("D", "E", "F", "G")) {
    p <- get(paste0("p", nm))
    ggsave(file.path(out_dir, sprintf("SupplFig2%s.png", nm)), p,
           width = if (nm == "G") 9 else 4.5, height = 4, dpi = 300, bg = "white")
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------------------------
# Workflow
# ---------------------------------------------------------------------------------------------

#' @param obj_path Malignant T-cell object carrying CTstrict, patient and the state column,
#'   built with combineTCR(removeNA = FALSE)
#' @param out_dir Output directory
main <- function(obj_path, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  obj <- if (grepl("\\.qs2$", obj_path)) qs2::qs_read(obj_path, nthreads = 1)
         else qs::qread(obj_path, nthreads = 1)

  obj <- flag_tcr_recovery(obj)
  rec <- recovery_by_state(obj)
  cx  <- complexity_by_state(obj)
  tcr <- tcr_gene_expression(obj)

  make_plots(rec, cx, tcr, out_dir)

  write.csv(rec$per_patient, file.path(out_dir, "tcr_recovery_per_patient.csv"), row.names = FALSE)
  write.csv(rec$per_state,   file.path(out_dir, "tcr_recovery_per_state.csv"),   row.names = FALSE)
  write.csv(cx,              file.path(out_dir, "complexity_per_state.csv"),     row.names = FALSE)
  write.csv(tcr$summary,     file.path(out_dir, "tcr_gene_expression.csv"),      row.names = FALSE)
  saveRDS(rec$test, file.path(out_dir, "recovery_paired_test.rds"))
  cat("\nDone.\n")
}

if (sys.nframe() == 0 && !interactive()) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 2) stop("Usage: 21_TCR_recovery_by_malignant_state.R <object> <out_dir>")
  main(a[1], a[2])
}
