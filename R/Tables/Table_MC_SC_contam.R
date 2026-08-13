# MC result - Table of estimated synthetic control contamination

library(dplyr)
library(tibble)
library(scales)

setwd("C:\\Users\\frede\\Desktop\\Master Thesis\\Simulation_R")

RESULTS_RDS     <- "data_ex_rho_2.rds"
RESULTS_RDS_SCM <- "data_ex_rho_scm.rds"
TABLE_DIR       <- "tables"

LOG_STEEPNESS   <- 20 #for log color scale

W_COL   <- "w_hat"
S_COL   <- "s_abs"
IDX_COL <- "donor_idx"

dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

weighted_contam <- function(w, s, idx) {
  if (is.null(w) || is.null(s)) return(NA_real_)
  if (length(s) != length(w) && !is.null(idx) && length(idx) == length(w))
    s <- s[idx]
  if (length(s) != length(w)) return(NA_real_)
  sum(w * s)
}

load_contam <- function(path, w_col = W_COL, s_col = S_COL, idx_col = IDX_COL) {
  d <- readRDS(path)
  stopifnot(w_col %in% names(d), s_col %in% names(d))
  idx <- if (idx_col %in% names(d)) d[[idx_col]] else vector("list", nrow(d))
  d$contam_sc <- mapply(weighted_contam, d[[w_col]], d[[s_col]], idx)
  d[, !vapply(d, is.list, logical(1)), drop = FALSE]
}

res <- load_contam(RESULTS_RDS)
res$estimator[res$estimator == "plain_est"] <- "bsynth_sc"

if (file.exists(RESULTS_RDS_SCM)) {
  res_scm <- load_contam(RESULTS_RDS_SCM)
  est_in  <- sort(unique(res_scm$estimator))

  keep <- if ("sc_classic" %in% est_in) "sc_classic"
  else if ("plain_est" %in% est_in) "plain_est"
  else if (length(est_in) == 1L) est_in
  else stop("Cannot identify the SCM estimator in ", RESULTS_RDS_SCM,
            " - found: ", paste(est_in, collapse = ", "), ". Set `keep` manually.")

  res_scm <- res_scm[res_scm$estimator == keep, , drop = FALSE]
  res_scm$estimator <- "sc_classic"

  res <- res[res$estimator != "sc_classic", , drop = FALSE]
  res <- dplyr::bind_rows(res, res_scm)
}

if (length(unique(res$ws_k)) > 1L) {
  k0 <- sort(unique(res$ws_k))[1]
  res <- res[res$ws_k == k0, , drop = FALSE]
}

res$T <- as.numeric(res$T)
res <- res[res$status == "ok" & is.finite(res$contam_sc), , drop = FALSE]

perf <- res %>%
  dplyr::group_by(ws_p, rho, T, estimator) %>%
  dplyr::summarise(
    contamination = mean(contam_sc, na.rm = TRUE),
    n_ok          = dplyr::n(),
    .groups = "drop"
  )

palette_hex <- c("#74ADD1FF", "#ABD9E9FF", "#E0F3F8FF", "#FFFFBFFF", "#FEE090FF", "#FDAE61FF")

shade_values <- function(vals) {
  out <- rep("FFFFFF", length(vals)); ok <- is.finite(vals)
  if (sum(ok) < 2 || diff(range(vals[ok])) == 0) return(out)
  rank01  <- scales::rescale(vals[ok], to = c(0, 1))
  rank01  <- log1p(LOG_STEEPNESS * rank01) / log1p(LOG_STEEPNESS)
  cols    <- scales::col_numeric(palette_hex, domain = c(0, 1))(rank01)
  out[ok] <- toupper(sub("^#", "", cols)); out
}

fmt_num    <- function(x, d) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "")
escape_tex <- function(x) gsub("_", "\\\\_", x)

report_estimators <- c("bc_reg_est", "bc_est", "reg_est", "bsynth_sc", "sc_classic")
report_labels     <- c(bc_reg_est = "NASC", bc_est = "BC", reg_est = "CR",
                       bsynth_sc  = "BSCM", sc_classic = "SCM")

perf <- perf[perf$estimator %in% report_estimators, , drop = FALSE]
perf$estimator_label <- unname(report_labels[perf$estimator])
est_levels <- unname(report_labels[report_estimators[report_estimators %in% perf$estimator]])


build_contam_table <- function(perf, est_levels,
                               output_tex = file.path(TABLE_DIR, "nasc_estimated_contamination.tex"),
                               digits = 3L) {

  p_levels   <- sort(unique(perf$ws_p))
  rho_levels <- sort(unique(perf$rho))
  T_levels   <- sort(unique(perf$T))
  n_est <- length(est_levels); n_rho <- length(rho_levels); n_T <- length(T_levels)
  hrule <- "\\hline"

  perf$.hex <- rep("FFFFFF", nrow(perf))
  for (Tv in T_levels) {
    idx <- which(perf$T == Tv)
    perf$.hex[idx] <- shade_values(perf$contamination[idx])
  }

  get_val <- function(pv, rv, Tv, est) {
    v <- perf$contamination[perf$ws_p == pv & perf$rho == rv &
                              perf$T == Tv & perf$estimator_label == est]
    if (length(v) == 0) NA_real_ else v[1]
  }
  get_hex <- function(pv, rv, Tv, est) {
    h <- perf$.hex[perf$ws_p == pv & perf$rho == rv &
                     perf$T == Tv & perf$estimator_label == est]
    if (length(h) == 0) "FFFFFF" else h[1]
  }

  est_block <- paste(sapply(est_levels, escape_tex), collapse = " & ")

  header_cells <- sapply(T_levels, function(Tv)
    sprintf("\\multicolumn{%d}{c}{$T = %d$}", n_est, as.integer(Tv)))
  top_row <- paste0(" & & ", paste(header_cells, collapse = " & & "), " \\\\")

  cline_row <- paste(sapply(seq_len(n_T), function(i) {
    start <- 3L + (i - 1L) * (n_est + 1L); end <- start + n_est - 1L
    sprintf("\\cline{%d-%d}", start, end)
  }), collapse = "")

  est_row <- paste0(" & $\\rho$ & ", paste(rep(est_block, n_T), collapse = " & & "), " \\\\")

  body <- character(0)
  for (pi in seq_along(p_levels)) {
    pv <- p_levels[pi]
    for (k in seq_along(rho_levels)) {
      rv <- rho_levels[k]
      cells <- sapply(T_levels, function(Tv) {
        blk <- sapply(est_levels, function(est)
          sprintf("\\cellcolor[HTML]{%s}%s",
                  get_hex(pv, rv, Tv, est),
                  fmt_num(get_val(pv, rv, Tv, est), digits)))
        paste(blk, collapse = " & ")
      })
      first_col <- if (k == 1L)
        sprintf("\\multirow{%d}{*}{\\rotatebox[origin=c]{90}{$p = %s$}}",
                n_rho, formatC(pv, digits = 2, format = "f"))
      else ""
      body <- c(body, paste0(first_col, " & ",
                             formatC(rv, digits = 1, format = "f"), " & ",
                             paste(cells, collapse = " & & "), " \\\\"))
    }
    if (pi < length(p_levels)) body <- c(body, hrule)
  }

  block_spec <- strrep("c", n_est)
  col_spec   <- paste0("ll", block_spec, strrep(paste0("c", block_spec), n_T - 1L))

  tex <- c(
    "\\begin{sidewaystable}[htbp]",
    "\\centering",
    "\\small",
    paste0("\\caption{title}"),
    "\\label{tab:nasc_estimated_contamination}",
    "\\resizebox{\\textheight}{!}{%",
    sprintf("\\begin{tabular}{%s}", col_spec),
    hrule, hrule,
    top_row, cline_row, est_row, hrule,
    body,
    hrule, hrule,
    "\\end{tabular}%",
    "}",
    "\\end{sidewaystable}",
    ""
  )

  writeLines(tex, output_tex)
  invisible(output_tex)
}

build_contam_table(perf, est_levels,
                   output_tex = file.path(TABLE_DIR, "MC_result_exogen_rho_contam.tex"))
