# MC result - Tables of Bias, RMSE, Coverage

library(dplyr)
library(tibble)
library(scales)

drop_list_cols <- function(df) df[, !vapply(df, is.list, logical(1)), drop = FALSE]

results_all <- drop_list_cols(readRDS(results_NASC))

results_all$estimator[results_all$estimator == "plain_est"] <- "bsynth_sc"

if (file.exists(results_SCM)) {
  res_scm <- drop_list_cols(readRDS(results_SCM))
  est_in  <- sort(unique(res_scm$estimator))
  
  keep <- if ("sc_classic" %in% est_in) {
    "sc_classic"
  } else if ("plain_est" %in% est_in) {
    "plain_est"
  } else if (length(est_in) == 1L) {
    est_in
  } else {
    stop("Cannot identify the SCM estimator in ", results_SCM,
         " - found: ", paste(est_in, collapse = ", "), ". Set `keep` manually.")
  }
  
  res_scm <- res_scm[res_scm$estimator == keep, , drop = FALSE]
  res_scm$estimator <- "sc_classic"
  
  results_all <- results_all[results_all$estimator != "sc_classic", , drop = FALSE]
  results_all <- dplyr::bind_rows(results_all, res_scm)
}

CONTAM_COL <- NULL
if (is.null(CONTAM_COL)) {
  contam_candidates <- c("contamination", "sc_contamination", "avg_contamination",
                         "contam", "contam_sc", "sc_contam", "w_contamination",
                         "weighted_contamination", "contam_weighted",
                         "contamination_share", "donor_contamination")
  hit <- contam_candidates[contam_candidates %in% names(results_all)]
  CONTAM_COL <- if (length(hit) > 0L) hit[1] else NA_character_
}
if (is.na(CONTAM_COL)) CONTAM_COL <- NULL

COVERAGE_COL     <- NULL
CI_LO_COL        <- NULL
CI_HI_COL        <- NULL
NOMINAL_COVERAGE <- 0.95

.pick_col <- function(cands, cols) { h <- cands[cands %in% cols]; if (length(h)) h[1] else NA_character_ }
.nm <- names(results_all)
if (is.null(COVERAGE_COL))
  COVERAGE_COL <- .pick_col(c("covered", "coverage", "is_covered", "in_ci", "cover", "covers"), .nm)
if (is.null(CI_LO_COL))
  CI_LO_COL <- .pick_col(c("ci_lo", "ci_low", "ci_lower", "att_lo", "att_low", "att_lower",
                           "lower", "conf_lo", "conf_low", "q_lo", "q025", "ci_l"), .nm)
if (is.null(CI_HI_COL))
  CI_HI_COL <- .pick_col(c("ci_hi", "ci_high", "ci_upper", "att_hi", "att_high", "att_upper",
                           "upper", "conf_hi", "conf_high", "q_hi", "q975", "ci_u"), .nm)

if (is.na(COVERAGE_COL)) COVERAGE_COL <- NULL
if (is.na(CI_LO_COL))    CI_LO_COL    <- NULL
if (is.na(CI_HI_COL))    CI_HI_COL    <- NULL

LOG_STEEPNESS <- 50 # log color scale


summarise_perf <- function(df, contam_col = NULL,
                           cov_col = NULL, ci_lo = NULL, ci_hi = NULL) {
  has_c   <- !is.null(contam_col) && contam_col %in% names(df)
  has_cov <- !is.null(cov_col)    && cov_col    %in% names(df)
  has_ci  <- !is.null(ci_lo) && !is.null(ci_hi) &&
    ci_lo %in% names(df) && ci_hi %in% names(df)
  
  df <- df %>% dplyr::filter(status == "ok")
  df$.contam <- if (has_c) df[[contam_col]] else NA_real_
  df$.cov <- if (has_cov)     as.numeric(df[[cov_col]])
  else if (has_ci) as.numeric(df$att_true >= df[[ci_lo]] & df$att_true <= df[[ci_hi]])
  else             NA_real_
  
  df %>%
    dplyr::group_by(ws_k, ws_p, N, T, rho, estimator) %>%
    dplyr::summarise(
      n_ok             = dplyr::n(),
      bias_vs_realized = mean(att_hat - att_true),
      rmse_vs_realized = sqrt(mean((att_hat - att_true)^2)),
      contamination    = if (has_c) mean(.contam, na.rm = TRUE) else NA_real_,
      coverage         = if (has_cov || has_ci) mean(.cov, na.rm = TRUE) else NA_real_,
      mc_sd            = sd(att_hat),
      .groups = "drop"
    )
}


generate_perf_tables <- function(perf,
                                 output_tex       = file.path(dir_table, "nasc_mc_performance_ws_report.tex"),
                                 keep_estimators  = NULL,
                                 estimator_labels = NULL,
                                 t_per_row        = 3L) {
  
  target_cols <- c(bias = "bias_vs_realized", rmse = "rmse_vs_realized")
  
  required_cols <- c(unname(target_cols), "ws_k", "ws_p", "T", "rho", "estimator")
  missing_cols  <- setdiff(required_cols, names(perf))
  if (length(missing_cols) > 0L)
    stop("generate_perf_tables(): missing columns: ", paste(missing_cols, collapse = ", "))
  
  if (!is.null(keep_estimators)) {
    perf <- perf[perf$estimator %in% keep_estimators, , drop = FALSE]
    if (nrow(perf) == 0L)
      stop("generate_perf_tables(): none of keep_estimators found in perf.")
  }
  
  if (!"contamination" %in% names(perf)) perf$contamination <- NA_real_
  if (!"coverage"      %in% names(perf)) perf$coverage      <- NA_real_
  
  perf <- perf %>%
    mutate(
      bias            = .data[[target_cols[["bias"]]]],
      rmse            = .data[[target_cols[["rmse"]]]],
      abs_bias        = abs(.data[[target_cols[["bias"]]]]),
      cover_dev       = abs(coverage - NOMINAL_COVERAGE),
      estimator_label = if (is.null(estimator_labels)) estimator
      else dplyr::coalesce(unname(estimator_labels[estimator]), estimator)
    )
  
  metrics <- tibble::tribble(
    ~key,   ~header, ~print_col, ~color_col, ~digits, ~group,
    "bias", "Bias",  "bias",     "abs_bias", 3L,      1L,
    "rmse", "RMSE",  "rmse",     "rmse",     3L,      1L
  )
  if (any(is.finite(perf$contamination))) {
    metrics <- dplyr::bind_rows(metrics, tibble::tibble(
      key = "contam", header = "Contamination",
      print_col = "contamination", color_col = "contamination", digits = 3L, group = 1L))
  }
  if (any(is.finite(perf$coverage))) {
    metrics <- dplyr::bind_rows(metrics, tibble::tibble(
      key = "cover", header = "Coverage",
      print_col = "coverage", color_col = "cover_dev", digits = 3L, group = 2L))
  }
  
  palette_hex <- c("#74ADD1FF", "#ABD9E9FF", "#E0F3F8FF", "#FFFFBFFF", "#FEE090FF", "#FDAE61FF")
  
  shade_values <- function(vals) {
    out <- rep("FFFFFF", length(vals)); ok <- is.finite(vals)
    if (sum(ok) < 2 || diff(range(vals[ok])) == 0) return(out)
    rank01   <- scales::rescale(vals[ok], to = c(0, 1))
    rank01   <- log1p(LOG_STEEPNESS * rank01) / log1p(LOG_STEEPNESS)
    cols     <- scales::col_numeric(palette_hex, domain = c(0, 1))(rank01)
    out[ok]  <- toupper(sub("^#", "", cols)); out
  }
  
  fmt_num    <- function(x, d) ifelse(is.finite(x), formatC(x, digits = d, format = "f"), "")
  escape_tex <- function(x) gsub("_", "\\\\_", x)
  
  rho_levels <- sort(unique(perf$rho))
  T_levels   <- sort(unique(perf$T))
  tab_cfgs   <- perf %>% distinct(ws_k, ws_p) %>% arrange(ws_k, ws_p)
  
  est_tbl <- perf %>% distinct(estimator, estimator_label)
  if (!is.null(keep_estimators)) {
    est_tbl <- est_tbl %>%
      mutate(.ord = match(estimator, keep_estimators)) %>%
      filter(!is.na(.ord)) %>% arrange(.ord)
  } else {
    est_tbl <- est_tbl %>% arrange(estimator)
  }
  est_levels <- est_tbl$estimator_label
  
  n_T   <- length(T_levels)
  n_rho <- length(rho_levels)
  n_est <- length(est_levels)
  n_m   <- nrow(metrics)
  hrule <- "\\hline"
  
  t_per_row  <- max(1L, as.integer(t_per_row))
  n_col      <- min(t_per_row, n_T)
  T_sections <- split(T_levels, ceiling(seq_along(T_levels) / n_col))
  
  build_table <- function(df_tab, ws_k_v, ws_p_v) {
    
    get_val <- function(Tv, rv, est, col) {
      v <- df_tab[[col]][df_tab$T == Tv & df_tab$rho == rv & df_tab$estimator_label == est]
      if (length(v) == 0) NA_real_ else v
    }
    
    hex_lookup <- list()
    for (m_i in seq_len(n_m)) {
      cc <- metrics$color_col[m_i]
      hx <- rep("FFFFFF", nrow(df_tab))
      for (Tv in unique(df_tab$T)) {
        idx     <- which(df_tab$T == Tv)
        hx[idx] <- shade_values(df_tab[[cc]][idx])
      }
      hex_lookup[[ metrics$key[m_i] ]] <- tibble::tibble(
        T = df_tab$T, rho = df_tab$rho,
        est_label = df_tab$estimator_label, hex = hx
      )
    }
    get_hex <- function(mk, Tv, rv, est) {
      tbl <- hex_lookup[[mk]]
      h   <- tbl$hex[tbl$T == Tv & tbl$rho == rv & tbl$est_label == est]
      if (length(h) == 0) "FFFFFF" else h
    }
    
    est_block   <- paste(sapply(est_levels, escape_tex), collapse = " & ")
    empty_block <- paste(rep("", n_est), collapse = " & ")
    
    build_section <- function(T_cols) {
      n_here <- length(T_cols)
      pad    <- n_col - n_here
      
      header_cells <- c(
        sapply(T_cols, function(Tv)
          sprintf("\\multicolumn{%d}{c}{$T = %d$}", n_est, as.integer(Tv))),
        rep(sprintf("\\multicolumn{%d}{c}{}", n_est), pad)
      )
      top_row <- paste0(" & & ", paste(header_cells, collapse = " & & "), " \\\\")
      
      cline_row <- paste(sapply(seq_len(n_here), function(i) {
        start <- 3L + (i - 1L) * (n_est + 1L)
        end   <- start + n_est - 1L
        sprintf("\\cline{%d-%d}", start, end)
      }), collapse = "")
      
      est_cells <- c(rep(est_block, n_here), rep(empty_block, pad))
      est_row   <- paste0(" & $\\rho$ & ", paste(est_cells, collapse = " & & "), " \\\\")
      
      body <- character(0)
      for (m_i in seq_len(n_m)) {
        mk       <- metrics$key[m_i]
        mlbl     <- metrics$header[m_i]
        pcol     <- metrics$print_col[m_i]
        digits_m <- metrics$digits[m_i]
        
        for (k in seq_along(rho_levels)) {
          rv <- rho_levels[k]
          block_strs <- c(
            sapply(T_cols, function(Tv) {
              blk <- sapply(est_levels, function(est)
                sprintf("\\cellcolor[HTML]{%s}%s",
                        get_hex(mk, Tv, rv, est),
                        fmt_num(get_val(Tv, rv, est, pcol), digits_m)))
              paste(blk, collapse = " & ")
            }),
            rep(empty_block, pad)
          )
          first_col <- if (k == 1L)
            sprintf("\\multirow{%d}{*}{\\rotatebox[origin=c]{90}{%s}}", n_rho, mlbl)
          else ""
          body <- c(body, paste0(first_col, " & ",
                                 formatC(rv, digits = 1, format = "f"), " & ",
                                 paste(block_strs, collapse = " & & "), " \\\\"))
        }
        if (m_i < n_m) body <- c(body, hrule)
      }
      
      c(top_row, cline_row, est_row, hrule, body)
    }
    
    sections <- lapply(T_sections, build_section)
    body_all <- character(0)
    for (s in seq_along(sections)) {
      body_all <- c(body_all, sections[[s]])
      if (s < length(sections)) body_all <- c(body_all, hrule, hrule)
    }
    
    block_spec <- strrep("c", n_est)
    col_spec   <- paste0("ll", block_spec, strrep(paste0("c", block_spec), n_col - 1L))
    p_tag      <- gsub("\\.", "p", formatC(ws_p_v, digits = 2, format = "f"))
    label      <- sprintf("tab:nasc_mc_ws_p%s_k%d_realized", p_tag, ws_k_v)
    
    c(
      "\\begin{sidewaystable}[htbp]",
      "\\centering",
      "\\small",
      sprintf(paste0("\\caption{Title $p = %s$}"),
              formatC(ws_p_v, digits = 2, format = "f")),
      sprintf("\\label{%s}", label),
      "\\resizebox{\\textheight}{!}{%",
      sprintf("\\begin{tabular}{%s}", col_spec),
      hrule,
      hrule,
      body_all,
      hrule,
      hrule,
      "\\end{tabular}%",
      "}",
      "\\end{sidewaystable}",
      ""
    )
  }
  
  all_tables <- unlist(lapply(seq_len(nrow(tab_cfgs)), function(i) {
    build_table(perf %>% filter(ws_k == tab_cfgs$ws_k[i], ws_p == tab_cfgs$ws_p[i]),
                tab_cfgs$ws_k[i], tab_cfgs$ws_p[i])
  }))
  
  writeLines(all_tables, output_tex)
  invisible(output_tex)
}


report_estimators <- c("bc_reg_est", "bc_est", "reg_est", "bsynth_sc", "sc_classic")
report_labels     <- c(bc_reg_est = "NASC",
                       bc_est     = "BC",
                       reg_est    = "CR",
                       bsynth_sc  = "BSCM",
                       sc_classic = "SCM")

perf_report <- summarise_perf(results_all,
                              contam_col = CONTAM_COL,
                              cov_col    = COVERAGE_COL,
                              ci_lo      = CI_LO_COL,
                              ci_hi      = CI_HI_COL)

generate_perf_tables(perf_report,
                     output_tex       = file.path(dir_table, "MC_result_exogen_rho.tex"),
                     keep_estimators  = report_estimators,
                     estimator_labels = report_labels,
                     t_per_row        = length(unique(perf_report$T)))
