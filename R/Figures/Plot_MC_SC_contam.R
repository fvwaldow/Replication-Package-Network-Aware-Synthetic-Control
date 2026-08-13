# MC results - Box Plots of estimated synthetic control contamination

library(dplyr)

W_COL       <- "w_hat"
S_COL       <- "s_abs"
IDX_COL     <- "donor_idx"
STAR_COL    <- "w_star_donor"

WHISKER_COEF  <- 1.5
SHOW_OUTLIERS <- FALSE
SHOW_ORACLE   <- FALSE

if (!exists("PLOT_FONT")) PLOT_FONT <- "Times"

weighted_contam <- function(w, s, idx) {
  if (is.null(w) || is.null(s)) return(NA_real_)
  if (length(s) != length(w) && !is.null(idx) && length(idx) == length(w))
    s <- s[idx]
  if (length(s) != length(w)) return(NA_real_)
  sum(w * s)
}

load_contam <- function(path, w_col = W_COL, s_col = S_COL,
                        idx_col = IDX_COL, star_col = STAR_COL) {
  d <- readRDS(path)
  stopifnot(w_col %in% names(d), s_col %in% names(d))
  idx <- if (idx_col %in% names(d)) d[[idx_col]] else vector("list", nrow(d))
  d$contam_sc <- mapply(weighted_contam, d[[w_col]], d[[s_col]], idx)
  if (star_col %in% names(d))
    d$contam_oracle <- mapply(weighted_contam, d[[star_col]], d[[s_col]], idx)
  d[, !vapply(d, is.list, logical(1)), drop = FALSE]
}

res <- load_contam(results_NASC)
res$estimator[res$estimator == "plain_est"] <- "bsynth_sc"

if (file.exists(results_SCM)) {
  rs <- load_contam(results_SCM)
  est_in <- sort(unique(rs$estimator))
  keep <- if ("sc_classic" %in% est_in) "sc_classic"
  else if ("plain_est" %in% est_in) "plain_est"
  else if (length(est_in) == 1L) est_in
  else stop(".")
  rs <- rs[rs$estimator == keep, , drop = FALSE]; rs$estimator <- "sc_classic"
  res <- res[res$estimator != "sc_classic", , drop = FALSE]
  res <- dplyr::bind_rows(res, rs)
}

res$T_val <- as.numeric(res$T)

est_levels <- c("sc_classic","bsynth_sc","reg_est","bc_est","bc_reg_est")
est_labels <- c(sc_classic="SCM", bc_est="BC", reg_est="CR",
                bc_reg_est="NASC", bsynth_sc="BSCM")

res <- res[res$status == "ok" & is.finite(res$contam_sc) &
             res$estimator %in% est_levels, , drop = FALSE]

bx <- function(x) boxplot.stats(x, coef = WHISKER_COEF)$stats

agg <- res %>%
  dplyr::group_by(ws_p, estimator, rho, T_val) %>%
  dplyr::summarise(stats = list(bx(contam_sc)), n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    wlo = vapply(stats, `[`, numeric(1), 1),
    q1  = vapply(stats, `[`, numeric(1), 2),
    med = vapply(stats, `[`, numeric(1), 3),
    q3  = vapply(stats, `[`, numeric(1), 4),
    whi = vapply(stats, `[`, numeric(1), 5)
  ) %>%
  dplyr::select(-stats) %>%
  as.data.frame()

has_oracle <- SHOW_ORACLE && "contam_oracle" %in% names(res)
if (has_oracle) {
  agg_or <- res %>%
    dplyr::group_by(ws_p, rho, T_val) %>%
    dplyr::summarise(m = mean(contam_oracle, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
}

est_present <- est_levels[est_levels %in% agg$estimator]
n_est       <- length(est_present)

cols <- c(sc_classic = "#E8442A", bsynth_sc  = "#E6850F",
          bc_reg_est = "#08519C", bc_est     = "#3182BD", reg_est = "#6BAED6")

rho_levels <- sort(unique(agg$rho))

dx      <- 0.06
offsets <- if (n_est > 1) seq(-dx, dx, length.out = n_est) else 0
names(offsets) <- est_present

gap_x  <- if (n_est > 1) 2 * dx / (n_est - 1) else 0.16
box_hw <- 0.38 * gap_x
cap_hw <- 0.55 * box_hw

y_all    <- c(agg$wlo, agg$whi, if (has_oracle) agg_or$m)
c_lim    <- range(y_all, finite = TRUE)
c_lim    <- c_lim + c(-1, 1) * 0.06 * diff(c_lim)
y_at     <- pretty(c_lim, n = 7)

rho_max  <- max(abs(rho_levels))
rho_lim  <- c(-1, 1) * (rho_max + dx + box_hw + 0.08)

smooth_line <- function(x, y, n = 400) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3L) return(list(x = x, y = y))
  spline(x, y, n = n, method = "natural")
}

draw_box <- function(x, wlo, q1, med, q3, whi, col) {
  fill <- adjustcolor(col, alpha.f = 0.35)
  segments(x, q3, x, whi, col = col, lwd = 1.1)
  segments(x, wlo, x, q1, col = col, lwd = 1.1)
  segments(x - cap_hw, whi, x + cap_hw, whi, col = col, lwd = 1.1)
  segments(x - cap_hw, wlo, x + cap_hw, wlo, col = col, lwd = 1.1)
  rect(x - box_hw, q1, x + box_hw, q3, col = fill, border = col, lwd = 1.1)
  segments(x - box_hw, med, x + box_hw, med, col = col, lwd = 1.8)
}

legend_glyph <- function(x, y, col) {
  fill <- adjustcolor(col, alpha.f = 0.35)
  bw <- 0.005; qh <- 0.11; wh <- 0.11; cw <- 0.003
  segments(x, y + qh, x, y + qh + wh, col = col, lwd = 1.1)
  segments(x, y - qh, x, y - qh - wh, col = col, lwd = 1.1)
  segments(x - cw, y + qh + wh, x + cw, y + qh + wh, col = col, lwd = 1.1)
  segments(x - cw, y - qh - wh, x + cw, y - qh - wh, col = col, lwd = 1.1)
  rect(x - bw, y - qh, x + bw, y + qh, col = fill, border = col, lwd = 1.1)
  segments(x - bw, y, x + bw, y, col = col, lwd = 1.8)
  bw
}

draw_panel <- function(p_val, t_val, show_legend = TRUE) {
  d <- agg[agg$ws_p == p_val & agg$T_val == t_val, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  
  if (show_legend) layout(matrix(c(1, 2), nrow = 2), heights = c(6, 0.7))
  par(mar = c(2.8, 3.8, 0.6, 0.8), mgp = c(2, 0.7, 0), family = PLOT_FONT)
  
  plot(NA, xlim = rho_lim, ylim = c_lim,
       xlab = expression(rho),
       ylab = expression(hat(bolditalic(gamma)) * minute * abs(bolditalic(s))),
       xaxt = "n", yaxt = "n")
  
  axis(1, at = rho_levels, labels = formatC(rho_levels, format = "g"))
  
  y_labeled <- c(0.0, 0.1, 0.2, 0.3, 0.4)
  y_unlabeled <- setdiff(y_at, y_labeled)
  
  axis(2, at = y_labeled, labels = y_labeled, las = 1, tcl = -0.5)
  
  if (length(y_unlabeled) > 0) {
    axis(2, at = y_unlabeled, labels = FALSE, tcl = -0.25)
  }
  
  abline(v = rho_levels, h = y_at, col = "grey90", lwd = 0.6)
  
  if (has_oracle) {
    o <- agg_or[agg_or$ws_p == p_val & agg_or$T_val == t_val, ]
    o <- o[order(o$rho), ]
    sm <- smooth_line(o$rho, o$m)
    lines(sm$x, sm$y, col = "grey45", lwd = 1.4, lty = 2)
  }
  
  for (est in est_present) {
    di <- d[d$estimator == est, , drop = FALSE]
    di <- di[order(di$rho), , drop = FALSE]
    if (nrow(di) == 0L) next
    xx <- di$rho + offsets[[est]]
    sm <- smooth_line(xx, di$med)
    lines(sm$x, sm$y, col = adjustcolor(cols[[est]], alpha.f = 0),
          lwd = 1.4, lty = 1)
  }
  
  for (est in est_present) {
    di <- d[d$estimator == est, , drop = FALSE]
    if (nrow(di) == 0L) next
    for (i in seq_len(nrow(di))) {
      x <- di$rho[i] + offsets[[est]]
      draw_box(x, di$wlo[i], di$q1[i], di$med[i], di$q3[i], di$whi[i],
               cols[[est]])
      if (SHOW_OUTLIERS) {
        o <- res$contam_sc[res$ws_p == p_val & res$T_val == t_val &
                             res$estimator == est & res$rho == di$rho[i]]
        o <- boxplot.stats(o, coef = WHISKER_COEF)$out
        if (length(o))
          points(rep(x, length(o)), o, pch = 16, cex = 0.35,
                 col = adjustcolor(cols[[est]], 0.6))
      }
    }
  }
  
  if (show_legend) {
    par(mar = c(0, 0, 0, 0), family = PLOT_FONT)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    leg_order <- c("sc_classic", "bsynth_sc", "reg_est", "bc_est", "bc_reg_est")
    leg_order <- leg_order[leg_order %in% est_present]
    n_leg     <- length(leg_order)
    gbw <- 0.005; g1 <- 0.02; g2 <- 0.06
    lab_w   <- sapply(leg_order, function(e) strwidth(est_labels[[e]], cex = 0.85))
    entry_w <- 2 * gbw + g1 + lab_w
    x <- (1 - (sum(entry_w) + (n_leg - 1) * g2)) / 2
    for (i in seq_len(n_leg)) {
      est <- leg_order[i]
      cx  <- x + gbw
      legend_glyph(cx, 0.5, cols[[est]])
      text(cx + gbw + g1, 0.5, labels = est_labels[[est]],
           adj = c(0, 0.5), cex = 0.85)
      x <- x + entry_w[i] + g2
    }
  }
}

save_panel <- function(p_val, t_val, show_legend = FALSE) {
  if (!any(agg$ws_p == p_val & agg$T_val == t_val)) return(invisible(NULL))
  fname <- file.path(dir_fig,
                     sprintf("contam_sc_box_p%.2f_T%d.pdf", p_val, t_val))
  plot_h <- 3
  dev_h  <- if (show_legend) plot_h * (6 + 0.7) / 6 else plot_h
  pdf(fname, width = 6, height = dev_h, pointsize = 9, family = PLOT_FONT)
  on.exit(dev.off(), add = TRUE)
  draw_panel(p_val, t_val, show_legend = show_legend)
}

save_panel(0,   30, show_legend = TRUE)
save_panel(0.1, 30)
save_panel(0.4, 30, show_legend = TRUE)
