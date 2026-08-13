suppressPackageStartupMessages(library(dplyr))

setwd("C:\\Users\\frede\\Desktop\\Master Thesis\\Simulation_R")

RESULTS_RDS <- "data_ex_rho_2.rds"
DATA_DIR    <- "data_ex_rho_2"
FIG_DIR     <- "figures"
dgp_type    <- "SAR"
treated_idx <- 1L
PLOT_FONT   <- "Times"

WHISKER_COEF  <- 1.5
SHOW_OUTLIERS <- FALSE

BC_EST_LEVELS <- c("bc_est", "bc_reg_est")
est_labels <- c(sc_classic = "SCM", bc_est = "BC", reg_est = "CR",
                bc_reg_est = "NASC", bsynth_sc = "BSCM")
cols <- c(sc_classic = "#E8442A", bsynth_sc = "#E6850F",
          bc_reg_est = "#08519C", bc_est = "#3182BD", reg_est = "#6BAED6")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

.net_tag <- function(dgp_type, ws_k, ws_p, N)
  sprintf("%s_k%d_p%.2f_N%d", dgp_type, ws_k, ws_p, N)

W_cache <- new.env(parent = emptyenv())
get_W <- function(ws_k, ws_p, N) {
  key <- .net_tag(dgp_type, ws_k, ws_p, N)
  if (!is.null(W_cache[[key]])) return(W_cache[[key]])
  f <- file.path(DATA_DIR, key, "network_W.rds")
  W <- tryCatch(as.matrix(readRDS(f)), error = function(e) NULL)
  W_cache[[key]] <- W
  W
}

signed_s <- function(W, rho, donor_idx, treated_idx) {
  J    <- length(donor_idx)
  W_J  <- W[donor_idx, donor_idx, drop = FALSE]
  w_J1 <- W[donor_idx, treated_idx]
  as.numeric(rho * solve(diag(J) - rho * W_J, w_J1))
}

d <- readRDS(RESULTS_RDS)
d$T_val <- suppressWarnings(as.numeric(d$T))

n      <- nrow(d)
bc_vec <- rep(NA_real_, n)
for (i in seq_len(n)) {
  W_i <- get_W(d$ws_k[i], d$ws_p[i], d$N[i])
  if (is.null(W_i)) next
  s_i <- tryCatch(signed_s(W_i, d$rho[i], d$donor_idx[[i]], treated_idx),
                  error = function(e) NULL)
  if (is.null(s_i)) next
  w_i <- suppressWarnings(as.numeric(d$w_hat[[i]]))
  if (length(w_i) == length(s_i) && !anyNA(w_i))
    bc_vec[i] <- 1 / (1 - sum(w_i * s_i))
}
d$bc_factor <- bc_vec

keep <- !is.na(d$status) & d$status == "ok" & is.finite(d$bc_factor) &
  d$estimator %in% BC_EST_LEVELS
d <- d[keep, , drop = FALSE]

bx <- function(x) boxplot.stats(x, coef = WHISKER_COEF)$stats

agg <- d %>%
  dplyr::group_by(ws_p, estimator, rho, T_val) %>%
  dplyr::summarise(stats = list(bx(bc_factor)), .groups = "drop") %>%
  dplyr::mutate(
    wlo = vapply(stats, `[`, numeric(1), 1),
    q1  = vapply(stats, `[`, numeric(1), 2),
    med = vapply(stats, `[`, numeric(1), 3),
    q3  = vapply(stats, `[`, numeric(1), 4),
    whi = vapply(stats, `[`, numeric(1), 5)
  ) %>%
  dplyr::select(-stats) %>%
  as.data.frame()

est_present <- BC_EST_LEVELS[BC_EST_LEVELS %in% agg$estimator]
n_est       <- length(est_present)

n_ref  <- 5
dx     <- 0.06
gap_x  <- 4 * dx / (n_ref - 1)
box_hw <- 0.38 * gap_x
cap_hw <- 0.55 * box_hw
offsets <- (seq_len(n_est) - (n_est + 1) / 2) * gap_x
names(offsets) <- est_present

rho_levels <- sort(unique(agg$rho))
off_max <- if (length(offsets)) max(abs(offsets)) else 0
rho_lim <- c(-1, 1) * (max(abs(rho_levels)) + off_max + box_hw + 0.08)

c_lim <- c(0.9, 1.45)
y_major <- round(seq(ceiling(c_lim[1] / 0.1) * 0.1,
                     floor(c_lim[2] / 0.1) * 0.1, by = 0.1), 1)

smooth_line <- function(x, y, n = 400) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 3L) return(list(x = x, y = y))
  stats::spline(x, y, n = n, method = "natural")
}

draw_box <- function(x, wlo, q1, med, q3, whi, col) {
  fill <- grDevices::adjustcolor(col, alpha.f = 0.35)
  graphics::segments(x, q3, x, whi, col = col, lwd = 1.1)
  graphics::segments(x, wlo, x, q1, col = col, lwd = 1.1)
  graphics::segments(x - cap_hw, whi, x + cap_hw, whi, col = col, lwd = 1.1)
  graphics::segments(x - cap_hw, wlo, x + cap_hw, wlo, col = col, lwd = 1.1)
  graphics::rect(x - box_hw, q1, x + box_hw, q3, col = fill, border = col, lwd = 1.1)
  graphics::segments(x - box_hw, med, x + box_hw, med, col = col, lwd = 1.8)
}

legend_glyph <- function(x, y, col) {
  fill <- grDevices::adjustcolor(col, alpha.f = 0.35)
  bw <- 0.005; qh <- 0.11; wh <- 0.11; cw <- 0.003
  graphics::segments(x, y + qh, x, y + qh + wh, col = col, lwd = 1.1)
  graphics::segments(x, y - qh, x, y - qh - wh, col = col, lwd = 1.1)
  graphics::segments(x - cw, y + qh + wh, x + cw, y + qh + wh, col = col, lwd = 1.1)
  graphics::segments(x - cw, y - qh - wh, x + cw, y - qh - wh, col = col, lwd = 1.1)
  graphics::rect(x - bw, y - qh, x + bw, y + qh, col = fill, border = col, lwd = 1.1)
  graphics::segments(x - bw, y, x + bw, y, col = col, lwd = 1.8)
  bw
}

draw_panel <- function(p_val, t_val, show_legend = TRUE) {
  dd <- agg[agg$ws_p == p_val & agg$T_val == t_val, , drop = FALSE]
  if (nrow(dd) == 0) return(invisible(NULL))

  if (show_legend)
    graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(6, 0.7))
  graphics::par(mar = c(2.8, 3.8, 0.6, 0.8), mgp = c(2, 0.7, 0),
                family = PLOT_FONT)

  plot(NA, xlim = rho_lim, ylim = c_lim,
       xlab = expression(rho),
       ylab = expression((1 - hat(bolditalic(gamma)) * minute * bolditalic(s))^-1),
       xaxt = "n", yaxt = "n")
  graphics::axis(1, at = rho_levels, labels = formatC(rho_levels, format = "g"))
  graphics::axis(2, at = y_major, labels = formatC(y_major, format = "f", digits = 1),
                 las = 1, tcl = -0.5)
  graphics::abline(v = rho_levels, h = y_major, col = "grey90", lwd = 0.6)
  graphics::abline(h = 1, col = "grey55", lwd = 0.9, lty = 2)

  for (est in est_present) {
    di <- dd[dd$estimator == est, , drop = FALSE]
    di <- di[order(di$rho), , drop = FALSE]
    if (nrow(di) == 0L) next
    sm <- smooth_line(di$rho + offsets[[est]], di$med)
    graphics::lines(sm$x, sm$y,
                    col = grDevices::adjustcolor(cols[[est]], alpha.f = 0),
                    lwd = 1.4)
  }

  for (est in est_present) {
    di <- dd[dd$estimator == est, , drop = FALSE]
    if (nrow(di) == 0L) next
    for (i in seq_len(nrow(di))) {
      x <- di$rho[i] + offsets[[est]]
      draw_box(x, di$wlo[i], di$q1[i], di$med[i], di$q3[i], di$whi[i],
               cols[[est]])
      if (SHOW_OUTLIERS) {
        oo <- d$bc_factor[d$ws_p == p_val & d$T_val == t_val &
                            d$estimator == est & d$rho == di$rho[i]]
        oo <- boxplot.stats(oo, coef = WHISKER_COEF)$out
        if (length(oo))
          graphics::points(rep(x, length(oo)), oo, pch = 16, cex = 0.35,
                           col = grDevices::adjustcolor(cols[[est]], 0.6))
      }
    }
  }
  graphics::box()

  if (show_legend) {
    graphics::par(mar = c(0, 0, 0, 0), family = PLOT_FONT)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    leg_order <- BC_EST_LEVELS[BC_EST_LEVELS %in% est_present]
    n_leg <- length(leg_order)
    gbw <- 0.005; g1 <- 0.02; g2 <- 0.06
    lab_w <- vapply(leg_order,
                    function(e) graphics::strwidth(est_labels[[e]], cex = 0.85),
                    numeric(1))
    entry_w <- 2 * gbw + g1 + lab_w
    x <- (1 - (sum(entry_w) + (n_leg - 1) * g2)) / 2
    for (i in seq_len(n_leg)) {
      est <- leg_order[i]
      cx  <- x + gbw
      legend_glyph(cx, 0.5, cols[[est]])
      graphics::text(cx + gbw + g1, 0.5, labels = est_labels[[est]],
                     adj = c(0, 0.5), cex = 0.85)
      x <- x + entry_w[i] + g2
    }
  }
  invisible(NULL)
}

save_panel <- function(p_val, t_val, show_legend = FALSE) {
  if (!any(agg$ws_p == p_val & agg$T_val == t_val)) return(invisible(NULL))
  fname <- file.path(FIG_DIR,
                     sprintf("bc_factor_box_p%.2f_T%d.pdf", p_val, t_val))
  plot_h <- 2
  dev_h  <- if (show_legend) plot_h * (6 + 0.7) / 6 else plot_h
  grDevices::pdf(fname, width = 6, height = dev_h, pointsize = 9,
                 family = PLOT_FONT)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_panel(p_val, t_val, show_legend = show_legend)
}

invisible(save_panel(0,   30))
invisible(save_panel(0.1, 30))
invisible(save_panel(0.4, 30, show_legend = TRUE))
