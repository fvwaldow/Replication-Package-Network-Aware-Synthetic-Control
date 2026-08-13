# Single MC iteration - SCM, BSCM, CR, BC, NASC

library(nasc)
library(dplyr)
library(igraph)
library(Synth)

N        <- 15
T_0      <- 30
T        <- T_0 + 10
rho_true <- 0.6
ws_p     <- 0.0

dgp_type       <- "SAR"
treated_idx    <- 1
beta           <- c(1.0, 0.5)
theta          <- c(0.3, 0.2)
delta_mean     <- 5
delta_sd       <- 0
sigma_u        <- 0.3
alpha_sd       <- 0.3
x_sd           <- 0.3
X_mean         <- c(0.0, 0.0)
twin_target    <- "cleanest"
weight_profile <- c(0.075, 0.075, 0.150, 0.150, 0.250, 0.300)

lambda_cv_grid    <- c(0, exp(seq(log(0.05), log(50), length.out = 20)))
lambda_train_frac <- 0.8

stan_iter    <- 1000L
stan_warmup  <- 500L
stan_control <- list(adapt_delta = 0.95, max_treedepth = 10)

if (!exists("PLOT_FONT")) PLOT_FONT <- "Times"
PANEL_MAR   <- c(2.8, 3.8, 0.6, 0.8)
PANEL_MGP   <- c(2, 0.7, 0)
LEGEND_H    <- 0.7
LEGEND_CEX  <- 0.85

LW_OBS   <- 1.8
LW_PATH  <- 1.6
LW_DENS  <- 1.1
LW_TICK  <- 1.8
LW_GRID  <- 0.6
LW_RULE  <- 0.8

nasc_palette <- c(
  sc_classic = "#E8442A",
  bsynth_sc  = "#E6850F",
  bc_reg_est = "#08519C",
  bc_est     = "#3182BD",
  reg_est    = "#6BAED6"
)

# Plotting functions

.nasc_model_colors <- function(model_names, colors = NULL) {
  out <- vapply(model_names, function(nm) {
    if (!is.null(colors) && nm %in% names(colors)) return(unname(colors[nm]))
    n <- tolower(nm)
    if (grepl("full", n) ||
        (grepl("nasc", n) && !grepl("pen|cr|reg", n)) ||
        (grepl("bc|bias", n) && grepl("cr|pen|reg", n))) {
      nasc_palette[["bc_reg_est"]]
    } else if (grepl("bsynth|bscm|bayes", n)) {
      nasc_palette[["bsynth_sc"]]
    } else if (grepl("cr|pen|reg", n)) {
      nasc_palette[["reg_est"]]
    } else if (grepl("bc|bias", n)) {
      nasc_palette[["bc_est"]]
    } else if (grepl("synth|scm|classic|sc", n)) {
      nasc_palette[["sc_classic"]]
    } else {
      NA_character_
    }
  }, character(1))
  if (any(is.na(out))) {
    fb <- grDevices::hcl.colors(max(sum(is.na(out)), 2), palette = "Dark 3")
    out[is.na(out)] <- fb[seq_len(sum(is.na(out)))]
  }
  unname(out)
}

.resolve_nasc_internal <- function(nm) {
  if (exists(nm, mode = "function")) return(get(nm, mode = "function"))
  tryCatch(utils::getFromNamespace(nm, "nasc"), error = function(e) NULL)
}
.resolve_indirect_matrix <- function() {
  .resolve_nasc_internal(".nasc_indirect_matrix")
}

fit_synth_unit <- function(df, treated_id, donor_ids, intervention_time,
                           use_covariates = TRUE) {
  
  dfS <- as.data.frame(df)
  dfS$id <- as.character(dfS$id)
  
  key <- unique(dfS$id)
  num_try <- suppressWarnings(as.numeric(key))
  key <- if (!any(is.na(num_try))) key[order(num_try)] else sort(key)
  
  dfS$..unit_num  <- match(dfS$id, key)
  dfS$..unit_name <- dfS$id
  
  all_times <- sort(unique(dfS$time))
  pre_times <- all_times[all_times < intervention_time]
  
  sp    <- lapply(pre_times, function(t) list("Y", t, "mean"))
  preds <- if (use_covariates) c("X1", "X2") else NULL
  n_pred <- length(sp) + length(preds)
  
  suppressMessages(invisible(utils::capture.output({
    dp <- Synth::dataprep(
      foo                   = dfS,
      predictors            = preds,
      predictors.op         = "mean",
      special.predictors    = sp,
      dependent             = "Y",
      unit.variable         = "..unit_num",
      unit.names.variable   = "..unit_name",
      time.variable         = "time",
      treatment.identifier  = match(as.character(treated_id), key),
      controls.identifier   = match(as.character(donor_ids), key),
      time.predictors.prior = pre_times,
      time.optimize.ssr     = pre_times,
      time.plot             = all_times
    )
    out <- Synth::synth(data.prep.obj = dp,
                        custom.v = rep(1 / n_pred, n_pred))
  })))
  
  w <- as.numeric(out$solution.w)
  names(w) <- key[as.integer(rownames(out$solution.w))]
  
  y_obs   <- as.numeric(dp$Y1plot)
  y_synth <- as.numeric(dp$Y0plot %*% out$solution.w)
  
  list(
    unit    = as.character(treated_id),
    times   = all_times,
    y_obs   = y_obs,
    y_synth = y_synth,
    gap     = y_obs - y_synth,
    w       = w
  )
}

.rmspe <- function(x) sqrt(mean(x^2))

synth_permutation_inference <- function(fit_tr, placebo_fits,
                                        intervention_time,
                                        ci_width = 0.95) {
  
  keep <- !vapply(placebo_fits, is.null, logical(1))
  placebo_fits <- placebo_fits[keep]
  
  times <- fit_tr$times
  pre   <- times < intervention_time
  
  r1_pre  <- .rmspe(fit_tr$gap[pre])
  r1_post <- .rmspe(fit_tr$gap[!pre])
  ratio1  <- r1_post / r1_pre
  
  G <- vapply(placebo_fits, `[[`, numeric(length(times)), "gap")
  rj_pre  <- apply(G[pre, , drop = FALSE],  2, .rmspe)
  rj_post <- apply(G[!pre, , drop = FALSE], 2, .rmspe)
  ratios  <- rj_post / rj_pre
  names(ratios) <- vapply(placebo_fits, `[[`, character(1), "unit")
  
  p_value <- mean(c(ratios, ratio1) >= ratio1)
  
  Gs   <- sweep(G, 2, rj_pre, "/")
  a    <- (1 - ci_width) / 2
  q_lo <- apply(Gs, 1, stats::quantile, probs = a,     names = FALSE)
  q_hi <- apply(Gs, 1, stats::quantile, probs = 1 - a, names = FALSE)
  
  list(
    p_value        = p_value,
    ratio_treated  = ratio1,
    ratios_placebo = sort(ratios, decreasing = TRUE),
    rmspe_pre      = r1_pre,
    rmspe_post     = r1_post,
    tau_LB         = fit_tr$gap - q_hi * r1_pre,
    tau_UB         = fit_tr$gap - q_lo * r1_pre,
    placebo_gaps   = G,
    ci_width       = ci_width
  )
}

synthSCM <- function(df, treated_id, intervention_time,
                     ci_width = 0.95, use_covariates = TRUE) {
  
  ids <- as.character(unique(df$id))
  num_try <- suppressWarnings(as.numeric(ids))
  ids <- if (!any(is.na(num_try))) ids[order(num_try)] else sort(ids)
  donor_ids <- setdiff(ids, as.character(treated_id))
  
  fit_tr <- fit_synth_unit(df, treated_id, donor_ids,
                           intervention_time, use_covariates)
  
  placebo_fits <- lapply(donor_ids, function(j) {
    tryCatch(
      fit_synth_unit(df, j, setdiff(donor_ids, j),
                     intervention_time, use_covariates),
      error = function(e) NULL
    )
  })
  
  inf <- synth_permutation_inference(fit_tr, placebo_fits,
                                     intervention_time, ci_width)
  
  pd <- data.frame(
    time    = fit_tr$times,
    Y       = fit_tr$y_obs,
    y_synth = fit_tr$y_synth,
    tau     = fit_tr$gap,
    tau_LB  = inf$tau_LB,
    tau_UB  = inf$tau_UB
  )
  pd$LB <- pd$Y - pd$tau_UB
  pd$UB <- pd$Y - pd$tau_LB
  pd <- pd[, c("time", "Y", "y_synth", "LB", "UB",
               "tau", "tau_LB", "tau_UB")]
  
  obj <- list(
    plotData         = pd,
    interventionTime = intervention_time,
    weights          = fit_tr$w,
    fit              = fit_tr,
    placebos         = placebo_fits,
    inference        = inf
  )
  class(obj) <- c("synthSCM", "list")
  obj
}

nascWeightSCM <- function(models,
                          synth_w,
                          true_w      = NULL,
                          synth_label = "SCM",
                          true_label  = "True Weights",
                          colors      = NULL,
                          font_family = PLOT_FONT,
                          show_legend = TRUE,
                          overlap     = 0.3,
                          scale       = 2,
                          fill_alpha  = 0.2,
                          max_donors  = NULL) {
  
  per_model <- lapply(names(models), function(m_name) {
    mod  <- models[[m_name]]
    priv <- mod$.__enclos_env__$private
    w_mat <- priv$y_synth_draws$w
    donor_ids <- priv$donor_ids
    colnames(w_mat) <- donor_ids
    list(name = m_name, w_mat = w_mat, donors = donor_ids)
  })
  
  all_donors <- unique(c(unlist(lapply(per_model, `[[`, "donors")),
                         names(synth_w), names(true_w)))
  num_attempt <- suppressWarnings(as.numeric(all_donors))
  donor_order <- if (!any(is.na(num_attempt))) {
    all_donors[order(num_attempt)]
  } else {
    all_donors[order(as.character(all_donors))]
  }
  if (!is.null(max_donors) && length(donor_order) > max_donors) {
    keep_n <- as.integer(max_donors)
    donor_order <- donor_order[(length(donor_order) - keep_n + 1L):
                                 length(donor_order)]
  }
  donor_plot <- donor_order
  n_d <- length(donor_plot)
  n_m <- length(per_model)
  
  dens_grid <- vector("list", n_d)
  for (i in seq_len(n_d)) {
    donor <- donor_plot[i]
    dens_grid[[i]] <- lapply(per_model, function(pm) {
      if (donor %in% pm$donors) {
        stats::density(pm$w_mat[, donor], na.rm = TRUE)
      } else NULL
    })
  }
  
  all_x <- unlist(lapply(dens_grid, function(row)
    unlist(lapply(row, function(d) if (!is.null(d)) d$x else NULL))))
  x_range <- range(c(all_x, synth_w, true_w, 0), na.rm = TRUE)
  max_y <- max(vapply(dens_grid, function(row) {
    ys <- vapply(row, function(d) if (!is.null(d)) max(d$y) else 0,
                 numeric(1))
    if (length(ys)) max(ys) else 0
  }, numeric(1)))
  if (!is.finite(max_y) || max_y <= 0) max_y <- 1
  
  ridge_h  <- scale * max_y
  step     <- ridge_h * (1 - overlap)
  ylim_top <- step * n_d + ridge_h *0.8
  ylim_bot <- 20
  
  cols  <- .nasc_model_colors(names(models), colors)
  fills <- grDevices::adjustcolor(cols, alpha.f = fill_alpha)
  synth_col <- "#E8442A"
  true_col  <- "gray40"
  
  op <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(op); graphics::layout(1) })
  
  if (isTRUE(show_legend)) {
    graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(6, LEGEND_H))
  }
  graphics::par(mar = PANEL_MAR, mgp = PANEL_MGP, family = font_family)
  
  plot(NA, xlim = x_range, ylim = c(ylim_bot, ylim_top),
       xlab = expression(hat(gamma[j])), ylab = "Donor Unit j", yaxt = "n")
  graphics::axis(2, at = step * seq_len(n_d), labels = donor_plot, las = 1)
  
  for (i in seq(n_d, 1L, by = -1L)) {
    baseline <- step * i
    graphics::segments(x_range[1], baseline, x_range[2], baseline,
                       col = "gray60", lwd = LW_RULE)
    for (m in seq_len(n_m)) {
      d <- dens_grid[[i]][[m]]
      if (is.null(d)) next
      y <- baseline + d$y * (ridge_h / max_y)
      graphics::polygon(x = c(d$x, rev(d$x)),
                        y = c(y, rep(baseline, length(d$x))),
                        col = fills[m], border = cols[m], lwd = LW_DENS)
    }
    donor <- donor_plot[i]
    if (!is.null(true_w) && donor %in% names(true_w)) {
      graphics::segments(true_w[donor], baseline,
                         true_w[donor], baseline + 0.55 * ridge_h,
                         col = grDevices::adjustcolor(true_col, alpha.f = 0.9),
                         lwd = LW_TICK, lty = 1)
    }
    if (donor %in% names(synth_w)) {
      graphics::segments(synth_w[donor], baseline,
                         synth_w[donor], baseline + 0.55 * ridge_h,
                         col = grDevices::adjustcolor(synth_col, alpha.f = 0.9),
                         lwd = LW_TICK, lty = 1)
    }
  }
  graphics::box()
  
  extra_lab <- synth_label
  extra_lty <- 1
  extra_col <- synth_col
  if (!is.null(true_w)) {
    extra_lab <- c(extra_lab, true_label)
    extra_lty <- c(extra_lty, 1)
    extra_col <- c(extra_col, true_col)
  }
  n_e <- length(extra_lab)
  
  if (isTRUE(show_legend)) {
    graphics::par(mar = c(0, 0, 0, 0), family = font_family)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    graphics::legend("center",
                     legend  = c(names(models), extra_lab),
                     fill    = c(fills, rep(NA, n_e)),
                     border  = c(cols, rep(NA, n_e)),
                     lty     = c(rep(NA, n_m), extra_lty),
                     lwd     = c(rep(NA, n_m), rep(LW_TICK, n_e)),
                     col     = c(rep(NA, n_m), extra_col),
                     horiz   = TRUE, bty = "n", cex = LEGEND_CEX)
  }
  
  rows <- list()
  for (donor in rev(donor_order)) {
    for (pm in per_model) {
      if (donor %in% pm$donors) {
        col <- pm$w_mat[, donor]
        rows[[length(rows) + 1L]] <- data.frame(
          donor = donor, model = pm$name,
          mean = mean(col, na.rm = TRUE),
          sd   = stats::sd(col, na.rm = TRUE),
          share_gt0 = mean(col > 1e-3, na.rm = TRUE))
      }
    }
    if (donor %in% names(synth_w)) {
      rows[[length(rows) + 1L]] <- data.frame(
        donor = donor, model = synth_label,
        mean = unname(synth_w[donor]), sd = NA_real_,
        share_gt0 = as.numeric(synth_w[donor] > 1e-3))
    }
    if (!is.null(true_w) && donor %in% names(true_w)) {
      rows[[length(rows) + 1L]] <- data.frame(
        donor = donor, model = true_label,
        mean = unname(true_w[donor]), sd = NA_real_,
        share_gt0 = as.numeric(true_w[donor] > 1e-3))
    }
  }
  invisible(do.call(rbind, rows))
}

nascPlot <- function(models, show_ci = FALSE,
                     ci_models       = c("SCM", "NASC"),
                     colors          = NULL,
                     font_family     = PLOT_FONT,
                     panels          = c("path", "effect"),
                     show_legend     = TRUE,
                     indirect_models = c("BC", "NASC"),
                     indirect_ylim   = NULL,
                     indirect_pre    = c("zero", "drop", "placebo"),
                     indirect_legend = c("drawn", "all"),
                     true_att        = NULL,
                     true_indirect   = NULL) {
  
  indirect_pre    <- match.arg(indirect_pre)
  indirect_legend <- match.arg(indirect_legend)
  
  panels <- match.arg(panels, c("path", "effect", "indirect"),
                      several.ok = TRUE)
  
  mod1 <- models[[1]]
  
  intervention_time <- mod1$interventionTime
  
  combined_list <- lapply(names(models), function(m_name) {
    df <- models[[m_name]]$plotData
    df <- as.data.frame(df)
    colnames(df)[1:2] <- c("time_var", "outcome_var")
    df$Model <- m_name
    df <- df[order(df$time_var), , drop = FALSE]
    df
  })
  names(combined_list) <- names(models)
  
  n_models <- length(models)
  cols <- .nasc_model_colors(names(models), colors)
  
  ci_draw <- if (is.null(ci_models)) rep(TRUE, n_models) else
    names(models) %in% ci_models
  
  op <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(op); graphics::layout(1) }, add = TRUE)
  
  obs_df <- combined_list[[1]][, c("time_var", "outcome_var")]
  xrng   <- range(obs_df$time_var, na.rm = TRUE)
  
  if ("path" %in% panels) {
    
    y_vals <- c(obs_df$outcome_var,
                unlist(lapply(combined_list, function(d) d$y_synth)))
    if (isTRUE(show_ci)) {
      y_vals <- c(y_vals,
                  unlist(lapply(combined_list[ci_draw],
                                function(d) c(d$LB, d$UB))))
    }
    yrng <- range(y_vals, na.rm = TRUE)
    
    if (isTRUE(show_legend)) {
      graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(6, LEGEND_H))
    }
    graphics::par(mar = PANEL_MAR, mgp = PANEL_MGP, family = font_family)
    
    plot(obs_df$time_var, obs_df$outcome_var, type = "n",
         xlim = xrng, ylim = yrng,
         xlab = "t", ylab = "y")
    graphics::grid(lty = 1, col = "grey90", lwd = LW_GRID)
    
    if (isTRUE(show_ci)) {
      for (i in seq_len(n_models)) {
        if (!ci_draw[i]) next
        d <- combined_list[[i]]
        graphics::polygon(c(d$time_var, rev(d$time_var)),
                          c(d$LB, rev(d$UB)),
                          col = grDevices::adjustcolor(cols[i], alpha.f = 0.15),
                          border = NA)
      }
    }
    
    graphics::abline(v = intervention_time, lty = 2, col = "gray30",
                     lwd = LW_RULE)
    graphics::lines(obs_df$time_var, obs_df$outcome_var,
                    lwd = LW_OBS, col = "black", lty = 1)
    for (i in seq_len(n_models)) {
      d <- combined_list[[i]]
      graphics::lines(d$time_var, d$y_synth,
                      col = cols[i], lwd = LW_PATH, lty = 1)
    }
    graphics::box()
    
    if (isTRUE(show_legend)) {
      graphics::par(mar = c(0, 0, 0, 0), family = font_family)
      graphics::plot.new()
      graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
      graphics::legend("center",
                       legend = c("Observed", names(models)),
                       col    = c("black", cols),
                       lty    = 1,
                       lwd    = c(LW_OBS, rep(LW_PATH, n_models)),
                       horiz  = TRUE, bty = "n", cex = LEGEND_CEX)
    }
  }
  
  if ("effect" %in% panels) {
    y_vals2 <- unlist(lapply(combined_list, function(d) d$tau))
    if (isTRUE(show_ci)) {
      y_vals2 <- c(y_vals2,
                   unlist(lapply(combined_list[ci_draw],
                                 function(d) c(d$tau_LB, d$tau_UB))))
    }
    
    effect_y_by <- 1
    yrng2 <- c(-1.5, 6.5)
    y_at2 <- seq(ceiling(yrng2[1] / effect_y_by) * effect_y_by,
                 floor(yrng2[2]   / effect_y_by) * effect_y_by,
                 by = effect_y_by)
    
    if (isTRUE(show_legend)) {
      graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(6, LEGEND_H))
    }
    graphics::par(mar = PANEL_MAR, mgp = PANEL_MGP, family = font_family)
    
    plot(combined_list[[1]]$time_var, combined_list[[1]]$tau, type = "n",
         xlim = xrng, ylim = yrng2,
         xlab = "t", ylab = expression(hat(tau)),
         yaxt = "n")
    graphics::axis(2, at = y_at2, las = 1)
    graphics::grid(nx = NULL, ny = NA, lty = 1, col = "grey90", lwd = LW_GRID)
    graphics::abline(h = y_at2, col = "grey90", lwd = LW_GRID)
    if (isTRUE(show_ci)) {
      for (i in seq_len(n_models)) {
        if (!ci_draw[i]) next
        d <- combined_list[[i]]
        graphics::polygon(c(d$time_var, rev(d$time_var)),
                          c(d$tau_LB, rev(d$tau_UB)),
                          col = grDevices::adjustcolor(cols[i], alpha.f = 0.15),
                          border = NA)
      }
    }
    
    graphics::abline(h = 0, lty = 1, col = "grey70", lwd = LW_RULE)
    graphics::abline(v = intervention_time, lty = 2, col = "gray30",
                     lwd = LW_RULE)
    if (!is.null(true_att)) {
      graphics::segments(xrng[1], true_att, xrng[2], true_att,
                         col = "gray30", lty = 2, lwd = LW_RULE)
    }
    for (i in seq_len(n_models)) {
      d <- combined_list[[i]]
      graphics::lines(d$time_var, d$tau, col = cols[i], lwd = LW_PATH)
    }
    graphics::box()
    
    if (isTRUE(show_legend)) {
      graphics::par(mar = c(0, 0, 0, 0), family = font_family)
      graphics::plot.new()
      graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
      graphics::legend("center",
                       legend = names(models),
                       col    = cols,
                       lty    = 1, lwd = LW_PATH,
                       horiz  = TRUE, bty = "n", cex = LEGEND_CEX)
    }
  }
  
  if ("indirect" %in% panels) {
    get_ind <- .resolve_indirect_matrix()
    
    ind_draw <- if (is.null(indirect_models)) rep(TRUE, n_models) else
      names(models) %in% indirect_models
    
    ind_list <- vector("list", n_models)
    for (i in seq_len(n_models)) {
      if (!ind_draw[i]) next
      m <- models[[i]]
      if (!inherits(m, "nascSynth")) next
      ind_list[i] <- list(get_ind(m, pre = indirect_pre))
    }
    has_ind <- !vapply(ind_list, is.null, logical(1))
    
    if (any(has_ind)) {
      
      yrng3 <- if (!is.null(indirect_ylim)) {
        indirect_ylim
      } else {
        range(c(0, unlist(lapply(ind_list[has_ind],
                                 function(d) as.numeric(d$mean))),
                if (!is.null(true_indirect)) as.numeric(true_indirect)),
              na.rm = TRUE)
      }
      
      if (isTRUE(show_legend)) {
        graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(6, LEGEND_H))
      }
      graphics::par(mar = PANEL_MAR, mgp = PANEL_MGP, family = font_family)
      
      plot(NA, xlim = xrng, ylim = yrng3,
           xlab = "t", ylab = expression(hat(delta)[j]))
      graphics::grid(lty = 1, col = "grey90", lwd = LW_GRID)
      graphics::abline(h = 0, lty = 1, col = "grey70", lwd = LW_RULE)
      graphics::abline(v = intervention_time, lty = 2, col = "gray30",
                       lwd = LW_RULE)
      if (!is.null(true_indirect)) {
        for (k in seq_along(true_indirect)) {
          lvl <- true_indirect[[k]]
          graphics::segments(xrng[1], lvl, xrng[2], lvl,
                             col = "gray30", lty = 2, lwd = LW_RULE)
        }
      }
      
      for (i in seq_len(n_models)) {
        d <- ind_list[[i]]
        if (is.null(d)) next
        ltype <- if (length(d$time) == 1L) "p" else "l"
        faded <- grDevices::adjustcolor(cols[i], alpha.f = 1)
        for (j in seq_along(d$donors)) {
          graphics::lines(d$time, d$mean[, j], col = faded,
                          lwd = LW_DENS, type = ltype, pch = 16)
        }
      }
      
      graphics::box()
      
      if (isTRUE(show_legend)) {
        leg_keep <- if (identical(indirect_legend, "all")) {
          rep(TRUE, n_models)
        } else {
          has_ind
        }
        graphics::par(mar = c(0, 0, 0, 0), family = font_family)
        graphics::plot.new()
        graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
        graphics::legend("center",
                         legend = names(models)[leg_keep],
                         col    = cols[leg_keep],
                         lty    = 1, lwd = LW_PATH,
                         horiz  = TRUE, bty = "n", cex = LEGEND_CEX)
      }
    }
  }
  
  invisible(NULL)
}

# Single MC Iterations

W_ws <- generate_watts_strogatz_matrix(N = N, k = 2, p = ws_p, seed = 13)

set.seed(1111)
sim <- generate_data_ws_planted(
  W              = W_ws,
  type           = dgp_type,
  N              = N,
  T              = T,
  T_0            = T_0,
  treated_idx    = treated_idx,
  beta           = beta,
  theta          = theta,
  delta_mean     = delta_mean,
  delta_sd       = delta_sd,
  rho            = rho_true,
  sigma_u        = sigma_u,
  X_mean         = X_mean,
  x_sd           = x_sd,
  alpha_sd       = alpha_sd,
  twin_target    = twin_target,
  weight_profile = weight_profile
)

df <- sim$df
W  <- sim$W
covariates <- df %>% dplyr::select(time, id, X1, X2)

fit_plain <- nascSynth$new(
  data            = df,
  time            = time,
  id              = id,
  treated         = D,
  outcome         = Y,
  covariates      = covariates,
  W               = W,
  spatial_model   = "exogenous",
  rho             = rho_true,
  bias_correction = FALSE,
  nasc_penalty    = FALSE,
  ci_width        = 0.95
)
invisible(utils::capture.output(suppressMessages(
  fit_plain$fit(n_samples = 50, cores = N_WORKERS,
                iter = stan_iter, warmup = stan_warmup, control = stan_control))))

fit_pen <- nascSynth$new(
  data              = df,
  time              = time,
  id                = id,
  treated           = D,
  outcome           = Y,
  covariates        = covariates,
  W                 = W,
  spatial_model     = "exogenous",
  rho               = rho_true,
  bias_correction   = FALSE,
  nasc_penalty      = TRUE,
  ci_width          = 0.95,
  lambda_cv_grid    = lambda_cv_grid,
  lambda_train_frac = lambda_train_frac
)
invisible(utils::capture.output(suppressMessages(
  fit_pen$fit(n_samples = 50, cores = N_WORKERS,
              iter = stan_iter, warmup = stan_warmup, control = stan_control))))

fit_bc <- nascSynth$new(
  data            = df,
  time            = time,
  id              = id,
  treated         = D,
  outcome         = Y,
  covariates      = covariates,
  W               = W,
  spatial_model   = "exogenous",
  rho             = rho_true,
  bias_correction = TRUE,
  nasc_penalty    = FALSE,
  ci_width        = 0.95
)
invisible(utils::capture.output(suppressMessages(
  fit_bc$fit(n_samples = 25, cores = N_WORKERS,
             iter = stan_iter, warmup = stan_warmup, control = stan_control))))

fit_nasc <- nascSynth$new(
  data              = df,
  time              = time,
  id                = id,
  treated           = D,
  outcome           = Y,
  covariates        = covariates,
  W                 = W,
  spatial_model     = "exogenous",
  rho               = rho_true,
  bias_correction   = TRUE,
  nasc_penalty      = TRUE,
  ci_width          = 0.95,
  lambda_cv_grid    = lambda_cv_grid,
  lambda_train_frac = lambda_train_frac
)
invisible(utils::capture.output(suppressMessages(
  fit_nasc$fit(n_samples = 25, cores = N_WORKERS,
               iter = stan_iter, warmup = stan_warmup, control = stan_control))))

models_list <- list(
  "Conventional SC" = fit_plain,
  "NASC Penalty"    = fit_pen,
  "Bias Correction" = fit_bc,
  "Full NASC"       = fit_nasc
)

scm_synth <- synthSCM(df, treated_id = treated_idx,
                      intervention_time = fit_plain$interventionTime,
                      ci_width = 0.95)

att_true      <- delta_mean
ind_donor_ids <- as.character(setdiff(seq_len(N), treated_idx))
M_mult        <- solve(diag(N) - rho_true * W)
s_true        <- M_mult[, treated_idx] / M_mult[treated_idx, treated_idx]
names(s_true) <- as.character(seq_len(N))
s_true        <- s_true[ind_donor_ids]
ind_true      <- att_true * s_true

nascPlot(list("SCM" = scm_synth, "BSCM" = fit_plain,
              "CR" = fit_pen, "BC" = fit_bc, "NASC" = fit_nasc),
         show_ci = TRUE, ci_models = c("SCM", "NASC"),
         true_att = att_true)

nascPlot(list("SCM" = scm_synth, "BSCM" = fit_plain,
              "CR" = fit_pen, "BC" = fit_bc, "NASC" = fit_nasc),
         panels = "indirect", indirect_models = c("BC", "NASC"),
         true_indirect = ind_true)

all_ids       <- as.character(unique(df$id))
num_try       <- suppressWarnings(as.numeric(all_ids))
all_ids       <- if (!any(is.na(num_try))) all_ids[order(num_try)] else sort(all_ids)
donor_ids_all <- setdiff(all_ids, as.character(treated_idx))
true_w        <- setNames(rep(0, length(donor_ids_all)), donor_ids_all)
true_w[as.character(sim$planted_donors)] <- weight_profile

w_table <- nascWeightSCM(list("BSCM" = fit_plain, "NASC" = fit_nasc),
                         synth_w = scm_synth$weights,
                         true_w  = true_w)

save_fig <- function(fname, expr, plot_h = 2, width = 6,
                     show_legend = TRUE) {
  dev_h <- if (show_legend) plot_h * (6 + LEGEND_H) / 6 else plot_h
  pdf(file.path(dir_fig, fname), width = width, height = dev_h,
      pointsize = 9, family = PLOT_FONT)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  invisible(NULL)
}



mods <- list("SCM" = scm_synth, "BSCM" = fit_plain,
             "CR" = fit_pen, "BC" = fit_bc, "NASC" = fit_nasc)

save_fig("example_synth_paths.pdf",
         nascPlot(mods, show_ci = TRUE, ci_models = c("SCM", "NASC"),
                  panels = "path", show_legend = FALSE), show_legend = FALSE)
save_fig("example_effects.pdf",
         nascPlot(mods, show_ci = TRUE, ci_models = c("SCM", "NASC"),
                  panels = "effect", show_legend = FALSE,
                  true_att = att_true),
         show_legend = FALSE)
save_fig("example_indirect.pdf",
         nascPlot(mods, panels = "indirect",
                  indirect_models = c("BC", "NASC", "CR"),
                  indirect_legend = "all",
                  true_indirect = ind_true))
save_fig("example_weights.pdf",
         nascWeightSCM(list("BSCM" = fit_plain, "NASC" = fit_nasc),
                       synth_w = scm_synth$weights, true_w = true_w),
         plot_h = 5)
