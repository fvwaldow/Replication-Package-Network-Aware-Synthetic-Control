# Monte Carlo Simulation - SCM

library(Synth)
library(dplyr)
library(tidyr)
library(scales)
library(tibble)
library(future.apply)
library(igraph)

B <- 1000
T_seq <- c(30,60,90)
N_seq <- c(15)
rho_seq <- c(0.8, 0.6, 0.4, 0.2, 0, -0.2, -0.4, -0.6, -0.8)

# Watts-Strogatz network
k_seq                    <- c(2)
p_seq                    <- c(0, 0.1, 0.4)
seed_w                   <- 13

dgp_type                 <- "SAR"

treated_idx              <- 1
beta                     <- c(1.0, 0.5)
theta                    <- c(0.3, 0.2)
delta_mean               <- 5
delta_sd                 <- 1
sigma_u                  <- 0.5
alpha_sd                 <- 1
x_sd                     <- 1
X_mean                   <- c(0.0, 0.0)
twin_target              <- "cleanest"

# weights, low |s| -> high |s|
weight_profile           <- c(0.075, 0.075, 0.150, 0.150, 0.250, 0.300)

# SCM benchmark settings
scm_v_mode               <- "equal"
scm_placebo              <- TRUE
scm_ci                   <- 0.95
scm_placebo_prefit_mult  <- Inf

SCM_INFERENCE_VERSION    <- "adh_rmspe_standardized_v2"

# dataset export
DATA_DIR       <- file.path(directory, "output", "data_ex_rho_scm")
RESULTS_RDS    <- file.path(directory, "output", "data_ex_rho_scm.rds")
save_datasets  <- TRUE
save_csv       <- TRUE
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# estimators
estimators <- list(
  sc_classic = list(engine = "scm")
)


# MC simulation

# SCM helpers

.synth_one <- function(foo, treated_id, control_ids, pre_times, all_times,
                       first_post, v_mode = "equal") {
  sp <- lapply(pre_times, function(tt) list("Y", tt, "mean"))
  dprep <- NULL
  invisible(utils::capture.output(
    dprep <- Synth::dataprep(
      foo                   = foo,
      predictors            = c("X1", "X2"),
      predictors.op         = "mean",
      special.predictors    = sp,
      dependent             = "Y",
      unit.variable         = "id",
      time.variable         = "time",
      treatment.identifier  = treated_id,
      controls.identifier   = as.numeric(control_ids),
      time.predictors.prior = pre_times,
      time.optimize.ssr     = pre_times,
      time.plot             = all_times
    )
  ))
  
  sout <- NULL
  if (identical(v_mode, "equal")) {
    np <- nrow(dprep$X1)
    invisible(utils::capture.output(
      sout <- Synth::synth(dprep, custom.v = rep(1 / np, np))
    ))
  } else {
    invisible(utils::capture.output(sout <- Synth::synth(dprep)))
  }
  
  w   <- sout$solution.w
  eff <- as.numeric(dprep$Y1plot[, 1]) - as.numeric(dprep$Y0plot %*% w)
  tp  <- as.numeric(rownames(dprep$Y1plot))
  pre_m  <- tp <  first_post
  post_m <- tp >= first_post
  list(
    att       = mean(eff[post_m]),
    pre_rmse  = sqrt(mean(eff[pre_m]^2)),
    post_rmse = sqrt(mean(eff[post_m]^2)),
    w         = w
  )
}

# SCM with in-space placebo inference (ADH 2015; Firpo-Possebom CI)
fit_conventional_sc <- function(df, donor_idx, treated_id, ci_width = 0.95,
                                placebo = TRUE, placebo_prefit_mult = Inf,
                                v_mode = "equal") {
  foo <- as.data.frame(df)
  foo$id   <- as.numeric(as.character(foo$id))
  foo$time <- as.numeric(as.character(foo$time))
  
  td <- foo[foo$id == treated_id, c("time", "D")]
  td <- td[order(td$time), , drop = FALSE]
  post_times <- td$time[td$D == 1]
  if (length(post_times) == 0L) stop("Synth SCM: no post periods (D==1).")
  first_post <- min(post_times)
  
  all_times <- sort(unique(foo$time))
  pre_times <- all_times[all_times < first_post]
  if (length(pre_times) < 2L) stop("Synth SCM: too few pre periods.")
  
  donor_ids <- as.numeric(donor_idx)
  
  # treated-unit fit
  main      <- .synth_one(foo, treated_id, donor_ids, pre_times, all_times,
                          first_post, v_mode = v_mode)
  att_hat   <- main$att
  pre_rmse  <- main$pre_rmse
  post_rmse <- main$post_rmse
  
  w_named   <- as.numeric(main$w); names(w_named) <- rownames(main$w)
  w_aligned <- as.numeric(w_named[as.character(donor_ids)])
  
  ratio_hat <- if (is.finite(pre_rmse) && pre_rmse > 0) post_rmse / pre_rmse else NA_real_
  
  # placebos
  att_lower     <- NA_real_; att_upper     <- NA_real_; att_sd <- NA_real_
  att_lower_raw <- NA_real_; att_upper_raw <- NA_real_
  p_value       <- NA_real_
  n_placebo     <- 0L
  if (isTRUE(placebo)) {
    pl_att  <- numeric(0)
    pl_pre  <- numeric(0)
    pl_post <- numeric(0)
    for (d in donor_ids) {
      ctrl <- setdiff(donor_ids, d)
      if (length(ctrl) < 2L) next
      pl <- tryCatch(
        .synth_one(foo, d, ctrl, pre_times, all_times, first_post,
                   v_mode = v_mode),
        error = function(e) NULL
      )
      if (is.null(pl)) next
      if (!is.finite(pl$att) || !is.finite(pl$pre_rmse) ||
          !is.finite(pl$post_rmse)) next
      if (pl$pre_rmse <= sqrt(.Machine$double.eps)) next
      if (is.finite(placebo_prefit_mult) &&
          pl$pre_rmse > placebo_prefit_mult * pre_rmse) next
      pl_att  <- c(pl_att,  pl$att)
      pl_pre  <- c(pl_pre,  pl$pre_rmse)
      pl_post <- c(pl_post, pl$post_rmse)
    }
    n_placebo <- length(pl_att)
    
    if (n_placebo >= 2L) {
      a <- (1 - ci_width) / 2
      
      if (is.finite(pre_rmse) && pre_rmse > 0) {
        gs <- pl_att / pl_pre
        qs <- as.numeric(stats::quantile(gs, probs = c(a, 1 - a),
                                         names = FALSE, type = 7))
        att_lower <- att_hat - qs[2] * pre_rmse
        att_upper <- att_hat - qs[1] * pre_rmse
        att_sd    <- stats::sd(gs) * pre_rmse
      }
      
      qr <- as.numeric(stats::quantile(pl_att, probs = c(a, 1 - a),
                                       names = FALSE, type = 7))
      att_lower_raw <- att_hat - qr[2]
      att_upper_raw <- att_hat - qr[1]
    }
    
    if (n_placebo >= 1L && is.finite(ratio_hat)) {
      ratios <- pl_post / pl_pre
      ratios <- ratios[is.finite(ratios)]
      if (length(ratios) >= 1L) {
        p_value <- mean(c(ratios, ratio_hat) >= ratio_hat)
      }
    }
  }
  
  list(
    att = c(mean = att_hat, sd = att_sd, lower = att_lower, upper = att_upper,
            lower_raw = att_lower_raw, upper_raw = att_upper_raw),
    weights = data.frame(donor = as.character(donor_ids),
                         mean  = w_aligned,
                         stringsAsFactors = FALSE),
    pre_rmse    = pre_rmse,
    post_rmse   = post_rmse,
    rmspe_ratio = ratio_hat,
    p_value     = p_value,
    n_placebo   = n_placebo
  )
}

# functions
extract_weights_aligned <- function(fit_summary, donor_idx) {
  w_tbl <- fit_summary$weights
  if (is.null(w_tbl) || nrow(w_tbl) == 0L) return(NULL)
  donor_ids_chr <- as.character(donor_idx)
  ord <- match(donor_ids_chr, w_tbl$donor)
  if (anyNA(ord)) {
    w_hat <- rep(NA_real_, length(donor_idx))
    keep  <- !is.na(ord)
    w_hat[keep] <- w_tbl$mean[ord[keep]]
  } else {
    w_hat <- w_tbl$mean[ord]
  }
  w_hat
}

weight_recovery_metrics <- function(w_hat, w_star_donor, s_abs) {
  if (is.null(w_hat) || any(is.na(w_hat))) {
    return(data.frame(
      w_bias_l1        = NA_real_,
      w_bias_l2        = NA_real_,
      w_bias_max       = NA_real_,
      w_sum            = NA_real_,
      contam_mass_hat  = NA_real_,
      contam_mass_star = NA_real_
    ))
  }
  diff_v <- w_hat - w_star_donor
  data.frame(
    w_bias_l1        = sum(abs(diff_v)),
    w_bias_l2        = sqrt(sum(diff_v^2)),
    w_bias_max       = max(abs(diff_v)),
    w_sum            = sum(w_hat),
    contam_mass_hat  = if (length(s_abs) == length(w_hat)) sum(abs(w_hat) * s_abs) else NA_real_,
    contam_mass_star = if (length(s_abs) == length(w_star_donor)) sum(abs(w_star_donor) * s_abs) else NA_real_
  )
}

# dataset path helpers
.net_tag <- function(dgp_type, ws_k, ws_p, N) {
  sprintf("%s_k%d_p%.2f_N%d", dgp_type, ws_k, ws_p, N)
}
cell_data_dir <- function(data_dir, dgp_type, ws_k, ws_p, N, rho) {
  d <- file.path(data_dir, .net_tag(dgp_type, ws_k, ws_p, N),
                 sprintf("rho_%.2f", rho))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
dataset_basename <- function(T, rep) {
  sprintf("dataset_T%d_rep%03d", T, rep)
}

# worker: one replicate
run_one_rep <- function(b, cell, W, estimators, dgp_type = "SAR",
                        data_dir = NULL, save_datasets = FALSE, save_csv = TRUE) {
  
  sim <- generate_data_ws_planted(
    W              = W,
    type           = dgp_type,
    N              = cell$N,
    T              = cell$T + 10,
    T_0            = cell$T,
    treated_idx    = cell$treated_idx,
    beta           = cell$beta,
    theta          = cell$theta,
    delta_mean     = cell$delta_mean,
    delta_sd       = cell$delta_sd,
    rho            = cell$rho,
    sigma_u        = cell$sigma_u,
    X_mean         = cell$X_mean,
    x_sd           = cell$x_sd,
    alpha_sd       = cell$alpha_sd,
    twin_target    = cell$twin_target,
    weight_profile = cell$weight_profile
  )
  df <- sim$df
  
  att_true_realized <- mean(sim$true_att)
  
  w_star_donor <- sim$w_star_donor
  donor_idx    <- sim$donor_idx
  s_abs        <- if (!is.null(sim$contam)) sim$contam$s_abs else rep(0, length(donor_idx))
  
  rows <- vector("list", length(estimators))
  for (i in seq_along(estimators)) {
    est_name <- names(estimators)[i]
    cfg      <- estimators[[i]]
    
    t0 <- proc.time()[["elapsed"]]
    fit_ok <- tryCatch({
      s_obj <- switch(cfg$engine,
                      "scm" = {
                        fit_conventional_sc(
                          df                  = df,
                          donor_idx           = donor_idx,
                          treated_id          = cell$treated_idx,
                          ci_width            = cell$scm_ci,
                          placebo             = cell$scm_placebo,
                          placebo_prefit_mult = cell$scm_placebo_prefit_mult,
                          v_mode              = cell$scm_v_mode
                        )
                      },
                      stop("unknown estimator: ", cfg$engine)
      )
      list(s = s_obj, ok = TRUE, msg = NA_character_)
    },
    error = function(e) list(s = NULL, ok = FALSE, msg = conditionMessage(e)))
    rt <- proc.time()[["elapsed"]] - t0
    
    base <- data.frame(
      ws_k         = cell$ws_k,
      ws_p         = cell$ws_p,
      N            = cell$N,
      T            = cell$T,
      rho          = cell$rho,
      rep          = b,
      estimator    = est_name,
      att_true     = att_true_realized,
      att_struct   = cell$delta_mean,
      runtime_s    = rt,
      stringsAsFactors = FALSE
    )
    
    if (fit_ok$ok) {
      att    <- fit_ok$s$att
      lo     <- as.numeric(att["lower"])
      hi     <- as.numeric(att["upper"])
      lo_raw <- as.numeric(att["lower_raw"])
      hi_raw <- as.numeric(att["upper_raw"])
      
      w_hat <- extract_weights_aligned(fit_ok$s, donor_idx)
      wmet  <- weight_recovery_metrics(w_hat, w_star_donor, s_abs)
      
      rows[[i]] <- cbind(base, data.frame(
        att_hat           = as.numeric(att["mean"]),
        att_sd            = as.numeric(att["sd"]),
        att_lower         = lo,
        att_upper         = hi,
        covers            = (att_true_realized >= lo) & (att_true_realized <= hi),
        covers_struct     = (cell$delta_mean    >= lo) & (cell$delta_mean    <= hi),
        ci_width          = hi - lo,
        att_lower_raw     = lo_raw,
        att_upper_raw     = hi_raw,
        covers_raw        = (att_true_realized >= lo_raw) & (att_true_realized <= hi_raw),
        covers_struct_raw = (cell$delta_mean    >= lo_raw) & (cell$delta_mean    <= hi_raw),
        ci_width_raw      = hi_raw - lo_raw,
        p_value           = fit_ok$s$p_value,
        pre_rmse      = fit_ok$s$pre_rmse,
        post_rmse     = fit_ok$s$post_rmse,
        rmspe_ratio   = fit_ok$s$rmspe_ratio,
        n_placebo     = fit_ok$s$n_placebo,
        status        = "ok",
        err_msg       = NA_character_,
        stringsAsFactors = FALSE
      ), wmet)
      rows[[i]]$w_hat        <- I(list(w_hat))
      rows[[i]]$w_star_donor <- I(list(w_star_donor))
      rows[[i]]$donor_idx    <- I(list(donor_idx))
      rows[[i]]$s_abs        <- I(list(s_abs))
    } else {
      empty_wmet <- weight_recovery_metrics(NULL, NULL, NULL)
      rows[[i]] <- cbind(base, data.frame(
        att_hat           = NA_real_,
        att_sd            = NA_real_,
        att_lower         = NA_real_,
        att_upper         = NA_real_,
        covers            = NA,
        covers_struct     = NA,
        ci_width          = NA_real_,
        att_lower_raw     = NA_real_,
        att_upper_raw     = NA_real_,
        covers_raw        = NA,
        covers_struct_raw = NA,
        ci_width_raw      = NA_real_,
        p_value           = NA_real_,
        pre_rmse      = NA_real_,
        post_rmse     = NA_real_,
        rmspe_ratio   = NA_real_,
        n_placebo     = NA_integer_,
        status        = "error",
        err_msg       = fit_ok$msg,
        stringsAsFactors = FALSE
      ), empty_wmet)
      rows[[i]]$w_hat        <- I(list(NA))
      rows[[i]]$w_star_donor <- I(list(w_star_donor))
      rows[[i]]$donor_idx    <- I(list(donor_idx))
      rows[[i]]$s_abs        <- I(list(s_abs))
    }
  }
  out <- dplyr::bind_rows(rows)
  
  if (isTRUE(save_datasets) && !is.null(data_dir)) {
    cdir  <- cell_data_dir(data_dir, dgp_type, cell$ws_k, cell$ws_p,
                           cell$N, cell$rho)
    bname <- dataset_basename(cell$T, b)
    saveRDS(list(sim = sim, results = out,
                 scm_inference = SCM_INFERENCE_VERSION),
            file.path(cdir, paste0(bname, ".rds")))
    if (isTRUE(save_csv)) {
      write.csv(df, file.path(cdir, paste0(bname, ".csv")), row.names = FALSE)
      res_scalar <- out[, !vapply(out, is.list, logical(1)), drop = FALSE]
      write.csv(res_scalar, file.path(cdir, paste0(bname, "_results.csv")),
                row.names = FALSE)
    }
  }
  
  out
}

# run
try(future:::ClusterRegistry("stop"), silent = TRUE)
gc()

max_workers <- 6L
n_workers <- min(max_workers, max(1L, parallel::detectCores() - 1L))
options(future.globals.maxSize = 2 * 1024^3)

plan(multisession, workers = n_workers)
cat(sprintf("%d cores\n", n_workers))

# build weight matrices
weights_list <- array(list(),
                      dim = c(length(k_seq), length(p_seq), length(N_seq)))

for (ki in seq_along(k_seq)) {
  for (pi in seq_along(p_seq)) {
    for (j in seq_along(N_seq)) {
      W_ws <- generate_watts_strogatz_matrix(
        N    = N_seq[j],
        k    = k_seq[ki],
        p    = p_seq[pi],
        seed = seed_w
      )
      weights_list[[ki, pi, j]] <- W_ws
      
      if (save_datasets) {
        net_dir <- file.path(DATA_DIR, .net_tag(dgp_type, k_seq[ki],
                                                p_seq[pi], N_seq[j]))
        dir.create(net_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(W_ws, file.path(net_dir, "network_W.rds"))
        if (save_csv) {
          write.csv(as.data.frame(as.matrix(W_ws)),
                    file.path(net_dir, "network_W.csv"), row.names = FALSE)
        }
        for (rr in rho_seq) {
          dir.create(file.path(net_dir, sprintf("rho_%.2f", rr)),
                     recursive = TRUE, showWarnings = FALSE)
        }
      }
    }
  }
}

# main loop
results <- list()

mc_t_start <- proc.time()[["elapsed"]]

for (ki in seq_along(k_seq)) {
  for (pi in seq_along(p_seq)) {
    for (j in seq_along(N_seq)) {
      W_cell <- weights_list[[ki, pi, j]]
      
      for (t in seq_along(T_seq)) {
        for (r in seq_along(rho_seq)) {
          cell <- list(
            ws_k                    = k_seq[ki],
            ws_p                    = p_seq[pi],
            N                       = N_seq[j],
            T                       = T_seq[t],
            rho                     = rho_seq[r],
            treated_idx             = treated_idx,
            beta                    = beta,
            theta                   = theta,
            delta_mean              = delta_mean,
            delta_sd                = delta_sd,
            sigma_u                 = sigma_u,
            alpha_sd                = alpha_sd,
            X_mean                  = X_mean,
            x_sd                    = x_sd,
            twin_target             = twin_target,
            weight_profile          = weight_profile,
            scm_v_mode              = scm_v_mode,
            scm_placebo             = scm_placebo,
            scm_ci                  = scm_ci,
            scm_placebo_prefit_mult = scm_placebo_prefit_mult
          )
          
          cat(sprintf("Cell k=%d p=%.2f j=%d t=%d r=%d : %d reps on %d cores\n",
                      k_seq[ki], p_seq[pi], j, t, r, B, n_workers))
          t_cell <- proc.time()[["elapsed"]]
          
          cell_rows <- tryCatch(
            future_lapply(
              seq_len(B),
              function(b) {
                if (!exists("generate_data_ws_planted", envir = globalenv(), inherits = FALSE)) {
                  source(dgp_path)
                }
                if (isTRUE(save_datasets)) {
                  rp <- file.path(
                    cell_data_dir(DATA_DIR, dgp_type, cell$ws_k, cell$ws_p, cell$N, cell$rho),
                    paste0(dataset_basename(cell$T, b), ".rds")
                  )
                  if (file.exists(rp)) {
                    prev <- tryCatch(readRDS(rp), error = function(e) NULL)
                    if (is.list(prev) && !is.null(prev$results) &&
                        identical(prev$scm_inference, SCM_INFERENCE_VERSION)) {
                      return(prev$results)
                    }
                  }
                }
                run_one_rep(b, cell, W_cell, estimators, dgp_type = dgp_type,
                            data_dir = DATA_DIR, save_datasets = save_datasets,
                            save_csv = save_csv)
              },
              future.seed     = 1234L,
              future.packages = c("Synth", "dplyr", "tidyr"),
              future.globals  = c("run_one_rep", "cell", "W_cell",
                                  "estimators", "dgp_type", "dgp_path",
                                  "extract_weights_aligned",
                                  "weight_recovery_metrics",
                                  ".synth_one", "fit_conventional_sc",
                                  "DATA_DIR", "save_datasets", "save_csv",
                                  "SCM_INFERENCE_VERSION",
                                  ".net_tag", "cell_data_dir", "dataset_basename")
            ),
            error = function(e) {
              plan(sequential)
              try(future:::ClusterRegistry("stop"), silent = TRUE)
              stop("MC error", call. = FALSE)
            }
          )
          
          cat(sprintf("%.1fs\n", proc.time()[["elapsed"]] - t_cell))
          
          results[[length(results) + 1L]] <- dplyr::bind_rows(cell_rows)
          saveRDS(dplyr::bind_rows(results), RESULTS_RDS)
          gc()
        }
      }
    }
  }
}

plan(sequential)
try(future:::ClusterRegistry("stop"), silent = TRUE)

mc_elapsed <- proc.time()[["elapsed"]] - mc_t_start
cat(sprintf("\n%.1fs (%.2f min).\n",
            mc_elapsed, mc_elapsed / 60))
