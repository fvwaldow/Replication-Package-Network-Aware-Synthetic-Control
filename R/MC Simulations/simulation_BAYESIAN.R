# Monte Carlo Simulation - Bayesian (NASC, BC. CR, BSCM)

library(nasc)
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

lambda_cv_grid           <- c(0, exp(seq(log(0.05), log(50), length.out = 20)))
lambda_train_frac        <- 0.8

# sampler settings
stan_iter                <- 1000L
stan_warmup              <- 500L
stan_control             <- list(adapt_delta = 0.95, max_treedepth = 10)

# dataset export
DATA_DIR       <- file.path(directory, "output", "data_ex_rho_2")
RESULTS_RDS    <- file.path(directory, "output", "data_ex_rho_2.rds")
save_datasets  <- TRUE
save_csv       <- TRUE
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# estimators
estimators <- list(
  plain_est  = list(engine = "nasc", bias_correction = FALSE, nasc_penalty = FALSE),
  bc_est     = list(engine = "nasc", bias_correction = TRUE,  nasc_penalty = FALSE),
  reg_est    = list(engine = "nasc", bias_correction = FALSE, nasc_penalty = TRUE),
  bc_reg_est = list(engine = "nasc", bias_correction = TRUE,  nasc_penalty = TRUE)
)
# MC simulation

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
  
  covariates <- df %>% dplyr::select(time, id, X1, X2)
  
  rows <- vector("list", length(estimators))
  for (i in seq_along(estimators)) {
    est_name <- names(estimators)[i]
    cfg      <- estimators[[i]]
    
    t0 <- proc.time()[["elapsed"]]
    fit_ok <- tryCatch({
      s_obj <- switch(cfg$engine,
                      "nasc" = {
                        ctor <- list(
                          data            = df,
                          time            = quote(time),
                          id              = quote(id),
                          treated         = quote(D),
                          outcome         = quote(Y),
                          covariates      = covariates,
                          W               = W,
                          spatial_model   = "exogenous",
                          rho             = cell$rho,
                          bias_correction = cfg$bias_correction,
                          nasc_penalty    = cfg$nasc_penalty,
                          ci_width        = 0.95
                        )
                        if (isTRUE(cfg$nasc_penalty)) {
                          ctor$lambda_cv_grid    <- cell$lambda_cv_grid
                          ctor$lambda_train_frac <- cell$lambda_train_frac
                        }
                        fit <- do.call(nascSynth$new, ctor)
                        fit$fit(cores   = 1,
                                iter    = cell$stan_iter,
                                warmup  = cell$stan_warmup,
                                control = cell$stan_control)
                        fit$summary(print = FALSE)
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
      att <- fit_ok$s$att
      lo  <- as.numeric(att["lower"])
      hi  <- as.numeric(att["upper"])
      
      w_hat <- extract_weights_aligned(fit_ok$s, donor_idx)
      wmet  <- weight_recovery_metrics(w_hat, w_star_donor, s_abs)
      
      lam_used <- NA_real_
      if (isTRUE(cfg$nasc_penalty)) {
        lam_used <- tryCatch(as.numeric(fit$lambdaCV()$lambda),
                             error = function(e) NA_real_)
      }
      
      rows[[i]] <- cbind(base, data.frame(
        att_hat       = as.numeric(att["mean"]),
        att_sd        = as.numeric(att["sd"]),
        att_lower     = lo,
        att_upper     = hi,
        covers        = (att_true_realized >= lo) & (att_true_realized <= hi),
        covers_struct = (cell$delta_mean    >= lo) & (cell$delta_mean    <= hi),
        ci_width      = hi - lo,
        pre_rmse      = fit_ok$s$pre_rmse,
        post_rmse     = fit_ok$s$post_rmse,
        rmspe_ratio   = fit_ok$s$rmspe_ratio,
        lambda_cv     = lam_used,
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
        att_hat       = NA_real_,
        att_sd        = NA_real_,
        att_lower     = NA_real_,
        att_upper     = NA_real_,
        covers        = NA,
        covers_struct = NA,
        ci_width      = NA_real_,
        pre_rmse      = NA_real_,
        post_rmse     = NA_real_,
        rmspe_ratio   = NA_real_,
        lambda_cv     = NA_real_,
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
    saveRDS(list(sim = sim, results = out),
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
            ws_k              = k_seq[ki],
            ws_p              = p_seq[pi],
            N                 = N_seq[j],
            T                 = T_seq[t],
            rho               = rho_seq[r],
            treated_idx       = treated_idx,
            beta              = beta,
            theta             = theta,
            delta_mean        = delta_mean,
            delta_sd          = delta_sd,
            sigma_u           = sigma_u,
            alpha_sd          = alpha_sd,
            X_mean            = X_mean,
            x_sd              = x_sd,
            twin_target       = twin_target,
            weight_profile    = weight_profile,
            lambda_cv_grid    = lambda_cv_grid,
            lambda_train_frac = lambda_train_frac,
            stan_iter         = stan_iter,
            stan_warmup       = stan_warmup,
            stan_control      = stan_control
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
                    if (is.list(prev) && !is.null(prev$results)) {
                      return(prev$results)
                    }
                  }
                }
                run_one_rep(b, cell, W_cell, estimators, dgp_type = dgp_type,
                            data_dir = DATA_DIR, save_datasets = save_datasets,
                            save_csv = save_csv)
              },
              future.seed     = 1234L,
              future.packages = c("nasc", "dplyr", "tidyr"),
              future.globals  = c("run_one_rep", "cell", "W_cell",
                                  "estimators", "dgp_type", "dgp_path",
                                  "extract_weights_aligned",
                                  "weight_recovery_metrics",
                                  "DATA_DIR", "save_datasets", "save_csv",
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
