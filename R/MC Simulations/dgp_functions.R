library(igraph)

# Watts-Strogatz network matrix
generate_watts_strogatz_matrix <- function(N = 15, k = 2, p = 0.1, seed = 13) {
  build <- function() {
    ws_graph <- sample_smallworld(dim = 1, size = N, nei = k, p = p)
    W <- as.matrix(as_adjacency_matrix(ws_graph))
    diag(W) <- 0
    rs <- rowSums(W); rs[rs == 0] <- 1
    W / rs
  }
  if (is.null(seed)) build() else withr::with_seed(seed, build())
}

# contamination vector
contamination_vector <- function(W, treated_idx, rho) {
  donor_idx <- setdiff(seq_len(nrow(W)), treated_idx)
  W_J  <- W[donor_idx, donor_idx]
  w_J1 <- W[donor_idx, treated_idx]
  J    <- length(donor_idx)
  s    <- rho * solve(diag(J) - rho * W_J, w_J1)
  list(
    s         = as.numeric(s),
    s_abs     = abs(as.numeric(s)),
    summary   = summary(abs(as.numeric(s))),
    cv        = if (mean(abs(s)) > 0) sd(abs(s)) / mean(abs(s)) else NA_real_,
    donor_idx = donor_idx
  )
}

# pre-specified synthetic weights
pick_planted_weights <- function(W, treated_idx, rho,
                                 weight_profile = c(0.075, 0.075,
                                                    0.150, 0.150,
                                                    0.250, 0.300)
) {
  n_planted <- length(weight_profile)
  cv        <- contamination_vector(W, treated_idx, rho)
  donor_idx <- cv$donor_idx
  s_abs     <- cv$s_abs
  J         <- length(donor_idx)

  if (n_planted > J)
    stop(sprintf("Need %d planted donors but only %d donors available.", n_planted, J))

  ord_asc    <- order(s_abs, donor_idx)
  picked_pos <- ord_asc[(J - n_planted + 1):J]
  picked_idx <- donor_idx[picked_pos]

  w_star <- numeric(nrow(W)); w_star[picked_idx] <- weight_profile
  w_star_donor <- numeric(J); w_star_donor[picked_pos] <- weight_profile

  list(w_star = w_star, w_star_donor = w_star_donor,
       planted_donors = picked_idx, weight_profile = weight_profile,
       donor_idx = donor_idx, s_abs_picked = s_abs[picked_pos],
       s_abs_all = s_abs, picked_pos = picked_pos,
       ord_asc = ord_asc, rank_picked = (J - n_planted + 1):J)
}

# Main DGP: TWIN donors across the contamination extremes
generate_data_ws_planted <- function(W,
                                     type           = "SAR",
                                     N              = 15,
                                     T              = 30,
                                     T_0            = 20,
                                     treated_idx    = 1,
                                     beta           = c(1.0, 0.5),
                                     theta          = c(0.3, 0.2),
                                     delta_mean     = 5.0,
                                     delta_sd       = 1.0,
                                     rho            = 0.6,
                                     sigma_u        = 0.30,
                                     X_mean         = c(0.0, 0.0),
                                     x_sd           = 0.50,
                                     alpha_sd       = 1.0,
                                     twin_target    = "modest",
                                     weight_profile = c(0.075, 0.075,
                                                        0.150, 0.150,
                                                        0.250, 0.300)
) {

  J0 <- nrow(W) - 1L
  n_twins <- 6L
  twin_gap <- 0.0

  if (length(X_mean) == 1) X_mean <- c(X_mean, X_mean)

  rho_for_pick <- if (type == "SLX") 0 else rho
  planted   <- pick_planted_weights(W, treated_idx, rho_for_pick,
                                    weight_profile = weight_profile)
  donor_idx <- planted$donor_idx
  w_star_J  <- planted$w_star_donor
  J         <- length(donor_idx)
  s_abs_J   <- planted$s_abs_all
  ord_asc   <- planted$ord_asc
  picked_pos <- planted$picked_pos

  n_planted <- length(weight_profile)
  if ((n_twins + n_planted) > J)
    stop(sprintf("Hardcoded n_twins (6) + planted (%d) exceeds J (%d).", n_planted, J))

  if (!twin_target %in% c("modest", "cleanest"))
    stop("'twin_target' must be \"modest\" or \"cleanest\".")

  avail_pos <- setdiff(ord_asc, picked_pos)
  if (n_twins > length(avail_pos))
    stop(sprintf("Hardcoded n_twins (6) exceeds the %d non-planted donors available.", length(avail_pos)))

  high_pos <- rev(ord_asc)[seq_len(n_twins)]
  low_pos  <- if (twin_target == "modest") {
    rev(avail_pos[(length(avail_pos) - n_twins + 1L):length(avail_pos)])
  } else {
    avail_pos[seq_len(n_twins)]
  }
  high_idx <- donor_idx[high_pos]
  low_idx  <- donor_idx[low_pos]

  W_J <- W[donor_idx, donor_idx]
  w_1 <- W[donor_idx, treated_idx]
  M   <- if (type %in% c("SAR", "SDM")) solve(diag(J) - rho * W_J) else diag(J)

  X1_mat <- matrix(rnorm(N * T, X_mean[1], x_sd), nrow = N, ncol = T)
  X2_mat <- matrix(rnorm(N * T, X_mean[2], x_sd), nrow = N, ncol = T)
  alpha  <- rnorm(N, 0, alpha_sd)

  for (i in seq_len(n_twins)) {
    h <- high_idx[i]; l <- low_idx[i]
    alpha[l]     <- alpha[h]
    X1_mat[l, ]  <- X1_mat[h, ]
    X2_mat[l, ]  <- X2_mat[h, ]
  }

  alpha_J  <- alpha[donor_idx]

  w_full <- numeric(N); w_full[donor_idx] <- w_star_J
  alpha_tr  <- sum(w_full[donor_idx] * alpha_J)
  X1_mat[treated_idx, ] <- as.numeric(crossprod(w_full, X1_mat))
  X2_mat[treated_idx, ] <- as.numeric(crossprod(w_full, X2_mat))

  X1_tr     <- X1_mat[treated_idx, ]
  X2_tr     <- X2_mat[treated_idx, ]

  Y_all <- numeric(N * T); D_all <- numeric(N * T); Y0_all <- numeric(N * T)
  tau_seq       <- rnorm(T - T_0, mean = delta_mean, sd = delta_sd)
  true_att_path <- numeric(0)

  for (t in seq_len(T)) {
    u_t  <- rnorm(J, 0, sigma_u)

    b_t  <- alpha_J + X1_mat[donor_idx, t] * beta[1] +
      X2_mat[donor_idx, t] * beta[2] + u_t

    if (type %in% c("SDM", "SLX")) {
      X_t  <- cbind(X1_mat[, t], X2_mat[, t])
      WX_t <- (W %*% X_t)[donor_idx, , drop = FALSE]
      b_t  <- b_t + as.numeric(WX_t %*% theta)
    }

    Y_tr_0 <- alpha_tr + X1_tr[t] * beta[1] + X2_tr[t] * beta[2]

    D_t <- rep(0, N); Y_t <- numeric(N); Y0_t <- numeric(N)

    if (t <= T_0) {
      Y_J <- as.numeric(M %*% (rho * w_1 * Y_tr_0 + b_t))
      Y_t[treated_idx]  <- Y_tr_0; Y_t[donor_idx]  <- Y_J
      Y0_t[treated_idx] <- Y_tr_0; Y0_t[donor_idx] <- Y_J
    } else {
      tau_t  <- tau_seq[t - T_0]
      Y_tr_1 <- Y_tr_0 + tau_t
      Y_J_1  <- as.numeric(M %*% (rho * w_1 * Y_tr_1 + b_t))
      Y_J_0  <- as.numeric(M %*% (rho * w_1 * Y_tr_0 + b_t))

      D_t[treated_idx]  <- 1
      Y_t[treated_idx]  <- Y_tr_1; Y_t[donor_idx]  <- Y_J_1
      Y0_t[treated_idx] <- Y_tr_0; Y0_t[donor_idx] <- Y_J_0
      true_att_path <- c(true_att_path, tau_t)
    }

    rng <- ((t - 1) * N + 1):(t * N)
    Y_all[rng] <- Y_t; D_all[rng] <- D_t; Y0_all[rng] <- Y0_t
  }

  df <- data.frame(time = rep(1:T, each = N), id = rep(1:N, times = T),
                   D = D_all, Y = Y_all, Y0 = Y0_all,
                   X1 = as.vector(X1_mat), X2 = as.vector(X2_mat))

  contam <- if (rho != 0) contamination_vector(W, treated_idx, rho) else NULL

  twin_df <- data.frame(
    pair          = seq_len(n_twins),
    contaminated  = high_idx,
    twin          = low_idx,
    s_abs_contam  = s_abs_J[high_pos],
    s_abs_twin    = s_abs_J[low_pos],
    target        = twin_target
  )

  list(
    df             = df,
    W              = W,
    alpha          = alpha,
    true_att       = true_att_path,
    contam         = contam,
    w_star         = w_full,
    w_star_donor   = planted$w_star_donor,
    planted_donors = planted$planted_donors,
    weight_profile = planted$weight_profile,
    donor_idx      = donor_idx,
    s_abs_picked   = planted$s_abs_picked,
    twin_tab       = twin_df,
    twins          = twin_df,
    settings       = list(N = N, T = T, T_0 = T_0, rho = rho,
                          sigma_u = sigma_u, delta_mean = delta_mean,
                          delta_sd = delta_sd, X_mean = X_mean,
                          x_sd = x_sd, alpha_sd = alpha_sd, n_twins = n_twins,
                          twin_gap = twin_gap, twin_target = twin_target,
                          weight_profile = weight_profile)
  )
}
