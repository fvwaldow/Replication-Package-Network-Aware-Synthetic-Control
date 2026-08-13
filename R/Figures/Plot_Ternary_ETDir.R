# Ternary plots for the network examples

A_1 <- matrix(c(
  0, 1, 1, 1,
  1, 0, 1, 1,
  1, 1, 0, 1,
  1, 1, 1, 0
), nrow = 4, byrow = TRUE)
A_2 <- matrix(c(
  0, 0, 1, 1,
  0, 0, 0, 1,
  1, 0, 0, 1,
  1, 1, 1, 0
), nrow = 4, byrow = TRUE)
A_3 <- matrix(c(
  0, 0, 1, 1,
  0, 0, 0, 0,
  1, 0, 0, 1,
  1, 0, 1, 0
), nrow = 4, byrow = TRUE)
A_4 <- matrix(c(
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 1,
  0, 0, 1, 0
), nrow = 4, byrow = TRUE)

A_list <- list(A_1, A_2, A_3, A_4)

rho    <- 0.2
lambda <- 10
alpha  <- c(1.2, 1.2, 1.2)

log_dens_etd <- function(x1, x2, alpha, theta) {
  x3 <- 1 - x1 - x2
  if (x1 <= 0 || x2 <= 0 || x3 <= 0) return(-Inf)
  sum((alpha - 1) * log(c(x1, x2, x3))) + sum(theta * c(x1, x2, x3))
}

v1 <- c(0, 0)             # x = (1,0,0)
v2 <- c(1, 0)             # x = (0,1,0)
v3 <- c(0.5, sqrt(3)/2)   # x = (0,0,1)


nx <- 1000; ny <- 1000
xg <- seq(0, 1, length.out = nx)
yg <- seq(0, sqrt(3)/2, length.out = ny)

xy_to_bary <- function(X, Y) {
  denom <- (v2[2] - v3[2]) * (v1[1] - v3[1]) + (v3[1] - v2[1]) * (v1[2] - v3[2])
  l1 <- ((v2[2] - v3[2]) * (X - v3[1]) + (v3[1] - v2[1]) * (Y - v3[2])) / denom
  l2 <- ((v3[2] - v1[2]) * (X - v3[1]) + (v1[1] - v3[1]) * (Y - v3[2])) / denom
  l3 <- 1 - l1 - l2
  list(x1 = l1, x2 = l2, x3 = l3)
}

make_Z <- function(theta) {
  Z <- matrix(NA, nx, ny)
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      b <- xy_to_bary(xg[i], yg[j])
      if (b$x1 > 0 && b$x2 > 0 && b$x3 > 0) {
        Z[i, j] <- exp(log_dens_etd(b$x1, b$x2, alpha, theta))
      }
    }
  }
  dxdy <- (xg[2] - xg[1]) * (yg[2] - yg[1])
  Z / sum(Z, na.rm = TRUE) / dxdy
}

draw_panel <- function(Z) {
  plot(NA, xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.0),
       asp = 1, axes = FALSE, xlab = "", ylab = "")
  
  n_colors <- 100
  colors   <- hcl.colors(n_colors, palette = "YlOrRd", rev = TRUE)
  image(xg, yg, Z, col = colors, add = TRUE, useRaster = TRUE)
  
  zmax <- max(Z, na.rm = TRUE)
  zmin <- min(Z, na.rm = TRUE)
  if (is.finite(zmax) && zmax - zmin > 1e-10 * zmax) {
    levels <- seq(0.05 * zmax, 0.95 * zmax, length.out = 6)
    contour(xg, yg, Z,
            levels = levels, drawlabels = FALSE, col = "black", add = TRUE)
  }
  
  polygon(c(v1[1], v2[1], v3[1]), c(v1[2], v2[2], v3[2]),
          border = "black", lwd = 2)
  
  text(v1[1] - 0.00, v1[2] - 0.05, expression("(1,0,0)"), cex = 2.2)
  text(v2[1] + 0.00, v2[2] - 0.05,  expression("(0,1,0)"), cex = 2.2)
  text(v3[1],        v3[2] + 0.05, expression("(0,0,1)"), cex = 2.2)
}

theta_from_A <- function(A) {
  rs      <- rowSums(A)
  rs_safe <- ifelse(rs == 0, 1, rs)
  W       <- A / rs_safe
  
  W_J  <- W[2:4, 2:4]
  w_J1 <- W[2:4, 1]
  
  if (all(w_J1 == 0)) {
    s <- rep(0, 3)
  } else {
    s <- as.numeric(rho * solve(diag(3) - rho * W_J) %*% w_J1)
  }
  list(s = s, theta = -lambda * abs(s))
}

outdir <- dir_fig


for (i in seq_along(A_list)) {
  A   <- A_list[[i]]
  out <- theta_from_A(A)
  Z   <- make_Z(out$theta)
  pdf_path <- file.path(outdir, sprintf("etd_A%d.pdf", i))
  pdf(pdf_path, width = 6, height = 5.5)
  par(mar = c(0.2, 0.2, 0, 0.2))
  draw_panel(Z)
  dev.off()
}