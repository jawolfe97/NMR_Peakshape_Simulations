# --- Header: Clear workspace, plots, and set working directory ---
{
  rm(list = ls())
  if (!is.null(dev.list())) dev.off()
  
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    current_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(current_path)) {
      setwd(dirname(current_path))
      message("Working directory set to: ", getwd())
    } else {
      warning("Cannot determine script path; working directory not changed.")
    }
  } else {
    warning("Package 'rstudioapi' not installed; working directory not changed.")
  }
}

fftshift <- function(x) {
  n <- length(x)
  p <- ceiling(n/2)
  c(x[(p+1):n], x[1:p])
}

SEQ <- seq(1, 5, by = 0.25)

MHz <- 750
ppm_center <- 8.0
delta_ppm <- 0.5
R2 <- (1/0.112)
n_points <- 2048
sw <- 4000
dt <- 1 / sw
t <- seq(0, (n_points - 1)) * dt

# Precompute frequency axis once (ppm scale)
freq_axis_Hz <- seq(-sw/2, sw/2, length.out = n_points)
freq_axis_ppm <- ppm_center + freq_axis_Hz / MHz

# Initialize dataframe with ppm axis as first column for reference
df <- data.frame(ppm = freq_axis_ppm)

for (n in 1:length(SEQ)) {
  
  kex <- 10^SEQ[n]
  delta_Hz <- delta_ppm * MHz
  omega_A <- (-delta_Hz / 2) * 2 * pi
  omega_B <- (delta_Hz / 2) * 2 * pi
  kAB <- kex / 2
  kBA <- kex / 2
  
  K <- matrix(c(-R2 - 1i * omega_A - kAB, kBA,
                kAB, -R2 - 1i * omega_B - kBA), nrow = 2, byrow = TRUE)
  
  pA <- kBA / (kAB + kBA)
  pB <- kAB / (kAB + kBA)
  M0 <- c(pA + 0i, pB + 0i)
  
  eigen_K <- eigen(K)
  
  exp_diag <- sapply(t, function(time) {
    exp(eigen_K$values * time)
  })
  
  V <- eigen_K$vectors
  V_inv <- solve(V)
  
  M_t <- matrix(0+0i, nrow = n_points, ncol = 2)
  for (i in seq_len(n_points)) {
    M_t[i, ] <- V %*% (exp_diag[, i] * (V_inv %*% M0))
  }
  
  signal <- rowSums(M_t)
  
  spectrum <- fft(signal)
  spectrum <- fftshift(spectrum)
  
  # Add real part of spectrum to dataframe with dynamic column name
  col_name <- paste0("kex = ", SEQ[n])
  df[[col_name]] <- Re(spectrum)
  
  # Optional: plot each spectrum (can comment out if not needed)
  plot(freq_axis_ppm, Re(spectrum), type = "l",
       xlab = "ppm", ylab = "Intensity",
       main = paste0("Two-State Bloch-McConnell NMR Spectrum\nkex = 10^", SEQ[n]),
       xlim = c(ppm_center + delta_ppm, ppm_center - delta_ppm))  # reversed ppm scale
}

# Save dataframe to CSV
write.csv(df, "Output.csv", row.names = FALSE)
