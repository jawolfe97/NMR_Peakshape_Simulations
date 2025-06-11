#Clear
{rm(list = ls())
if (!is.null(dev.list())) dev.off()

}

library(expm)
library(pracma)

# --- Helper: FFT Shift ---
fftshift <- function(x) {
  n <- length(x)
  p <- ceiling(n / 2)
  c(x[(p + 1):n], x[1:p])
}

# --- Bloch-McConnell Simulation ---
simulate_bloch_mcconnell <- function(kex, dw_Hz, pA, R2, t) {
  omega_A <- -dw_Hz / 2 * 2 * pi
  omega_B <-  dw_Hz / 2 * 2 * pi
  
  K <- matrix(c(
    -R2 - kex * (1 - pA),  kex * (1 - pA),
    kex * pA,            -R2 - kex * pA
  ), 2, 2, byrow = TRUE)
  
  M0 <- c(pA, 1 - pA)
  signal <- complex(length(t))
  
  for (i in seq_along(t)) {
    Mt <- expm(K * t[i]) %*% M0
    signal[i] <- Mt[1] * exp(1i * omega_A * t[i]) + Mt[2] * exp(1i * omega_B * t[i])
  }
  
  spectrum <- fftshift(Mod(fft(signal)))
  return(spectrum)
}

# --- Frequency Axis in ppm ---
frequency_axis_ppm <- function(n_points, sw, MHz) {
  df_Hz <- sw / n_points
  freqs_Hz <- seq(-sw / 2, sw / 2 - df_Hz, length.out = n_points)
  freqs_ppm <- freqs_Hz / MHz
  return(freqs_ppm)
}

# --- Estimate kex via Least-Squares Over Full Spectrum ---
estimate_kex <- function(observed_spectrum, x01, x02, R2,
                         pA, MHz = 700, n_points = 2048, sw = 4000,
                         kex_bounds = c(1, 10000), plot = TRUE) {
  dt <- 1 / sw
  t <- seq(0, (n_points - 1)) * dt
  dw_Hz <- abs(x02 - x01) * MHz
  
  objective_function <- function(kex) {
    simulated <- simulate_bloch_mcconnell(kex, dw_Hz, pA, R2, t)
    sum((simulated - observed_spectrum)^2)
  }
  
  fit <- tryCatch(
    optimize(objective_function, interval = kex_bounds),
    error = function(e) NA
  )
  
  if (is.na(fit)) {
    warning("Optimization failed.")
    return(NA)
  }
  
  kex_est <- fit$minimum
  
  # Optional plotting
  if (plot) {
    sim_fit <- simulate_bloch_mcconnell(kex_est, dw_Hz, pA, R2, t)
    ppm_axis <- frequency_axis_ppm(n_points, sw, MHz)
    
    plot(ppm_axis, observed_spectrum, type = "l", col = "black", lwd = 2,
         xlab = "ppm", ylab = "Signal Intensity", xlim = c(-0.5,0.75) , main = sprintf("Fit: True vs Estimated\nEstimated kex = %.0f | True kex = 3000", kex_est)
         )
    lines(ppm_axis, sim_fit, col = "red", lwd = 2, lty = 2)
    legend("topright", legend = c("Simulated (True)", "Fitted (Estimated)"), col = c("black", "red"), lwd = 2, lty = c(1, 2))
  }
  
  return(kex_est)
}

# --- Wrapper: Simulate and Estimate ---
simulate_and_estimate <- function(kex_true, dw_ppm = 0.5, pA = 0.5, R2 = 10,
                                  MHz = 700, n_points = 2048, sw = 4000) {
  dt <- 1 / sw
  t <- seq(0, (n_points - 1)) * dt
  dw_Hz <- dw_ppm * MHz
  
  # Simulate ground truth spectrum
  spectrum <- simulate_bloch_mcconnell(kex_true, dw_Hz, pA, R2, t)
  
  # Simulated peak centers
  x01 <- -dw_ppm / 2
  x02 <-  dw_ppm / 2
  
  # Estimate kex from spectrum
  kex_est <- estimate_kex(spectrum, x01, x02, R2, pA, MHz, n_points, sw, plot = TRUE)
  
  cat(sprintf("True kex: %g | Estimated kex: %.2f\n", kex_true, kex_est))
  return(kex_est)
}

# --- Example Usage ---
set.seed(1)
simulate_and_estimate(kex_true = 3000)
