# --- Clear Workspace ---
rm(list = ls())
if (!is.null(dev.list())) dev.off()

# --- Function Definitions ---
fftshift <- function(x) {
  n <- length(x)
  p <- ceiling(n/2)
  c(x[(p+1):n], x[1:p])
}

# --- Adjustable Parameters ---
log_kex_values <- seq(2, 4, by = 0.5)   # Exchange rates (log10 scale)
pA_values <- seq(0, 1, by = 0.1)        # Populations of state A
delta_ppm <- 0.5                        # Chemical shift difference (ppm)
res_num <- 200                          # Number of residues
ppm_center <- 8.0                       # Spectral center (ppm)
MHz <- 700                              # Spectrometer frequency (MHz)

# --- Derived Parameters ---
protein_kDa <- res_num * 110 / 1000
R2 <- 6 + 0.3 * protein_kDa

# --- Simulation Constants ---
n_points <- 2048
sw <- 4000
dt <- 1 / sw
t <- seq(0, (n_points - 1)) * dt

# --- Frequency Axis ---
freq_axis_Hz <- seq(-sw/2, sw/2, length.out = n_points)
freq_axis_ppm <- ppm_center + freq_axis_Hz / MHz
delta_Hz <- delta_ppm * MHz
omega_A <- (-delta_Hz / 2) * 2 * pi
omega_B <- (delta_Hz / 2) * 2 * pi

# --- Set Layout: 5 rows x 11 columns ---
par(mfrow = c(5, 11), mar = c(1.5, 1.5, 1.5, 0.5), oma = c(4, 4, 4, 1))

# --- Find y-axis limits across all spectra to standardize
get_spec_range <- function(log_kex, pA) {
  kex <- 10^log_kex
  pB <- 1 - pA
  kAB <- kex * pB
  kBA <- kex * pA
  
  K <- matrix(c(-R2 - 1i * omega_A - kAB, kBA,
                kAB, -R2 - 1i * omega_B - kBA), nrow = 2, byrow = TRUE)
  M0 <- c(pA + 0i, pB + 0i)
  eigen_K <- eigen(K)
  V <- eigen_K$vectors
  V_inv <- solve(V)
  
  exp_diag <- sapply(t, function(time) exp(eigen_K$values * time))
  M_t <- matrix(0+0i, nrow = n_points, ncol = 2)
  for (i in seq_len(n_points)) {
    M_t[i, ] <- V %*% (exp_diag[, i] * (V_inv %*% M0))
  }
  
  signal <- rowSums(M_t)
  spec_real <- Re(fftshift(fft(signal)))
  range(spec_real)
}

# --- Precompute global y-axis limits
all_ranges <- unlist(lapply(log_kex_values, function(lk)
  lapply(pA_values, function(pa) get_spec_range(lk, pa))))
ylim_fixed <- range(all_ranges)

# --- Plot All Subplots ---
for (log_kex in log_kex_values) {
  for (pA in pA_values) {
    kex <- 10^log_kex
    pB <- 1 - pA
    kAB <- kex * pB
    kBA <- kex * pA
    
    K <- matrix(c(-R2 - 1i * omega_A - kAB, kBA,
                  kAB, -R2 - 1i * omega_B - kBA), nrow = 2, byrow = TRUE)
    M0 <- c(pA + 0i, pB + 0i)
    eigen_K <- eigen(K)
    V <- eigen_K$vectors
    V_inv <- solve(V)
    
    exp_diag <- sapply(t, function(time) exp(eigen_K$values * time))
    M_t <- matrix(0+0i, nrow = n_points, ncol = 2)
    for (i in seq_len(n_points)) {
      M_t[i, ] <- V %*% (exp_diag[, i] * (V_inv %*% M0))
    }
    
    signal <- rowSums(M_t)
    spec_real <- Re(fftshift(fft(signal)))
    
    plot(freq_axis_ppm, spec_real, type = "l", col = "blue",
         axes = FALSE, xlab = "", ylab = "",
         main = paste0("log??????(kex) = ", log_kex, "\n[A] = ", round(pA, 1)),
         cex.main = 0.6,
         xlim = c(ppm_center + delta_ppm, ppm_center - delta_ppm),
         ylim = ylim_fixed)
    box()
  }
}

# --- Add outer axis labels ---
mtext("Varying pA\nShift (ppm)", side = 1, outer = TRUE, line = 2)
mtext("Real Spectrum (A.U)\n Varying kex", side = 2, outer = TRUE, line = 1)
mtext("Bloch-McConnell Simulation Grid", side = 3, outer = TRUE, line = 1, cex = 1.4)
