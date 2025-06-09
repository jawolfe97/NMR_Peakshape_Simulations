#Clear and Set Up WOrkspace
{
# --- Clear Workspace ---
rm(list = ls())
if (!is.null(dev.list())) dev.off()

# --- Function Definitions ---
fftshift <- function(x) {
  n <- length(x)
  p <- ceiling(n/2)
  c(x[(p+1):n], x[1:p])
}
}

## --- "Sliders": Adjustable Parameters ---
log_kex   <- 2       # Exchange rate (log10 scale in s???¹)
delta_ppm <- 0.5       # Chemical shift difference (ppm), usually < 2 ppm
res_num   <- 88      # Number of residues (e.g. 363 residues ~ 40 kDa protein)
ppm_center <- 8.0      # Spectral center (ppm)
MHz       <- 700       # Spectrometer frequency (MHz)
pA        <- 0.5       # Population of State A (0 = all B, 1 = all A)

#Derived Parameters and Simulation Constants
{
# --- Derived Parameters ---
protein_kDa <- res_num * 110 / 1000     # MW in kDa derived from residue count
pB  <- 1 - pA
R2  <- 6 + 0.3 * protein_kDa            # Empirical R2 for folded proteins (s???¹)

# --- Simulation Constants ---
n_points <- 2048
sw <- 4000
dt <- 1 / sw
t <- seq(0, (n_points - 1)) * dt
}
#Run and Visualize Simulation
{
# --- Frequency Axis ---
freq_axis_Hz <- seq(-sw/2, sw/2, length.out = n_points)
freq_axis_ppm <- ppm_center + freq_axis_Hz / MHz

# --- Bloch-McConnell Simulation ---
kex <- 10^log_kex
delta_Hz <- delta_ppm * MHz
omega_A <- (-delta_Hz / 2) * 2 * pi
omega_B <- (delta_Hz / 2) * 2 * pi
kAB <- kex * pB
kBA <- kex * pA

K <- matrix(c(-R2 - 1i * omega_A - kAB, kBA,
              kAB, -R2 - 1i * omega_B - kBA), nrow = 2, byrow = TRUE)

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

signal <- rowSums(M_t)                # complex FID
spectrum <- fftshift(fft(signal))     # FFT of FID

# --- Extract Components ---
fid_real <- Re(signal)
fid_imag <- Im(signal)
spec_real <- Re(spectrum)
spec_imag <- Im(spectrum)

# --- Exchange regime assessment ---
delta_omega <- 2 * pi * delta_Hz  # angular frequency diff in rad/s

{s
# # --- Define exchange regime based on ???? and kex ---
# exchange_regime <- if (kex < delta_omega / 10) {
#   "Slow exchange (kex << ????): Separate peaks are visible."
# } else if (kex > delta_omega * 10 && kex < 10 * R2) {
#   "Fast exchange (kex >> ????): Peaks coalesce to a single averaged signal."
# } else if (kex >= 10 * R2) {
#   "Too fast to observe (kex >> ???? and R2): Exchange is spectroscopically invisible."
# } else {
#   "Intermediate exchange (kex ??? ????): Broadening and coalescence effects observed."
# }
}

# Coalescence threshold derived from theoretical 2-site exchange
# ???? in Hz; coalescence occurs around kex ??? ?? * ???? / ???2 ??? 2.22 * ????
# ???? in rad/s = 2?? * ????
coalescence_kex <- pi * delta_Hz / sqrt(2)

exchange_regime <- if (kex < delta_omega / 10) {
  "Slow exchange (separate peaks)"
} else if (kex < coalescence_kex) {
  "Intermediate exchange (broad peaks, partial coalescence)"
} else if (kex < 10 * R2) {
  "Fast exchange (coalesced peak with minor shoulders)"
} else {
  "Very fast exchange (no exchange broadening observable)"
}


# --- Plot Layout ---
# Define layout matrix: 3 rows, 2 cols; last row spans both cols
layout_matrix <- matrix(c(1, 2,
                          3, 4,
                          5, 5), nrow = 3, byrow = TRUE)
layout(layout_matrix)

# Set margins for plots
par(mar = c(4, 4, 2, 1))

# Plot 1: Real FID
plot(t, fid_real, type = "l", col = "blue",
     xlab = "Time (s)", ylab = "FID Real",
     main = "Real Part of FID")

# Plot 2: Imaginary FID
plot(t, fid_imag, type = "l", col = "red",
     xlab = "Time (s)", ylab = "FID Imaginary",
     main = "Imaginary Part of FID")

# Plot 3: Real Spectrum
plot(freq_axis_ppm, spec_real, type = "l", col = "blue",
     xlab = "ppm", ylab = "Spectrum Real",
     main = "Real Spectrum",
     xlim = c(ppm_center + delta_ppm, ppm_center - delta_ppm))

# Plot 4: Imaginary Spectrum
plot(freq_axis_ppm, spec_imag, type = "l", col = "red",
     xlab = "ppm", ylab = "Spectrum Imaginary",
     main = "Imaginary Spectrum",
     xlim = c(ppm_center + delta_ppm, ppm_center - delta_ppm))

# Plot 5: Parameters and Exchange regime (spans bottom row)
plot.new()
text(0.5, 0.9, "Simulation Parameters:", cex = 1.5, font = 2)
text(0.5, 0.80, paste("Spectrometer Freq (MHz):", MHz), cex = 1.2)
text(0.5, 0.73, paste("???? (ppm):", delta_ppm), cex = 1.2)
text(0.5, 0.66, paste("???? (rad/s):", round(delta_omega, 1)), cex = 1.2)
text(0.5, 0.59, paste("Number of residues:", res_num), cex = 1.2)
text(0.5, 0.52, paste("MW (kDa):", round(protein_kDa, 2)), cex = 1.2)
text(0.5, 0.45, paste("R2 (s???¹):", round(R2, 2)), cex = 1.2)
text(0.5, 0.38, paste("Population A:", round(pA, 2)), cex = 1.2)
text(0.5, 0.31, paste0("kex (s???¹): 10^", round(log_kex, 2), " = ", signif(kex, 3)), cex = 1.2)
text(0.5, 0.24, paste("Exchange regime:", exchange_regime), cex = 1.2, font = 2)
}
