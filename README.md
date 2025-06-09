# NMR_Peakshape_Simulations
**NMR_SimulationMatrix_kex_vs_pa.R**
![Description](https://github.com/jawolfe97/NMR_Peakshape_Simulations/blob/main/pA%20vs%20kex.svg)
 ___________________________________________________________________________
```r
--- Adjustable Parameters ---
log_kex_values <- seq(2, 4, by = 0.5)   # Exchange rates (log10 scale)
pA_values <- seq(0, 1, by = 0.1)        # Populations of state A
delta_ppm <- 0.5                        # Chemical shift difference (ppm)
res_num <- 200                          # Number of residues
ppm_center <- 8.0                       # Spectral center (ppm)
MHz <- 700                              # Spectrometer frequency (MHz)
```
___________________________________________________________________________
_**NMR_SingleSimulation_Visualizer.R**_
![Description](https://github.com/jawolfe97/NMR_Peakshape_Simulations/blob/main/Single_Simulation.svg)
___________________________________________________________________________ 
```r
--- "Sliders": Adjustable Parameters ---
log_kex   <- 2       # Exchange rate (log10 scale in s???¹)
delta_ppm <- 0.5       # Chemical shift difference (ppm), usually < 2 ppm
res_num   <- 200      # Number of residues (e.g. 363 residues ~ 40 kDa protein)
ppm_center <- 8.0      # Spectral center (ppm)
MHz       <- 700       # Spectrometer frequency (MHz)
pA        <- 0.5       # Population of State A (0 = all B, 1 = all A)
```
