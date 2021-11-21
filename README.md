# SulfateSpectraJBO

## Overview
 Upload optical absorption spectrum data of CuSO4 and NiSO4 solutions at various dilutions and mixtures (As reported in https://doi.org/10.1117/1.JBO.22.12.125007)


## Spectra

<b>Figure. Optical absorption of solutions of either NiSO4 and CuSO4 at various levels of dilution. (a, c) Absorption spectra of four levels of dilution (error bars omitted from the figures as they could not be properly resolved—the standard deviation of three measurements per solution was below 1% of the mean for most cases and wavelengths); (b, d) plotting optical absorption versus concentration at a prominent wavelength to highlight that absorption linearity with concentration is respected. Three measurements are shown per concentration, but their overlap is such that they cannot be easily resolved when plotted.</b>
![Figure1](plots/JBO_22_12_125007_f001.png "Optical absorption of solutions of either NiSO4 and CuSO4 at various levels of dilution. (a, c) Absorption spectra of four levels of dilution (error bars omitted from the figures as they could not be properly resolved—the standard deviation of three measurements per solution was below 1% of the mean for most cases and wavelengths); (b, d) plotting optical absorption versus concentration at a prominent wavelength to highlight that absorption linearity with concentration is respected. Three measurements are shown per concentration, but their overlap is such that they cannot be easily resolved when plotted.")

<b>Figure. Absorption spectra of mother batch solutions of NiSO4 and CuSO4 and of three mixtures of these, at ratios 1:3, 2:2, and 3:1. The mixtures have a measured absorption spectrum (meas) that matches the predicted spectrum (pred). The predictions are made through simple linear combination of the measured spectra of the mother batch solutions of nickel sulfate, cNiSO4,b ([1:0]), and copper sulfate, cCuSO4,b ([0:1]), weighted by the expected ratios.</b>
![Figure1](plots/JBO_22_12_125007_f002.png "Absorption spectra of mother batch solutions of NiSO4 and CuSO4 and of three mixtures of these, at ratios 1:3, 2:2, and 3:1. The mixtures have a measured absorption spectrum (meas) that matches the predicted spectrum (pred). The predictions are made through simple linear combination of the measured spectra of the mother batch solutions of nickel sulfate, cNiSO4,b ([1:0]), and copper sulfate, cCuSO4,b ([0:1]), weighted by the expected ratios.")



## Underlying data
Data in Matlab .mat structure file:
- wl - vector of wavelengths
- uaA - matrix of mu_a [mm-1] of each measurement taken (33, 3 per physical sample/solution) - 33 cols, each row a wavelength
- uaAav - matrix of mu_a [mm-1] but takes simple average of the three measurements on the same solution - 11 cols, each row a wavelength
- sampleid - identifies the 33 measurements into the 11 solutions. `sampleid=[1*ones(1,3) 2*ones(1,3) 3*ones(1,3) 4*ones(1,3) 5*ones(1,3) 6*ones(1,3) 7*ones(1,3) 8*ones(1,3) 9*ones(1,3) 10*ones(1,3) 11*ones(1,3)]`
- conccuav - an 11 item vector giving the CuSO4.5H2O concentration in Molar for that sample id -e.g. sample id one has 1 molar CuSO4.5H2O but sample id two has nil
- concniav - an 11 item vector giving the NiSO4.6H2O concentration in Molar for that sample id -e.g. sample id one has nil of NiSO4.6H2O but sample id two has 2.2M .

The Matlab analysis script 'analysis_ft.m' is provided, which generated most optical absorption related analysis and figures in the published study. This is for guidance - it will not run end-to-end as the upstream raw data is not uploaded. The 'data_sulphates.mat' csv be loaded instead.
Note that for some figures, Inkscape was used for post-processing of SVG figures.

'Hbspec.mat' contains oxy and deoxyhemoglobin spectra used downstream in the script for illustrative purposes.
Data sourced from Scott Prahl's ['Tabulated Molar Extinction COefficient for Hemoglobin in Water' OMLC page](https://omlc.org/spectra/hemoglobin/summary.html). (tbc)
