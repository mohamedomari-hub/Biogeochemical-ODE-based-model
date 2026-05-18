# Biogeochemical Model

This repository contains an R/C implementation of a dynamic biogeochemical model
for aquatic ecosystem processes. The model simulates carbon, nitrogen,
phosphorus, dissolved oxygen, algae, and submerged aquatic plant dynamics using
environmental forcing data.

## Overview

The model is implemented as ordinary differential equations compiled from C and
called from R through `deSolve`. The main workflow:

1. Load model parameters from `ModelFinal_para.R`.
2. Compile and load the C model implementation in `ModelFinal.c`.
3. Read environmental forcing data from `data_forcing3.csv`.
4. Run the ODE model with `lsodes`.
5. Plot simulated state trajectories.
6. Optionally run local sensitivity analysis with `FME`.

## Model States

The model includes the following state variables:

```text
Algae, DWeed, DO,
POC, DOC, DIC, TC,
ON, NH4_N, NO3_N, TN,
OP, PO4, TP,
TPOC, TON, TOP
```

These represent phytoplankton/algae, duckweed or submerged plant biomass,
dissolved oxygen, organic/inorganic carbon pools, nitrogen pools, and phosphorus
pools.

## Repository Contents

```text
ModelFinal.c              C implementation of the ODE right-hand side
ModelFinal_R.c            R-generated C interface/helper code
ModelFinal_para.R         Model parameter definitions
ModelFinal_Control.R      Main model run script
Model_Sensitivity.R       Sensitivity analysis workflow
data_forcing3.csv         Environmental forcing and observation data
Published_Paper.pdf       Related published paper
Omari_paperdraft.pdf      Draft manuscript
```

Compiled files such as `ModelFinal.o` and `ModelFinal.so` can be regenerated
from `ModelFinal.c`.

## Requirements

R packages:

```r
install.packages(c(
  "deSolve",
  "zoo",
  "xlsx",
  "Hmisc",
  "scatterplot3d",
  "FME",
  "reshape2"
))
```

You also need a working C compiler available to R, because the model is compiled
with:

```r
R CMD SHLIB ModelFinal.c
```

## Running the Model

From the project directory, open R and run:

```r
source("ModelFinal_Control.R")
```

The control script compiles the C model, loads the shared library, reads the
forcing data, runs the ODE simulation, and plots the model states.

If needed, update the hard-coded paths in `ModelFinal_Control.R` from the
original local path to the current project directory.

## Running Sensitivity Analysis

Run:

```r
source("Model_Sensitivity.R")
```

This uses the `FME` package to compute local sensitivity functions for model
states with respect to model parameters.

## Forcing Data

`data_forcing3.csv` contains environmental drivers and observations, including:

```text
DateTime, radiation, water temperature, depth, air temperature,
wind speed, PO4, NO3, NH4, DOC, DIC, algae, duckweed, DO
```

The file uses semicolon separators and comma decimal notation.

## Notes

- The current scripts contain absolute paths from the original development
  machine. Before running on another computer, replace those paths with relative
  paths.
- Generated files are excluded through `.gitignore`.
- The model equations are in `ModelFinal.c`; parameters are in
  `ModelFinal_para.R`.

## Citation

If this model is used in research, cite the related publication included in this repository as `Published_Paper.pdf`.
