# MaxEnt species distribution model

This folder contains the inputs, configuration details, and outputs of the MaxEnt species distribution model used in the analysis.

## Software

The model was run using **MaxEnt v3.4.4** (https://github.com/mrmaxent/Maxent).

## Input data

The model used the following inputs:

* **Occurrence data**: large_mammals_thinned.csv 
* **Environmental predictors**: ASCII raster layers (`.asc`)

## Model settings

The log of the model run is provided in the repository and indicates all the selected parameters.

## Outputs

The model produced:

* Raster prediction: `Large_mammal.asc`
* Response plots: `plots/Large_mammal.png`

Outputs were written to the MaxEnt output directory during execution and are included here for reproducibility.

## Reproducibility note

Users should download independently the environmental predictors.

The file paths shown in the command line and log correspond to the original execution environment. 
Users should update all paths to match their local directory structure before running the model.
