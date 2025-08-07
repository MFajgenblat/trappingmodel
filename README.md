# Code for the temporal disaggregation approach through interval-integrated B-splines

This repository contains the code for the statistical analysis presented in the manuscript entitled "Temporal disaggregation through interval-integrated B-splines for the integrated analysis of trapping counts in ecology."

The main folder contains four R scripts that are used throughout the paper:
- [Script_A_simulation_study.R](Script_A_simulation_study.R) is needed to replicate the simulation study
- [Script_B_data_preparation.R](Script_B_data_preparation.R) is needed to prepare the data for the case studies
- [Script_C_varying_phenology.R](Script_C_varying_phenology.R) is needed to replicate the case study on inferring yearly phenological patterns and change thereof
- [Script_D_JSDM.R](Script_D_JSDM.R) is needed to replicate the cases study on estimating a Joint Species Distribution Model (JSDM) accounting for temporal aggregation

Additionally, the [Models](Models/) folder includes the four Stan models used and mentioned in the paper.

he [Data](Data/) folder contains the necessary data to replicate the case study. The pitfall trapping dataset has been desensitized by excluding personal names, by adding a normally distributed jitter (sd = 5 km) to the exact trapping locations, by randomly discarding 90% of the trapping events, and by only considering the 50 most abundantly trapped species.

The folders [Simulation](Simulation/), [Output](Output/) and [Plots](Plots/) are empty but get populated while running the R scripts.
