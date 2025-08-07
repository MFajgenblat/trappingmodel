# Code for the temporal disaggregation approach through interval-integrated B-splines

This repository contains the code for the statistical analysis presented in the manuscript entitled "Temporal disaggregation through interval-integrated B-splines for the integrated analysis of trapping counts in ecology."

The main folder contains four R scripts that are used throuhgout the paper. Additionally, the [Models](Models/) folder includes the four Stan models used and mentioned in the paper. The [Data](Data/) folder contains the necessary data to replicate the case study. The pitfall trapping dataset has been desensitized by excluding personal names, by adding a normally distributed jitter (sd = 5 km) to the exact trapping locations, by randomly discarding 90% of the trapping events, and by only considering the 50 most abundantly trapped species.

The folders [Simulation](Simulation/), [Output](Output/) and [Plots](Plots/) are empty but get populated while running the R scripts.
