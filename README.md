# Demographic Resilience of Cerrado Lizards
This repository provides all the data and code (to run in R) from the analysis of the study "Severe fire regimes decrease resilience of ectothermic populations", published in Journal of Animal Ecology.

Authors: Heitor Campos de Sousa, Adriana Malvasio, Guarino Rinaldi Colli, and Rob Salguero-Gómez.

## Ecophysio
The folder "Ecophysio" provides the data and code necessary to download the weather data (ERA5), estimate the microclimate measures, and derive the ecophysiological variables.

## Hierarchical Models
The folder "Hierarchical_Models" provides the data and code necessary to estimate the vital rates of three Cerrado lizard species (_Copeoglossum nigropunctatum_, _Micrablepharus atticolus_, and _Tropidurus itambere_). We used Cormack Jolly-Seber (CJS) to relate survival with the individuals' body size (snout-vent length – SVL), and Pradel Jolly-Seber (PJS) models to estimate survival, recruitment, and derive population growth considering environmental variation. We also used Generalized Linear Models (GLMs) to relate probability (prep) of reproduction and number of newborns (nb, fecundity) with the individuals' body size (SVL). We provide the R scripts with the code for each species (thus, three R scripts).

## IPMs
The folder "IPMs" provides the data and code necessary to build stochastic Integral Projection Models (IPMs) for three Cerrado lizard species. The code also provides the perturbation analyses and the calculation of life history traits and demographic resilience components.

## DGAMs
The folder "DGAMs" provides the data and code necessary to run the comparative analyses between the three Cerrado lizard species, regarding their demographic resilience components and life history traits using dynamic generalized additive models (DGAMs). The code provides the code to generate all the figures present in the original manuscript, supplementary material, and others not used (such as diagnostic plots and other comparative plots). Besides DGAMs, the code also provides code for multilevel models using the package brms, used in previous versions of the manuscript that reinforces the robustness of the results.
