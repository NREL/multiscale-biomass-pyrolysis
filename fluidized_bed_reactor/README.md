# 2D Fluidized bed pyrolysis reactor CFD model developed in OpenFOAM
Developers: Meagan Crowley, Hari Sitaraman, Hilary Egan, Peter Ciesielski 
contact: mcrowley2@nrel.gov

- This directory contains the 2D model of the 2FBR for biomass pyrolysis developed in OpenFOAM v11.
- The Miller-Bellan pyrolysis scheme was used to model the conversion of biomass to biofuel and side-products via fast-pyrolysis as reported in Xiong, Kong, and Passalacqua, Chem. Eng. Sci. 99, 2013, 305-313. (http://dx.doi.org/10.1016/j.ces.2013.06.017)
[![Bed volume fraction during the simulation of fast-pyrolysis in the 2FBR](https://raw.githubusercontent.com/NREL/multiscale-biomass-pyrolysis/tree/main/fluidized_bed_reactor/2FBR.gif)]

## System requirements
- This model was developed in OpenFOAM v11 and requires an OpenFOAM installation to run.
- An MPI library is required to run in parallel. 
- It has been tested four dual socket Intel Zeon Sapphire Rapids 52-core processors. It takes about 1.5 hours using 4 cores to run.

## How to run
The `run.sh` script contains commands used on NREL's Kestrel supercomputer to run the simulation in parallel using the multiphaseEulerFoam solver in OpenFOAM.

[![Bed volume fraction during the simulation of the 2FBR](https://raw.githubusercontent.com/nrel/multiscale-biomass-pyrolysis/fluidized_bed_reactor/2FBR.gif)]
