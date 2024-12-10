# 2D Pine particle pyrolysis CFD model developed in COMSOL
Developers: Meagan Crowley, Peter Ciesielski 

contact: mcrowley2@nrel.gov

- This directory contains the 2D model of a pine particle under reactor conditions of the 2FBR developed in COMSOL Multiphysics v 6.1.
- The Miller-Bellan pyrolysis scheme was also used in this model (see citation above).
- Pine and char transport properties were determined from X-ray computed tomography and microstructural analysis (see: Crowley et al, Front. Ener. Res. 10, 2022. https://doi.org/10.3389/fenrg.2022.850630)


## System requirements
- This model was developed in COMSOL v.6.1 and requires a COMSOL installation and software license to run. Additionally, the model uses the following COMSOL modules:
    - Chemical Reaction Engineering Module
    - Heat Transfer Module
- This model has been tested on a workstation with an Intel Xeon 81600 CPU with 24 cores. It takes about 3 hours to run per particle size.



## How to run
The model can be run from the COMSOL GUI using the parametric sweep study which runs fast-pyrolysis simulations for five different pine particle sizes.

<img src="./bigparticle_temperature.gif"/>
