# Multiscale Biomass Pyrolysis Modeling
Developers: Meagan Crowley, Hari Sitaraman, Hilary Egan, Peter Ciesielski 
contact: mcrowley2@nrel.gov

This repository contains multiscale reactive computational fluid dynamics models of biomass pyrolysis in the 2-fluidized bed reactor (2FBR) at NREL developed in OpenFOAM and COMSOL. 

The following sub-directories are included:
* reactor_model
- This directory contains the 2D model of the 2FBR for biomass pyrolysis developed in OpenFOAM v11.
- The Miller-Bellan pyrolysis scheme was used to model the conversion of biomass to biofuel and side-products via fast-pyrolysis as reported in Xiong, Kong, and Passalacqua, Chem. Eng. Sci. 99 2013 305-313 (http://dx.doi.org/10.1016/j.ces.2013.06.017)

* particle_model
- This directory contains the 2D model of a pine particle under reactor conditions of the 2FBR developed in COMSOL Multiphysics v 6.1.
-


Acknowledgements: The authors would like to thank Brennan Pecha, who helped develop early iterations of the particle models in COMSOL, Calvin Mukarakate, Kristiina Iisa, and others on their team in the Catalytic Carbon Transformation & Scale-Up Center at NREL who collected validation data from the 2FBR , NREL's high performance computing (HPC) kestrel super computer which provided computational resources for openFOAM simulations, the Feedstock Conversion Interface Consortium (FCIC) and the United States Department of Energy (DOE) Bioenergy Technologies Office (BETO) who provided funding for this research.