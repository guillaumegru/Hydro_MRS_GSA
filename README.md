# Sensitivity of time-lapse magnetic resonance sounding to vadose zone hydrodynamic parameters: monitoring of an intense meteorological event

This repository contains the supporting codes and data for the paper **Sensitivity of time-lapse magnetic resonance sounding to vadose zone hydrodynamic parameters: monitoring of an intense meteorological event** by Guillaume Gru, Jean-François Girard, Philippe Ackerer and Nolwenn Lesparre.

## Description

Folders list : 

* Hydro_MRS_GSA_numexp: routines and main programs for sampling the hydrodynamic parameters space and running Phydrus-MRS model simulations.
* Hydro_MRS_GSA_numexp_postprocess: jupyter notebooks used to train the PCE metamodels and compute Sobol indices, compute AMA sensitivity indices, and visualize the results to produce the figures (7, 8, 9 and 10) of the paper.
* Phydrus_MRS_modules-py: routines to compute MRS signals from water content distributions as described in the Background section of the paper.
* Plots_hydrometeo: processed data for hydrological model boundary conditions and jupyter notebook for visualization (Figure 2).
* Print_MRS_kernel: MRS kernel computed with Samovar and jupyter notebook for visualization (Figure 4).
* SA_module: routines to compute AMA sensitivity indices.


## Authors

Guillaume Gru (gru@unistra.fr)
Jean-François Girard
Philippe Ackerer
Nolwenn Lesparre (lesparre@unistra.fr)

## Help

To reduce the size of the folders to upload, we compressed some outputs to zip files. In order for the postprocessing codes to run properly, one needs to extract the files in /Hydro_MRS_GSA_numexp_postprocess/Res_sim_1.zip and Res_sim_2.zip and place them into a folder named /Hydro_MRS_GSA_numexp_postprocess/Res_sim. Likewise, the pkl files storing the data relative to the PCE were zipped into 3 files, which need to be extracted to the folder \Hydro_MRS_GSA_numexp_postprocess\Data for the postprocessing notebook to run succesfully.


## Version History

* 1.0
    * (01/28/2026) Initial Release: for submission to WRR



