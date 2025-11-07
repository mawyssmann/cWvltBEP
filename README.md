Author: Micah A. Wyssmann

This wavelet-based software in the folder named after the version (such as "V1.0.1") contains: 
- A class with all of the processing functions (cWvltBEP.m)
- A driver script that is used to call the function class (dr_cWvltBEP_BASIC.m). 

The recommended usage is to go to the driver script and call it in sections. This driver script is constructed specifically to import and process the example (synthetic) data that is also included in the repository. If you have your own data, you can load it into the workspace and assign them with the script. Alternatively, the driver script lays out the general workflow of calculations if you would like to write your own functions/scripts to call datasets (i.e., import data from folders, loop over datasets, etc.).

Also contained herein is the functions used for generating synthetic BEPs with dunes of known scales (in folder "SimDuneField"). This folder contains the function "simDuneField.m" that takes inputs with a vector of desired dune parameters and generates the corresponding synthetic BEP. The driver script "dr_SimDuneField_MO.m" synthesizes a BEP to mimic the Missouri (MO) River data. To run this script, you also need to load "ARIMAmdl_MORivJune_p8q6.mat" that contains the ARMA model fitted to the MO River data. 
