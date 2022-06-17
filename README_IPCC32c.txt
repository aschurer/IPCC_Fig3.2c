##########################################################################
# ---------------------------------------------------------------------------------------------------------------------
# This is a mix of Python code and shell scripts to produce IPCC AR6 WGI Figure 3.2c
# Creator: Andrew Schurer, University of Edinburgh
# Contact: a.schurer@ed.ac.uk
# Last updated on: April 26th, 2022
# --------------------------------------------------------------------------------------------------------------------
#
# - Code functionality: Create_plot.sh this is the main program and calls all other python scripts. To run the scripts requires the user to create a main directory which includes the directories: CMIP5, CMIP6, CESM_LME, PAGES2K, Volcanic, SavedModData. Once the directory structure is set-up and all the datasets are downloaded, create the figure with the command.
Create_plot.sh Path_to_main_directory. This calls: Process_CMIP5.sh: This script prepares the CMIP5 model files. Process_CMIP6.sh: This script prepares the CMIP6 model files. Process_CESM.py: This script prepares the CESM1 model files. Process_volcs.py: This script prepares the volcanic forcing files. Plot_past1000.py: This files reads and plots the pre-prepared data.
# - Input data: CCBY4.0 All the input data required is open access. All model data is distributed under the CCBY4.0 licence.
In each of the model directories (CMIP5, CMIP6, CESM_LME) the user must create a RAWDATA directory.
The RAWDATA directories need to be filled with the model data listed in MODELDATA.txt, . The PAGES2K directory must contain the 7 text reconstruction files downloaded from:
https://figshare.com/collections/Global_mean_temperature_reconstructions_over_the_Common_Era/4507043
BHM.txt, CPS.txt, DA.txt, M08.txt, OIE.txt, PAI.txt, PCR.txt
The Volcanic directory needs to contain three directories.
CMIP6, Gao, Crowley, and each must contain the relevant volcanic reconstruction. 
volcanic_sAOD_monthly_-50001-201912_new.csv From https://github.com/chrisroadmap/ar6/tree/main/data_output
IVI2LoadingLatHeight501-2000.txt From http://climate.envsci.rutgers.edu/IVI2/
ICI5_030N_AOD_c.txt, ICI5_030S_AOD_c.txt, ICI5_3090N_AOD_c.txt, ICI5_3090S_AOD_c.txt From https://wiki.lsce.ipsl.fr/pmip3/doku.php/pmip3:design:lm:final#volcanic_forcings 
# - Output variables: The code plots the figure as in the report and saves the global mean temperature model timeseries to the directory SavedModData
#
# ----------------------------------------------------------------------------------------------------
# Information on the software used
# - Software Version: iris2.2, python3.6.7
# - Landing page to access the software: https://scitools-iris.readthedocs.io/en/stable/
# - Operating System: Linux Scientific Linux 7.9 (Nitrogen)
# - Environment required to compile and run: See environment.yml
#  ----------------------------------------------------------------------------------------------------
#
#  License: Creative Commons Attribution 4.0 International License (http://creativecommons.org/licenses/by/4.0/)
#
# ----------------------------------------------------------------------------------------------------
# How to cite:
# When citing this code, please include both the code citation and the following citation for the related report component:
########################################################################

