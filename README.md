# Quenching, Mergers and Age Profiles for z=2 Galaxies in IllustrisTNG

In this repository, you can find resources to replicate the analysis of IllustrisTNG data presented in [Pathak et al. (2021)](https://arxiv.org/abs/2105.02234).

We also provide some resources for streamlining your analysis of relevant IllustrisTNG data with detailed instructions to download and process this data. 

We include **two additional figures** to complement the results presented in Pathak et al. (2020), and to facilitate further discussion. Scroll to the bottom of `Figures.ipynb` and look under **Extras** to view the new figures.


<br />

# Getting Started:

## 1. Download all code and data included in the repository.

All data acquisition and analysis in this project was done in Python 3.8.3. 

We suggest downloading all notebooks, modules, and data files included in this repository before starting your analysis.  

## 2. Run the notebook titled `Data_Acquisition.ipynb`.

Follow the steps outlined in the notebook. If this is your first time working with IllustrisTNG data, *do not skip Step 0.* 
Go to  www.tng-project.org to set up an account and get an API key. Detailed steps are included in the notebook. You will need this API key to download TNG data from the server.

Note that downloading data on individual galaxies from the TNG server requires a stable internet connection. Most of the subsequent analysis can be done offline.

This notebook also provides instructions for storing data on individual galaxies, calculating halo properties from downloaded data, and compiling data on selected populations.

## 3. Run the notebook titled `Figures.ipynb`.

This notebook analyzes the data that was downloaded and processed in `Data_Acquisition.ipynb`. 
This notebook will help you generate all five figures from the paper https://arxiv.org/abs/2105.02234.

<br />

# Data Sources:
This repository does not host any proprietary data. All data included in this repository has been processed before release. Appropriate citations for the original data (not included in repository) have been included below.

### IllustrisTNG Collaboration
Some post-processed data from the [IllustrisTNG public data release](https://www.tng-project.org/data/docs/background/#sec4) have been included.

Details can be found at:
  [(Nelson et al. 2019a)](http://arxiv.org/abs/1812.05609);
  [Pillepich et al. (2018b)](https://arxiv.org/abs/1707.03406);
  [Springel et al. (2018)](https://arxiv.org/abs/1707.03397);
  [Nelson et al. (2018a)](https://arxiv.org/abs/1707.03395);
  [Naiman et al. (2018)](https://arxiv.org/abs/1707.03401);
  [Marinacci et al. (2018)](https://arxiv.org/abs/1707.03396).

### Suess et al. (2020)
A huge thank you to **Wren Suess** for sharing her data!
Some data from [Suess et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020ApJ...899L..26S/abstract) was processed similar to our halo selection function, and population percentiles have been included in `Suess2020_data.hdf5` for reproducing the comparisons with observations in Figure 2. 

### TNG Stellar Assembly Files
Calculations for halo data in the the field `maximum_merger_ratio_30kpc_current_fraction` in `galaxy_population_data_2.hdf5` was done by post-processing data from Stellar Assembly Files provided by Vicente Rodriguez-Gomez. Details can be found at [Rodriguez-Gomez et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.2371R/abstract).

<br />

If you have any questions or comments, please reach out to pathakde@grinnell.edu!
