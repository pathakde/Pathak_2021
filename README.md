# Quenching, Mergers and Age Profiles for z=2 Galaxies in IllustrisTNG

In this repository, you can find resources to replicate the analysis of IllustrisTNG data presented in Pathak et al. (2021) https://arxiv.org/abs/2105.02234.

We also provide some resources for streamlining your analysis of relevant IllustrisTNG data with detailed instructions to download and process this data. 

We include **two additional figures** to complement the results presented in Pathak et al. (2020), and to facilitate further discussion. Scroll to the bottom of `Figures.ipynb` and look under **Extras** to view the new figures.

## Note: Selected Suess et al. (2020) data included:
Some processed data from Suess et al. (2020) has been included in `Suess2020_data.hdf5` for reproducing the comparisons with observations in Figure 2. A huge thank you to **Wren Suess** for sharing her data! You can find her paper at https://iopscience.iop.org/article/10.3847/2041-8213/abacc9 (arXiv: https://arxiv.org/abs/2008.02817).

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

If you have any questions or comments, please reach out to pathakde@grinnell.edu!
