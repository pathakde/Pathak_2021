# Quenching, Mergers and Age Profiles for z=2 Galaxies in IllustrisTNG

In this repository, you can find resources to replicate the analysis of IllustrisTNG data presented in Pathak et al. (2021) [insert link].

We also provide some resources for streamlining your analysis of relevant IllustrisTNG data with detailed instructions to download and process this data.

<br />

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
This notebook will help you generate all five figures from the paper [insert link].

<br />

If you have any questions or comments, please reach out to pathakde@grinnell.edu!
