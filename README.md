# **basil_dock - Accessible and Flexible Molecular Docking**
[**Description**](#description) | [**Requirements**](#requirements)| [**Installation**](#installation) | [**Limitations**](#limitations) | [**Citation**](#citation) | [**License**](#license)

## Description
basil_dock is a series of Jupyter notebooks designed to perform interactive and flexible molecular docking procedures, regardless of the user's skill level. <br>

**1. Docking preparation**
> Description
  - fetch PDB file from rcsb.org using PDB ID and download to computer
  - split PDB file into separate files for protein (PDB format) and ligand/s (PDB and MOL2 format)
  - add additional ligands as desired from local MOL2 files on a personal computer or by inputting ligand/s as a SMILES string
  - sanitize protein and ligand files and convert to required formats (protein - PDBQT format, ligand - MOL2/PDBQT formats)
  - determine likely pockets on the protein receptor
  - visualize protein and ligand/s <br>
  
**2. Molecular docking**
> Description
  - choose desired ligand/s to dock to protein receptor
  - choose between site-specific and blind docking
  - choose between docking engines to use in docking (Vina and Smina)
  - view interaction fingerprints between ligand/s and protein
  - visualize protein and ligand docking pose  <br>
  
**3. Data manipulation and collection**
> Description
  - create hypothetical ligands by modifying different functional groups on original ligand/s
  - choose desired modified ligand/s (derivatives) to dock to protein receptor
  - choose between site-specific and blind docking
  - choose between docking engines to use in docking (Vina and Smina)
  - view interaction fingerprints between derivative/s and protein
  - visualize protein and derivative docking pose <br>

**4. Machine learning analysis**
> Description
  - determine likelihood of ligands/derivatives being orally bioactive using Lipinski's Rule of Five <br>

## Requirements

## Installation

**1. Obtaining notebooks** <br>
To download this repository (containing all notebooks and an environment.yml file with all necessary libraries) to your computer, press the green code button, select "Local", and download the zip file. You can then unzip the file and move it into whatever computer folder you would like.
You can also create a virtual codespace, which will not download a folder to your computer but will instead create a virtual environment for you to run the notebooks.

**2. Creating a conda environment** <br>
To make library installation easier, an environment.yml file has been included with this repository. (explain anaconda and environments)
After installing Anaconda on your computer, (explain).
```
conda env create -f environment.yml
```
This will create a new conda environment called "basil_dock", which will be used to access all of the libraries that are needed to run this notebook series.

**3. Opening and using notebooks** <br>
This notebook series can be accessed through a few different applications. JupyterLab, Jupyter Notebooks, and VSCode (explain)
  - JupyterLab
    1. Open a new terminal window and change the conda environment to "basil_dock" <br>
    2. Open JupyterLab <br>
```
conda activate basil_dock
```
```
jupyter lab
```
  - Jupyter Notebooks
    1. Open a new terminal window and change the conda environment to "basil_dock" <br>
    2. Open Jupyter Notebook <br>
```
conda activate basil_dock
``` 
```
jupyter nbclassic
```
  - VSCode
    1. Download VSCode on their website
    2. Right-click on the notebook you would like to open. VSCode should be an option. Select VSCode.
    3. Change kernel to "basil_dock" in the app

## Limitations

## Citation

## License
