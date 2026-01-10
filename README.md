# **basil_dock - Accessible and Flexible Molecular Docking**
[**Description**](#description) | [**Requirements**](#requirements)| [**Installation**](#installation) | [**Usage**](#usage) | [**Limitations**](#limitations) | [**Citation**](#citation) | [**License**](#license)

## Description
basil_dock is a series of Jupyter notebooks designed to perform interactive and flexible molecular docking procedures, regardless of the user's skill level. <br>

**1. Docking preparation**
> Retrieve, sanitize, and convert beween file formats for protein receptor and ligands to be used.
  - fetch PDB file from rcsb.org using PDB ID and download to computer
  - split PDB file into separate files for protein (PDB format) and ligand/s (PDB and MOL2 format)
  - add additional ligands as desired from local MOL2 files on a personal computer or by inputting ligand/s as a SMILES string
  - sanitize protein and ligand files and convert to required formats (protein - PDBQT format, ligand - MOL2/PDBQT formats)
  - determine likely pockets on the protein receptor
  - visualize protein and ligand/s <br>
  
**2. Molecular docking**
> Model potential binding poses between a prepared protein receptor and a ligand
  - choose desired ligand/s to dock to protein receptor
  - choose between site-specific and blind docking
  - choose between docking engines to use in docking (Vina and Smina)
  - view interaction fingerprints between ligand/s and protein
  - visualize protein and ligand docking pose  <br>
  
**3. Data manipulation and collection**
> Determine how modifications to ligand composition affects binding affinity and complex formation
  - create hypothetical ligands by modifying different functional groups on original ligand/s
  - choose desired modified ligand/s (derivatives) to dock to protein receptor
  - choose between site-specific and blind docking
  - choose between docking engines to use in docking (Vina and Smina)
  - view interaction fingerprints between derivative/s and protein
  - visualize protein and derivative docking pose <br>

## Requirements
All package and software requirements for basil_dock can be found in the YML files included with this repository. 

## Installation

**1. Obtaining notebooks** <br>
To download this repository (containing all notebooks and an environment.yml file with all necessary libraries) to your computer, press the green code button, select "Local", and download the zip file.
You can also create a virtual codespace, which will not download a folder to your computer but will instead create a virtual environment for you to run the notebooks.

From this point, there are two ways to set up the notebooks; you can either use the notebook entitled "00-setup-docking.ipynb" (See "Guided Setup") or set up all environments unassisted (See "Unassisted Setup").

**2A. Guided Setup** <br>
For assisted setup, open the first notebook "00-setup-docking.ipynb" in VSCode and execute cells as described in the notebook. 

**2B. Unassisted Setup** <br>
To make library installation easier, an environment.yml file has been included with this repository. After installing Anaconda on your computer, the conda command can be used to install all necessary dependencies.

Navigate to evirnoment files in the basil_dock folder
```
cd basil_dock/conda-environments
```

1. MacOS Installation
```
conda env create -f environment_macos.yml
```
```
conda activate basil_dock
```

2. Windows Installation
```
conda env create -f environment_windows.yml
```
```
conda activate basil_dock
```

3. Linux Installation
```
conda env create -f environment_linux.yml
```
```
conda activate basil_dock
```

This will create a new conda environment called "basil_dock", which will be used to access all of the libraries that are needed to run this notebook series.

If you would rather install using pip, use the command below to install all dependencies. Please note that not all libraries used in basil_dock are able to be accessed via pip, so some functionalities may be missing.
```
cd basil_dock/pip-installation
```
```
pip install requirements.txt
```

## Usage
basil_dock can be executed as a series of jupyter notebooks using JupyterLab or VSCode, or as a Streamlit application

**1. Jupyter Notebooks** <br>
**. Opening and using notebooks** <br>
This notebook series can be accessed through a few different applications. JupyterLab, Jupyter Notebooks, and VSCode 
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

**2. Streamlit Application** <br>
basil_dock can also be accessed as a Streamlit application either locally or via the Streamlit Community Cloud. 

- Accessing Streamlit application locally
```
conda activate basil_dock
```
```
streamlit run basildock.py
```

- Accessing Streamlit application via Streamlit Community Cloud
  1. To access the current stable version of basil_dock, go to https://basildock.streamlit.app/
  2. To access the development version of basil_dock, go to https://basildock-dev.streamlit.app/

## Limitations
basil_dock is still in its early stages of development. Many ideas that were proposed have not yet been implemented.  
1. Docking preparation
  - Nucleic acid - protein binding simulations

2. Docking functionalities
  - Addition of more docking engines to increase customizability

3. Data manipulation
  - Integration of more ligand functional groups for derivative creation
  - Include site directed mutagenisis of protein receptor

Additionally, due to how new the program is, errors are inevitable. If you come across a bug, please submit a request under the "issues" tab.

## Citation
Schoneman, Lee, "basil_dock: Development of Accessible and Customizable Molecular Docking Procedures" (2025).

## License
basil_dock is licensed under GLP-3.0