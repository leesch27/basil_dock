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

**1. Obtaining notebooks**

**2. Creating a conda environment**

**3. Opening and using notebooks**

## Limitations

## Citation

## License
