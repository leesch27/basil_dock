import streamlit as st
import streamlit.components.v1 as components
from io import StringIO
import sys, os
import numpy as np
import pandas as pd
import numbers
import re
import glob
import subprocess

from Bio.PDB import PDBList
import MDAnalysis as mda 
from MDAnalysis.coordinates import PDB
from openbabel import pybel
from rdkit import Chem

sys.path.insert(1, 'utilities/')
from ligandsplitter.ligandsplitter.basefunctions import create_folders, convert_type

def save_keys(key):
    st.session_state[key] = st.session_state["_" + key]

current_dir = create_folders()
st.session_state._current_dir = current_dir
save_keys("current_dir")
# change formatting of sidebar, remove site specific page. only want it accessible if option is selected in form
pages = {
    "Set parameters to be used for docking": [
        st.Page("pages/set-parameters.py", title="Parameters")
    ],
    "Docking methods": [
        st.Page("pages/blind-docking.py", title="Dock ligands to predicted binding sites in a receptor"),
        st.Page("pages/site-specific-docking.py", title="Dock ligands to a user-defined area of a receptor"),
        st.Page("pages/docking-analysis.py", title="Analyze docking results")
    ],
    "Get ligands to use in docking": [
        st.Page("pages/advanced-search-ligands.py", title="Search for ligands using RCSB PDB Advanced Search"),
        st.Page("pages/create-ligands-from-smiles.py", title="Create ligands using SMILES strings"),
        st.Page("pages/derive-ligands.py", title="Derive ligands from existing ligands")
    ],
    "Get receptors to use in docking": [
        st.Page("pages/advanced-search-proteins.py", title="Search for receptors using RCSB PDB Advanced Search")
    ]
}

pg = st.navigation(pages)
pg.run()

