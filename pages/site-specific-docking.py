import streamlit as st
import streamlit.components.v1 as components
from io import StringIO
import sys, os
import numpy as np
import pandas as pd
import numbers
import re
import glob

from Bio.PDB import PDBList
import pdb2pqr
import MDAnalysis as mda 
from MDAnalysis.coordinates import PDB
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Draw
from vina import Vina
import py3Dmol

def save_keys(key):
    st.session_state[key] = st.session_state["_" + key]

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def create_viewer(pdb_id, center, size):
    viewer = py3Dmol.view()
    viewer.setViewStyle({'style':'outline','color':'black','width':0.1})

    # add receptor (protein) model to py3Dmol viewer
    viewer.addModel(open(f'data/PDB_files/{pdb_id}','r').read(),format='pdb')
    Prot=viewer.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})

    viewer.addBox({"center": dict(x = center[0], y = center[1], z= center[2]), "dimensions": dict(d = abs(size[0]), h = abs(size[1]), w = abs(size[2])), "color" : "red", "opacity" : 0.5})
    components.html(viewer._make_html(), height = 500,width=500)

ligs = []
filenames = []
load_keys("protein_upload")
load_keys("ligand_upload")
protein = st.session_state.protein_upload
receptor_name = protein.name
ligands = st.session_state.ligand_upload
for uploaded_file in ligands:
    ligand_name = uploaded_file.name
    ligand_id = ligand_name.split(".")[:-1]
    ligs.append(ligand_id)
    filenames.append(f"data/MOL2_files/{ligand_id}.mol2")

col1, col2 = st.columns(2)
with col1:
    centerx = st.number_input("Select Center of Docking Box (X)", key = "x_center")
    centery = st.number_input("Select Center of Docking Box (Y)", key = "y_center")
    centerz = st.number_input("Select Center of Docking Box (Z)", key = "z_center")
with col2:
    sizex = st.slider("Select Size of Docking Box (X)", 1, 100, 5, key = "x_size")
    sizey = st.slider("Select Size of Docking Box (Y)", 1, 100, 5, key = "y_size")
    sizez = st.slider("Select Size of Docking Box (Z)", 1, 100, 5, key = "z_size")

center_dims = [st.session_state.x_center, st.session_state.y_center, st.session_state.z_center]
size_dims = [st.session_state.x_size, st.session_state.y_size, st.session_state.z_size]
create_viewer(receptor_name, center_dims, size_dims)
if st.button("Dock"):
    pass