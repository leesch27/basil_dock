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
from rdkit.Chem import Draw
from vina import Vina
import py3Dmol

sys.path.insert(1, '../utilities/')
from basil_utils import get_ifps, get_scores, get_largest_array_column, expand_df, fill_df

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def create_viewer(pdb_id, center, size):
    viewer = py3Dmol.view()
    viewer.removeAllModels()
    viewer.setViewStyle({'style':'outline','color':'black','width':0.1})

    # add receptor (protein) model to py3Dmol viewer
    viewer.addModel(open(f'{pdb_id}','r').read(),format='pdb')
    Prot=viewer.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})

    viewer.addBox({"center": dict(x = center[0], y = center[1], z= center[2]), "dimensions": dict(d = abs(size[0]), h = abs(size[1]), w = abs(size[2])), "color" : "red", "opacity" : 0.5})
    viewer.zoomTo()
    components.html(viewer._make_html(), height = 500,width=500)

def view_ligands(ligand):
    lig = f"{current_dir}/data/MOL2_files/{ligand}_H.mol2"
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    view.addModel(open(lig,'r').read(),format='mol2')
    ref_m = view.getModel()
    ref_m.setStyle({},{'stick':{'colorscheme':'greenCarbon','radius':0.2}})
    view.zoomTo()
    components.html(view._make_html(), height = 500,width=500)

def dock_vina(pdb_id, ligand, centers, sizes, exhaust, pose):
    v = Vina(sf_name='vina')
    v.set_receptor(f'{current_dir}/data/PDBQT_files/{pdb_id}_protein.pdbqt')
    v.set_ligand_from_file(f"{current_dir}/data/PDBQT_files/{ligand}_H.pdbqt")
    v.compute_vina_maps(center = centers, box_size = sizes)
    v.dock(exhaustiveness=exhaust, n_poses=pose)
    v.write_poses(f"{current_dir}/data/vina_out/{ligand}.pdbqt", n_poses=pose, overwrite=True)

def dock_smina(pdb_id, ligand, centers, sizes, exhaust, pose):
    rec = f'{current_dir}/data/PDBQT_files/{pdb_id}_protein.pdbqt'
    lig = f'{current_dir}/data/MOL2_files/{ligand}_H.mol2'
    outfile = f'{current_dir}/data/smina_out/{ligand}_smina_out.sdf'
    try:
        smina = subprocess.run(["smina", "-r", rec, "-l", lig, "-o", outfile, "--center_x", centers[0], "--center_y", centers[1], "--center_z", centers[2], "--size_x", sizes[0], "--size_y", sizes[1], "--size_z", sizes[2], "--exhaustiveness", exhaust, "--num_modes", pose])
        mols = []
        with Chem.SDMolSupplier(f'{current_dir}/data/smina_out/{ligand}_smina_out.sdf') as suppl:
            for mol in suppl:
                if mol is not None:
                    Chem.MolToMolBlock(mol)
                    mols.append(mol)
        with Chem.SDWriter(f"{current_dir}/data/smina_out_2/{ligand}_smina_out_2.sdf") as w:
            for mol in mols:
                w.write(mol)    
    except subprocess.CalledProcessError as smina:
        print(smina.stderr, end="")

ligs = []
filenames = []
load_keys("current_dir")
load_keys("docking_engine")
load_keys("exhaust_val")
load_keys("poses_val")
load_keys("ligs")
load_keys("filenames")
load_keys("filenames_H")
load_keys("filenames_pdbqt")
load_keys("pdb_id")

current_dir = st.session_state.current_dir
pdb_id = st.session_state.pdb_id
receptor_name = f"{current_dir}/data/PDB_files/{pdb_id}_protein.pdb"
ligs = st.session_state.ligs
filenames = st.session_state.filenames
filenames_H = st.session_state.filenames_H
filenames_pdbqt = st.session_state.filenames_pdbqt
docking_engine = st.session_state.docking_engine
num_poses = st.session_state.poses_val
exhaustiveness = st.session_state.exhaust_val

st.title("Site-Specific Docking Parameters")
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

st.write("View ligands prior to docking")
docking_engine = st.selectbox("Select ligand to view", ligs, key = "_view_lig")
view_ligands(st.session_state._view_lig)

if st.button("Dock!"):
    if docking_engine == "Vina":
        dataPath = os.path.join(current_dir, "data")
        vina_out = os.path.join(current_dir, "data", "vina_out")
        vina_out_2 = os.path.join(current_dir, "data", "vina_out_2")

        bool_data = os.path.exists(dataPath)
        bool_vina = os.path.exists(vina_out)
        bool_vina2 = os.path.exists(vina_out_2)

        if(bool_data == False):
            print("ERROR: Cannot find 'data' folder. Make sure you are in the correct directory.")
        elif(bool_vina == False):
            print("ERROR: Cannot find 'vina_out' folder. Creating folder...")
            os.mkdir(vina_out)
        elif(bool_vina2 == False):
            print("ERROR: Cannot find 'vina_out_2' folder. Creating folder...")
            os.mkdir(vina_out_2)
        for ligand in ligs:
            dock_vina(pdb_id, ligand, center_dims, size_dims, exhaustiveness, num_poses)
    else:
        dataPath = os.path.join(current_dir, "data")
        smina_out = os.path.join(current_dir, "data", "smina_out")
        smina_out_2 = os.path.join(current_dir, "data", "smina_out_2")

        bool_data = os.path.exists(dataPath)
        bool_smina = os.path.exists(smina_out)
        bool_smina2 = os.path.exists(smina_out_2)

        if(bool_data == False):
            print("ERROR: Cannot find 'data' folder. Make sure you are in the correct directory.")
        elif(bool_smina == False):
            print("ERROR: Cannot find 'smina_out' folder. Creating folder...")
            os.mkdir(smina_out)
        elif(bool_smina2 == False):
            print("ERROR: Cannot find 'smina_out_2' folder. Creating folder...")
            os.mkdir(smina_out_2)
        for ligand in ligs:
            dock_smina(pdb_id, ligand, center_dims, size_dims, exhaustiveness, num_poses)

st.page_link("pages/docking-analysis.py", label="Analyze docking results", icon="📊")
st.page_link("pages/set-parameters.py", label="Return to parameter selection", icon="🏠")