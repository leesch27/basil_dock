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

sys.path.insert(1, 'utilities/')
from ligandsplitter.ligandsplitter.basefunctions import create_folders, convert_type
from ligandsplitter.ligandsplitter.ligandsplit import File_Info, Ligand, retrieve_pdb_file, get_mol2_info, get_ligands, find_ligands_unique, write_mol2, separate_mol2_ligs, isolate_by_method

# change formatting of sidebar, remove site specific page. only want it accessible if option is selected in form

def create_viewer(pdb_id, center, size):
    viewer = py3Dmol.view()
    viewer.setViewStyle({'style':'outline','color':'black','width':0.1})

    # add receptor (protein) model to py3Dmol viewer
    viewer.addModel(open(f'data/PDB_files/{pdb_id}','r').read(),format='pdb')
    Prot=viewer.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})

    viewer.addBox({"center": dict(x = center[0], y = center[1], z= center[2]), "dimensions": dict(d = abs(size[0]), h = abs(size[1]), w = abs(size[2])), "color" : "red", "opacity" : 0.5})
    components.html(viewer._make_html(), height = 500,width=500)

def save_keys(key):
    st.session_state[key] = st.session_state["_" + key]

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

with st.form("enter_docking_parameters"):
    header = st.columns([1,1])
    header[0].subheader("basil dock Settings")
    header[1].subheader("Receptor and Ligand Selection Settings")

    row1 = st.columns([1,1])
    docking_engine = row1[0].selectbox("Select Docking Engine to Use", ["Vina", "Smina"], key = "docking_engine")
    protein_retrieval = row1[1].file_uploader("Upload Protein Receptor To Use", accept_multiple_files= False, type=["cif", "pdb", "ent"], key = "_protein_upload")
    
    row2 = st.columns([1,1])
    docking_method = row2[0].selectbox("Select Docking Method to Use", ["Blind Docking", "Site-Specific Docking"], key = "docking_method")
    ligand_retrieval = row2[1].file_uploader("Upload Ligand/s To Use", accept_multiple_files= True, type=["mol2", "sdf"], key = "_ligand_upload")

    row3 = st.columns([1,1])
    exhaustiveness = row3[0].slider("Select Exhaustiveness Value", 1, 16, 5, key = "exhaust_val")

    row4 = st.columns([1,1])
    poses = row4[0].slider("Select Number of Poses to Generate", 1, 10, 5, key = "poses_val")

    row5 = st.columns([1,1])
    #wd = row5[0].file_uploader("Select Folder to Download Results To", accept_multiple_files="directory", key = "poses_val")
    submitted = row5[0].form_submit_button("Submit Parameters")
    
    if submitted:
        # create new pages (blind dock and site-specific dock) so forms can be used
        current_dir = create_folders()
        ligs = []
        filenames = []
        protein = st.session_state._protein_upload
        ligands = st.session_state._ligand_upload

        # get receptor data
        if protein is not None:
            receptor_name = protein.name
            stringio = StringIO(protein.getvalue().decode("utf-8"))
            bytes_data = stringio.read()
            with open(f"{current_dir}/data/PDB_files/{receptor_name}", "w+") as datafile:
                datafile.write(bytes_data)
        
        # get ligand data
        for uploaded_file in ligands:
            ligand_name = uploaded_file.name
            ligand_id = ligand_name.split(".")[:-1]
            ligand_extension = ligand_name.split(".")[-1]
            if ligand_extension == "sdf":
                pdb_mol2 = [m for m in pybel.readfile(filename = ligand_name, format='sdf')][0]
                out_mol2 = pybel.Outputfile(filename = f"data/MOL2_files/{ligand_id}.mol2", overwrite = True, format='mol2')
                out_mol2.write(pdb_mol2)
                ligs.append(ligand_id)
                filenames.append(f"data/MOL2_files/{ligand_id}.mol2")
            else:
                stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                bytes_data = stringio.read()
                with open(f"{current_dir}/data/PDB_files/{ligand_name}", "w+") as datafile:
                    datafile.write(bytes_data)
                ligs.append(ligand_id)
                filenames.append(f"data/MOL2_files/{ligand_id}.mol2")
        
        dock_method = st.session_state.docking_method
        if dock_method == "Blind Docking":
            # remake fpocket function, dock
            # separate blind docking page to choose drugability score cutoff?
            pass
        else:
            save_keys("protein_upload")
            save_keys("ligand_upload")
            st.switch_page("pages/site-specific-docking.py")

#st.page_link("advanced-search-ligands.py", label="")
#st.page_link("advanced-search-proteins.py", label="")
#st.page_link("create-ligands-from-smiles.py", label="")
#st.page_link("derive-ligands.py", label="")
