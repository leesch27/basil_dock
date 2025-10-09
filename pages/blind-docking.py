import streamlit as st
import streamlit.components.v1 as components
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
from IPython.display import display

sys.path.insert(1, '../utilities/')
from basil_utils import get_prot_pockets_data, get_ifps, get_scores, get_largest_array_column, expand_df, fill_df

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def create_viewer(pdb_id, revised_files):
    viewer = py3Dmol.view()
    viewer.removeAllModels()
    viewer.setViewStyle({'style':'outline','color':'black','width':0.1})

    # add receptor (protein) model to py3Dmol viewer
    viewer.addModel(open(f'data/PDB_files/{pdb_id}_protein.pdb','r').read(),format='pdb')
    Prot=viewer.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})

    #visualization ligands
    colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple', 'magenta']
    a = 0
    for file in revised_files:
        viewer.addModel(open(file,'r').read(),format = 'pqr')
        pockets = viewer.getModel()
        pockets.setStyle({},{'sphere':{'color':colors[a],'opacity':0.5}}) 
        a += 1
        if a > 6:
            a = 0
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

def dock_vina(pdb_id, ligand, pocket_center, pocket_size, exhaust, pose):
    v = Vina(sf_name='vina')
    v.set_receptor(f'data/PDBQT_files/{pdb_id}_protein.pdbqt')
    v.set_ligand_from_file(f"data/PDBQT_files/{ligand}_H.pdbqt")
    if len(pocket_center) == len(pocket_size):
        for pock_num, pocket in enumerate(pocket_center):
            v.compute_vina_maps(center = pocket, box_size = pocket_size[pock_num])
            v.dock(exhaustiveness=exhaust, n_poses=pose)
            v.write_poses(f"data/vina_out/{ligand}_vina_pocket_{pock_num}.pdbqt", n_poses=5, overwrite=True)

def dock_smina(pdb_id, ligand, pocket_center, pocket_size, exhaust, pose):
    if len(pocket_center) == len(pocket_size):
        for pock_num, pocket in enumerate(pocket_center):
            rec = f'data/PDBQT_files/{pdb_id}_protein.pdbqt'
            lig = f'data/MOL2_files/{ligand}_H.mol2'
            outfile = f'data/smina_out/{ligand}_pocket_{pock_num}_smina_out.sdf'
            try:
                smina = subprocess.run(["smina", "-r", rec, "-l", lig, "-o", outfile, "--center_x", pocket_center[0], "--center_y", pocket_center[1], "--center_z", pocket_center[2], "--size_x", pocket_size[0], "--size_y", pocket_size[1], "--size_z", pocket_size[2], "--exhaustiveness", exhaust, "--num_modes", pose])
                mols = []
                with Chem.SDMolSupplier(f'data/smina_out/{ligand}_pocket_{pock_num}_smina_out.sdf') as suppl:
                    for mol in suppl:
                        if mol is not None:
                            Chem.MolToMolBlock(mol)
                            mols.append(mol)
                with Chem.SDWriter(f"data/smina_out_2/{ligand}_pocket_{pock_num}_smina_out_2.sdf") as w:
                    for mol in mols:
                        w.write(mol)
            except subprocess.CalledProcessError as smina:
                print(smina.stderr, end="")

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
ligs = st.session_state.ligs
filenames = st.session_state.filenames
filenames_H = st.session_state.filenames_H
filenames_pdbqt = st.session_state.filenames_pdbqt
docking_engine = st.session_state.docking_engine
num_poses = st.session_state.poses_val
exhaustiveness = st.session_state.exhaust_val

st.title("Blind Docking Parameters")
try:
    out_csv = f"{current_dir}/data/PDB_files/pocket_descriptors_{pdb_id}.csv"
    fpocket = subprocess.run(["fpocket", "-f", f"{current_dir}/data/PDB_files/{pdb_id}_protein.pdb", f"-d > {out_csv}"])
except subprocess.CalledProcessError as fpocket:
    print(fpocket.stderr, end="")
try:
    prot_pockets = pd.read_csv(f'{current_dir}/data/PDB_files/pocket_descriptors_{pdb_id}.csv',sep=' ',index_col=[0])
    get_prot_pockets_data(current_dir, pdb_id, prot_pockets)
    st.write("Potential binding pockets found")
    st.dataframe(prot_pockets)
except:
    pass
value = [round(x * 0.05, 2) for x in range(20)]
value = value[1:]
st.write("To rre")
druggability_cutoff = st.selectbox("Select cutoff value", value, key = "_cutoff_val")
if st.button("Update Parameters"):
    prot_pockets2 = prot_pockets[prot_pockets['drug_score'] >= float(druggability_cutoff)]
    prot_pockets2.to_csv(f"{current_dir}/data/protein_pockets_id_{pdb_id}_cutoff_{druggability_cutoff}.csv")
    revised_files = []
    pocketPath = os.path.join(current_dir, "data", "PDB_files", str(pdb_id) + "_protein_out", "*.pqr")
    pocketFiles = glob.glob(pocketPath)
    for file in pocketFiles:
        split_1 = file.split("/")[-1]
        split_2 = split_1.split("_")[0]
        index_num = re.findall(r'\d+', split_2)
        index_num2 = ''.join(str(x) for x in index_num)
        if int(index_num2) in prot_pockets2.index:
            revised_files.append(file)
    st.dataframe(prot_pockets2)
    create_viewer(pdb_id, revised_files)

st.write("View ligands prior to docking")
docking_engine = st.selectbox("Select ligand to view", ligs, key = "_view_lig")
view_ligands(st.session_state._view_lig)

if st.button("Dock!"):
    pocket_center = []
    pocket_size = []
    for pocket in prot_pockets.index:
        c_x = prot_pockets.loc[pocket,'center_x']
        c_y = prot_pockets.loc[pocket,'center_y']
        c_z = prot_pockets.loc[pocket,'center_z']
        s_x = prot_pockets.loc[pocket,'size_x']
        s_y = prot_pockets.loc[pocket,'size_y']
        s_z = prot_pockets.loc[pocket,'size_z']
        pocket_center.append([c_x, c_y, c_z])
        pocket_size.append([s_x, s_y, s_z])
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

st.page_link("pages/set-parameters.py", label="Return to parameter selection", icon="🏠")