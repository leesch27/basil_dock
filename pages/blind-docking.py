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
import prolif as plf
import py3Dmol
from IPython.display import display

#filter warnings
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

sys.path.insert(1, '../utilities/')
from utils import pdbqt_to_sdf
from basil_utils import get_prot_pockets_data, get_ifps, get_scores, save_dataframe, get_largest_array_column, expand_df, fill_df

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def create_viewer(pdb_id, revised_files):
    # view selected pockets in 3Dmol.js viewer prior to docking
    viewer = py3Dmol.view()
    viewer.removeAllModels()
    viewer.setViewStyle({'style':'outline','color':'black','width':0.1})

    # add receptor (protein) model to py3Dmol viewer
    viewer.addModel(open(f'data/PDB_files/{pdb_id}_protein.pdb','r').read(),format='pdb')
    Prot=viewer.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})

    #visualization of pockets
    colors = ['red', 'orange', 'yellow', 'green', 'blue', 'purple', 'magenta']
    a = 0
    for file in revised_files:
        viewer.addModel(open(file,'r').read(),format = 'pqr')
        pockets = viewer.getModel()
        pockets.setStyle({},{'sphere':{'color':colors[a],'opacity':0.5}}) 
        a += 1
        if a > 6:
            a = 0
    viewer.zoomTo()
    components.html(viewer._make_html(), height = 500,width=500)

def view_ligands(ligand):
    # view ligand in 3Dmol.js viewer prior to docking
    lig = f"{current_dir}/data/MOL2_files/{ligand}_H.mol2"
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    view.addModel(open(lig,'r').read(),format='mol2')
    ref_m = view.getModel()
    ref_m.setStyle({},{'stick':{'colorscheme':'greenCarbon','radius':0.2}})
    view.zoomTo()
    components.html(view._make_html(), height = 500,width=500)

def dock_vina(pdb_id, ligand, cavity, pocket_center, pocket_size, exhaust, pose):
    # iterate through each pocket and dock for a given ligand
    v = Vina(sf_name='vina')
    v.set_receptor(f'data/PDBQT_files/{pdb_id}_protein.pdbqt')
    v.set_ligand_from_file(f"data/PDBQT_files/{ligand}_H.pdbqt")
    if len(pocket_center) == len(pocket_size):
        for cav_num, id in enumerate(cavity):
            v.compute_vina_maps(center = pocket_center[cav_num], box_size = pocket_size[cav_num])
            v.dock(exhaustiveness=exhaust, n_poses=pose)
            v.write_poses(f"data/vina_out/{ligand}_vina_pocket_{id}.pdbqt", n_poses=pose, overwrite=True)
            # write output to sdf
            pdbqt_to_sdf(pdbqt_file=f"data/vina_out/{ligand}_vina_pocket_{id}.pdbqt",output=f"data/vina_out_2/{ligand}_pocket_{id}_vina_out_2.sdf")

def dock_smina(pdb_id, ligand, cavity, pocket_center, pocket_size, exhaust, pose):
    # iterate through each pocket and dock for a given ligand
    if len(pocket_center) == len(pocket_size):
        for cav_num, id in enumerate(cavity):
            rec = f'data/PDBQT_files/{pdb_id}_protein.pdbqt'
            lig = f'data/MOL2_files/{ligand}_H.mol2'
            outfile = f'data/smina_out/{ligand}_pocket_{id}_smina_out.sdf'
            pocket_center_cavity = pocket_center[cav_num]
            pocket_size_cavity = pocket_size[cav_num]
            # run smina docking
            try:
                smina = subprocess.run(["smina", "-r", rec, "-l", lig, "-o", outfile, "--center_x", str(pocket_center_cavity[0]), "--center_y", str(pocket_center_cavity[1]), "--center_z", str(pocket_center_cavity[2]), "--size_x", str(pocket_size_cavity[0]), "--size_y", str(pocket_size_cavity[1]), "--size_z", str(pocket_size_cavity[2]), "--exhaustiveness", str(exhaust), "--num_modes", str(pose)], text=True)
                mols = []
                # Rewrite sdf output files to add 3D tag
                with Chem.SDMolSupplier(f'data/smina_out/{ligand}_pocket_{id}_smina_out.sdf') as suppl:
                    for mol in suppl:
                        if mol is not None:
                            Chem.MolToMolBlock(mol)
                            mols.append(mol)
                with Chem.SDWriter(f"data/smina_out_2/{ligand}_pocket_{id}_smina_out_2.sdf") as w:
                    for mol in mols:
                        w.write(mol)
            except subprocess.CalledProcessError as smina:
                print(smina.stderr, end="")

# load session state variables from set parameters page
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
with st.status("Running fpocket on the protein, searching for binding pockets...") as status:
    # run fpocket to find potential binding pockets on protein
    try:
        with open(f"data/pocket_descriptors_{pdb_id}.csv", "w+") as out_file:
            fpocket = subprocess.run(["fpocket", "-f", f"data/PDB_files/{pdb_id}_protein.pdb", "-d"], text= True, check=True, stdout=out_file)
            status.update(label="fpocket run completed!")
    except subprocess.CalledProcessError as fpocket:
        print(fpocket.stderr, end="")
        status.update(label="fpocket run failed! Please check log file for details.")


prot_pockets = pd.read_csv(f'data/pocket_descriptors_{pdb_id}.csv',sep=' ',index_col=[0])
get_prot_pockets_data(current_dir, pdb_id, prot_pockets)
# initialize session state variables for pocket list and initial pocket list (copy of pocket list)
if 'pocket_list' not in st.session_state:
    st.session_state.pocket_list = prot_pockets
    st.session_state.pocket_list_initial = prot_pockets.copy()
st.write("List of potential binding pockets found")
st.dataframe(prot_pockets)

value = [round(x * 0.05, 2) for x in range(20)]
value = value[1:]
st.write("Select druggability cutoff to filter pockets for docking")
druggability_cutoff = st.selectbox("Select cutoff value", value, key = "_cutoff_val")
if st.button("Update Parameters"):
    prot_pockets_csv = st.session_state.pocket_list_initial
    prot_pockets2 = prot_pockets_csv[prot_pockets_csv['drug_score'] >= float(druggability_cutoff)]
    prot_pockets2.to_csv(f"data/protein_pockets_id_{pdb_id}_cutoff_{druggability_cutoff}.csv")
    
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
    # update session state variable for pocket list
    st.session_state.pocket_list = prot_pockets2
    st.dataframe(prot_pockets2)
    create_viewer(pdb_id, revised_files)

st.write("View ligands prior to docking")
view_lig= st.selectbox("Select ligand to view", ligs, key = "_view_lig")
view_ligands(st.session_state._view_lig)

if st.button("Dock!"):
    pocket_center = []
    pocket_size = []
    cav_ids = []
    pockets = st.session_state.pocket_list
    for pocket in pockets.index:
        c_x = pockets.loc[pocket,'center_x']
        c_y = pockets.loc[pocket,'center_y']
        c_z = pockets.loc[pocket,'center_z']
        s_x = pockets.loc[pocket,'size_x']
        s_y = pockets.loc[pocket,'size_y']
        s_z = pockets.loc[pocket,'size_z']
        pocket_center.append([c_x, c_y, c_z])
        pocket_size.append([s_x, s_y, s_z])
        cav_ids.append(pocket)
    
    if docking_engine == "Vina":
        with st.status("Docking with Autodock Vina...") as status_vina:
            dataPath = os.path.join(current_dir, "data")
            vina_out = os.path.join(current_dir, "data", "vina_out")
            vina_out_2 = os.path.join(current_dir, "data", "vina_out_2")

            bool_data = os.path.exists(dataPath)
            bool_vina = os.path.exists(vina_out)
            bool_vina2 = os.path.exists(vina_out_2)

            if(bool_data == False):
                print("ERROR: Cannot find 'data' folder. Make sure you are in the correct directory.")
                status_vina.update(label="Cannot find 'data' folder. Make sure you are in the correct directory.")
            elif(bool_vina == False):
                print("ERROR: Cannot find 'vina_out' folder. Creating folder...")
                os.mkdir(vina_out)
            elif(bool_vina2 == False):
                print("ERROR: Cannot find 'vina_out_2' folder. Creating folder...")
                os.mkdir(vina_out_2)
            for ligand in ligs:
                st.write(f"Docking ligand {ligand}...")
                dock_vina(pdb_id, ligand, cav_ids, pocket_center, pocket_size, exhaustiveness, num_poses)
            status_vina.update(label="Docking with Autodock Vina completed!")
    else:
        with st.status("Docking with Smina...") as status_smina:
            dataPath = os.path.join(current_dir, "data")
            smina_out = os.path.join(current_dir, "data", "smina_out")
            smina_out_2 = os.path.join(current_dir, "data", "smina_out_2")

            bool_data = os.path.exists(dataPath)
            bool_smina = os.path.exists(smina_out)
            bool_smina2 = os.path.exists(smina_out_2)

            if(bool_data == False):
                print("ERROR: Cannot find 'data' folder. Make sure you are in the correct directory.")
                status_smina.update(label="Cannot find 'data' folder. Make sure you are in the correct directory.")
            elif(bool_smina == False):
                print("ERROR: Cannot find 'smina_out' folder. Creating folder...")
                os.mkdir(smina_out)
            elif(bool_smina2 == False):
                print("ERROR: Cannot find 'smina_out_2' folder. Creating folder...")
                os.mkdir(smina_out_2)
            for ligand in ligs:
                st.write(f"Docking ligand {ligand}...")
                dock_smina(pdb_id, ligand, cav_ids,pocket_center, pocket_size, exhaustiveness, num_poses)
            status_smina.update(label="Docking with Smina completed!")
    
    # save results to dataframe
    with st.status("Interpreting docking results...") as status_results:
        engine_name = docking_engine.lower()
        prot_mol = Chem.MolFromPDBFile(f"data/PDB_files/{pdb_id}_protein_H.pdb")
        protein_plf = plf.Molecule.from_rdkit(prot_mol)
        st.write(f"Obtaining interaction fingerprints...")
        all_df, all_ifps, all_ligand_plf, ligand_plf_descriptors = get_ifps("Blind docking", engine_name, ligs, protein_plf, pockets)
        st.write(f"Sorting through scores...")
        scores = get_scores("Blind docking", engine_name, ligs, pockets)
        ligs_in_order = []
        pocks_in_order = []
        for value in ligand_plf_descriptors:
            ligand_name = value.split(",")[0]
            ligs_in_order.append(ligand_name)
            pocket_name = value.split(",")[1]
            pock_name_isolated = pocket_name.split(" ")[-1]
            pocks_in_order.append(int(pock_name_isolated))
        st.write(f"Creating dataframe of results...")
        df = pd.concat([d for d in all_df], axis=0, ignore_index=False, sort=False).reset_index()
        df.insert(1, "Score", pd.Series(scores))
        df.insert(1, "Ligand", pd.Series(ligs_in_order))
        df.insert(1, "Pocket", pd.Series(pocks_in_order))
        df = df.fillna(0)
        # save csv
        save_dataframe(df, engine_name, pdb_id, ligs)
        st.write(f"Expanding dataframe...")
        df2 = df[["Frame", "Score", "Ligand", "Pocket", "UNL1"]].copy()
        largest_array_column = get_largest_array_column(df, "Blind docking")
        expand_df(all_ifps, df, df2, largest_array_column)
        st.write(f"Filling expanded dataframe...")
        fill_df(df2, all_ifps, all_ligand_plf, largest_array_column)
        save_dataframe(df2, engine_name, pdb_id, ligs, csv_name = f"{pdb_id}_{str(len(ligs))}_ligands_docking_information_{engine_name}_extended")
        status_results.update(label="Results saved! To analyze docking results, click the 'Analyze docking results' link below.")

st.page_link("pages/docking-analysis.py", label="Analyze docking results", icon="📊")
st.page_link("pages/set-parameters.py", label="Return to parameter selection", icon="🏠")