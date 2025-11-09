import streamlit as st
import streamlit.components.v1 as components
from io import StringIO, BytesIO
import sys, os
import numpy as np
import pandas as pd
import numbers
import re
import glob
import subprocess
import zipfile
import base64

from Bio.PDB import PDBList
import MDAnalysis as mda 
from MDAnalysis.coordinates import PDB
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Draw
from vina import Vina
import prolif as plf
import py3Dmol

#filter warnings
import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#sys.path.insert(1, '../utilities/')
from utilities.utils import pdbqt_to_sdf
from utilities.basil_utils import get_ifps, get_scores, save_dataframe, get_largest_array_column, expand_df, fill_df

cur_dir = os.getcwd()
local = True
if "mount/src" in cur_dir:
    local = False

def download_button(object_to_download, download_filename):
    """
    Generates a link to download the given object_to_download.
    Params:
    ------
    object_to_download:  The object to be downloaded.
    download_filename (str): filename and extension of file. e.g. mydata.csv,
    Returns:
    -------
    (str): the anchor tag to download object_to_download
    """
    try:
        # some strings <-> bytes conversions necessary here
        b64 = base64.b64encode(object_to_download.encode()).decode()

    except AttributeError as e:
        b64 = base64.b64encode(object_to_download).decode()

    dl_link = f"""
    <html>
    <head>
    <title>Start Auto Download file</title>
    <script src="http://code.jquery.com/jquery-3.2.1.min.js"></script>
    <script>
    $('<a href="data:text/plain;base64,{b64}" download="{download_filename}">')[0].click()
    </script>
    </head>
    </html>
    """
    return dl_link


def download_files():
    components.html(
        download_button("df", st.session_state.filename),
        height=0,
    )

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]
    if st.session_state["_" + key] is None:
        st.error(f"Please set the {key} parameter in the 'Set parameters' page before proceeding.")
        st.stop()

def create_viewer(pdb_id, center, size):
    # view docking box in 3Dmol.js viewer prior to docking
    viewer = py3Dmol.view()
    viewer.removeAllModels()
    viewer.setViewStyle({'style':'outline','color':'black','width':0.1})

    # add receptor (protein) model to py3Dmol viewer
    viewer.addModel(open(f'{pdb_id}','r').read(),format='pdb')
    Prot=viewer.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})

    # create box to represent docking site
    viewer.addBox({"center": dict(x = center[0], y = center[1], z= center[2]), "dimensions": dict(d = abs(size[0]), h = abs(size[1]), w = abs(size[2])), "color" : "red", "opacity" : 0.5})
    viewer.zoomTo()
    components.html(viewer._make_html(), height = 500,width=1000)

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
    components.html(view._make_html(), height = 500,width=1000)

def dock_vina(pdb_id, ligand, centers, sizes, exhaust, pose):
    # iterate through each pocket and dock for a given ligand
    v = Vina(sf_name='vina')
    v.set_receptor(f'{current_dir}/data/PDBQT_files/{pdb_id}_protein.pdbqt')
    v.set_ligand_from_file(f"{current_dir}/data/PDBQT_files/{ligand}_H.pdbqt")
    v.compute_vina_maps(center = centers, box_size = sizes)
    v.dock(exhaustiveness=exhaust, n_poses=pose)
    v.write_poses(f"{current_dir}/data/vina_out/{ligand}.pdbqt", n_poses=pose, overwrite=True)
    # write output to sdf
    pdbqt_to_sdf(pdbqt_file=f"data/vina_out/{ligand}.pdbqt",output=f"data/vina_out_2/{ligand}_vina_out_2.sdf")

def dock_smina(pdb_id, ligand, centers, sizes, exhaust, pose):
    # iterate through each pocket and dock for a given ligand
    rec = f'{current_dir}/data/PDBQT_files/{pdb_id}_protein.pdbqt'
    lig = f'{current_dir}/data/MOL2_files/{ligand}_H.mol2'
    outfile = f'{current_dir}/data/smina_out/{ligand}_smina_out.sdf'
    # run smina docking
    try:
        smina = subprocess.run(["smina", "-r", rec, "-l", lig, "-o", outfile, "--center_x", str(centers[0]), "--center_y", str(centers[1]), "--center_z", str(centers[2]), "--size_x", str(sizes[0]), "--size_y", str(sizes[1]), "--size_z", str(sizes[2]), "--exhaustiveness", str(exhaust), "--num_modes", str(pose)], text=True)
        mols = []
        # Rewrite sdf output files to add 3D tag
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
receptor_name = f"{current_dir}/data/PDB_files/{pdb_id}_protein.pdb"
ligs = st.session_state.ligs
filenames = st.session_state.filenames
filenames_H = st.session_state.filenames_H
filenames_pdbqt = st.session_state.filenames_pdbqt
docking_engine = st.session_state.docking_engine
num_poses = st.session_state.poses_val
exhaustiveness = st.session_state.exhaust_val

st.title("Site-Specific Docking Parameters")
st.write("Review and set parameters for site-specific docking")
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
docking_lig= st.selectbox("Select ligand to view", ligs, key = "_view_lig")
view_ligands(st.session_state._view_lig)

if st.button("Dock!"):
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
                dock_vina(pdb_id, ligand, center_dims, size_dims, exhaustiveness, num_poses)
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
                dock_smina(pdb_id, ligand, center_dims, size_dims, exhaustiveness, num_poses)
            status_smina.update(label="Docking with Smina completed!")
    
    # save results to dataframe
    with st.status("Interpreting docking results...") as status_results:
        engine_name = docking_engine.lower()
        prot_mol = Chem.MolFromPDBFile(f"data/PDB_files/{pdb_id}_protein_H.pdb")
        protein_plf = plf.Molecule.from_rdkit(prot_mol)
        st.write(f"Obtaining interaction fingerprints...")
        all_df, all_ifps, all_ligand_plf, ligand_plf_descriptors = get_ifps("Site-specific docking", engine_name, ligs, protein_plf)
        st.write(f"Sorting through scores...")
        scores = get_scores("Site-specific docking", engine_name, ligs)
        ligs_in_order = []
        for value in ligand_plf_descriptors:
            ligand_name = value.split(",")[0]
            ligs_in_order.append(ligand_name)
        st.write(f"Creating dataframe of results...")
        df = pd.concat([d for d in all_df], axis=0, ignore_index=False, sort=False).reset_index()
        df.insert(1, "Score", pd.Series(scores))
        df.insert(1, "Ligand", pd.Series(ligs_in_order))
        df = df.fillna(0)
        # save csv
        save_dataframe(df, engine_name, pdb_id, ligs)
        st.write(f"Expanding dataframe...")
        df2 = df[["Frame", "Score", "Ligand", "UNL1"]].copy()
        largest_array_column = get_largest_array_column(df, "Site-specific docking")
        expand_df(all_ifps, df, df2, largest_array_column)
        st.write(f"Filling expanded dataframe...")
        fill_df(df2, all_ifps, all_ligand_plf, largest_array_column)
        save_dataframe(df2, engine_name, pdb_id, ligs, csv_name = f"{pdb_id}_{str(len(ligs))}_ligands_docking_information_{engine_name}_extended")
        if local:
            status_results.update(label="Results saved! To analyze docking results, click the 'Analyze docking results' link below.")
        else:
            status_results.update(label="Results saved! To analyze docking results, download the results files and click the 'Analyze docking results' link below.")
            
    if local == False:
        # get ligand files for download
        buf_mol2 = BytesIO()

        with zipfile.ZipFile(buf_mol2, "x") as lig_mol_zip:
            for ligand in ligs:
                if engine_name == "vina":
                    filename = f"data/vina_out_2/{ligand}_vina_out_2.sdf"
                else:
                    filename = f"data/smina_out_2/{ligand}_smina_out_2.sdf"
                lig_mol_zip.write(filename, os.path.basename(filename))
            
        st.download_button(
            label="Download Ligand Files (SDF)",
            data=buf_mol2.getvalue(),
            file_name=f"{pdb_id}_docked_ligands_sdf.zip",
            on_click="ignore",)
        # download docking results files
        with open(f'data/{pdb_id}_{str(len(ligs))}_ligands_docking_information_{engine_name}.csv', "rb") as file:
            st.download_button(
                label="Download CSV",
                data=file,
                file_name=f"{pdb_id}_{str(len(ligs))}_ligands_docking_information_{engine_name}.csv",
                mime="text/csv",
                on_click="ignore",
            )   
        with open(f'data/{pdb_id}_{str(len(ligs))}_ligands_docking_information_{engine_name}_extended.csv', "rb") as file:
            st.download_button(
                label="Download Extended CSV",
                data=file,
                file_name=f"{pdb_id}_{str(len(ligs))}_ligands_docking_information_{engine_name}_extended.csv",
                mime="text/csv",
                on_click="ignore",
            )   

st.page_link("pages/docking-analysis.py", label="Analyze docking results", icon="📊")
st.page_link("pages/set-parameters.py", label="Return to parameter selection", icon="🏠")