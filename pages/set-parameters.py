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

from Bio.PDB import PDBList
import MDAnalysis as mda 
from MDAnalysis.coordinates import PDB
from openbabel import pybel
from rdkit import Chem

from ligandsplitter.basefunctions import create_folders
from ligandsplitter.ligandsplit import retrieve_pdb_file

def save_keys(key):
    st.session_state[key] = st.session_state["_" + key]

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

cur_dir = os.getcwd()
local = True
if "mount/src" in cur_dir:
    local = False

print(f"Local mode: {local}")
with st.form("enter_docking_parameters"):
    header = st.columns([1,1])
    header[0].subheader("basil dock Settings")
    header[1].subheader("Receptor and Ligand Selection Settings")

    row1 = st.columns([1,1])
    docking_engine = row1[0].selectbox("Select Docking Engine to Use", ["Vina", "Smina"], key = "_docking_engine")
    protein_retrieval = row1[1].file_uploader("Upload Protein Receptor To Use", accept_multiple_files= False, type=["cif", "pdb", "ent"], key = "_protein_upload")
    
    row2 = st.columns([1,1])
    docking_method = row2[0].selectbox("Select Docking Method to Use", ["Blind Docking", "Site-Specific Docking"], key = "docking_method")
    ligand_retrieval = row2[1].file_uploader("Upload Ligand/s To Use", accept_multiple_files= True, type=["mol2", "sdf"], key = "_ligand_upload")

    row3 = st.columns([1,1])
    exhaustiveness = row3[0].slider("Select Exhaustiveness Value", 1, 16, 5, key = "_exhaust_val")

    row4 = st.columns([1,1])
    poses = row4[0].slider("Select Number of Poses to Generate", 1, 10, 5, key = "_poses_val")

    row5 = st.columns([1,1])
    submitted = row5[1].form_submit_button("Submit Parameters")
    
    if submitted:
        ligs = []
        filenames = []
        protein = st.session_state._protein_upload
        ligands = st.session_state._ligand_upload

        if protein is None:
            st.error("Please upload a protein receptor to proceed.")
            st.stop()
        if len(ligands) == 0:
            st.error("Please upload at least one ligand to proceed.")
            st.stop()

        current_dir = create_folders()
        st.session_state._current_dir = current_dir

        # get receptor data
        with st.status("Retrieving protein and ligands...") as status:
            if protein is not None:
                st.write("Isolating protein from ligands...")
                receptor_name = protein.name
                pdb_id = receptor_name.split(".")[0]
                stringio = StringIO(protein.getvalue().decode("utf-8"))
                bytes_data = stringio.read()
                protein_filename = f"data/PDB_files/{receptor_name}"
                with open(protein_filename, "w+") as datafile:
                    datafile.write(bytes_data)
                retrieve_pdb_file(protein_filename, "local")
            
                # protein sanitization
                # create pdbqt file for receptor
                st.write("Sanitizing protein...")
                input_file = f"data/PDB_files/{pdb_id}_protein.pdb"
                pqr_file = f"data/PDB_files/{pdb_id}_protein.pqr"
                output_file = f"data/PDB_files/{pdb_id}_protein_H.pdb"
                try:
                    if local:
                        pqr = subprocess.run(["pdb2pqr", f"--pdb-output={output_file}", "--pH=7.4", input_file, pqr_file, "--whitespace", "--quiet"])
                    else:
                        with open(output_file, "w+") as out_file:
                            pqr = subprocess.run([f"{sys.executable}", "pdb2pqr", f"--pdb-output={output_file}", "--pH=7.4", input_file, pqr_file, "--whitespace", "--quiet"], text= True, check=True, stdout=out_file)
                            with open(output_file, "r") as out_check:
                                for line in out_check:
                                    print(line)
                            with open(pqr_file, "r") as pqr_check:
                                for line in pqr_check:
                                    print(line)
                    to_pdbqt = mda.Universe(pqr_file)
                    to_pdbqt.atoms.write(f"data/PDBQT_files/{pdb_id}_protein.pdbqt")

                    # remove "TITLE" and "CRYST1" labels with "REMARK" to reduce chance of errors later on
                    with open(f"data/PDBQT_files/{pdb_id}_protein.pdbqt", 'r') as file:
                        file_content = file.read()
                    file_content = file_content.replace('TITLE', 'REMARK').replace('CRYST1', 'REMARK')
                    with open(f"data/PDBQT_files/{pdb_id}_protein.pdbqt", 'w') as file:
                        file.write(file_content)
                except subprocess.CalledProcessError as pqr:
                    print(pqr.stderr, end="")
                    status.update(label="Errors encountered during protein sanitization! Please check log file for details.")
                st.session_state._pdb_id = pdb_id
        
            # get ligand data
            st.write("Writing ligand files...")
            for uploaded_file in ligands:
                ligand_name = uploaded_file.name
                ligand_id_proto = ligand_name.split(".")[:-1]
                ligand_id = "".join(str(x) for x in ligand_id_proto)
                ligand_extension = ligand_name.split(".")[-1]
                stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
                bytes_data = stringio.read()
                with open(f"data/MOL2_files/{ligand_name}", "w+") as datafile:
                    datafile.write(bytes_data)
                if ligand_extension == "sdf":
                    pdb_mol2 = [m for m in pybel.readfile(filename = f"data/MOL2_files/{ligand_name}", format='sdf')][0]
                    out_mol2 = pybel.Outputfile(filename = f"data/MOL2_files/{ligand_id}.mol2", overwrite = True, format='mol2')
                    out_mol2.write(pdb_mol2)
                    ligs.append(ligand_id)
                    filenames.append(f"data/MOL2_files/{ligand_id}.mol2")
                else:
                    ligs.append(ligand_id)
                    filenames.append(f"data/MOL2_files/{ligand_id}.mol2")
            # ligand sanitization
            # add hydrogens to ligands
            st.write("Sanitizing ligands...")
            filenames_H = []
            a = 0
            for i in filenames:
                mol= [m for m in pybel.readfile(filename= str(i),format='mol2')][0]
                mol.addh()
                s = f"data/MOL2_files/{ligs[a]}_H.mol2"
                filenames_H.append(s)
                out = pybel.Outputfile(filename= s,format='mol2',overwrite=True)
                out.write(mol)
                out.close()
                a += 1
            # ligand sanitization
            # convert to pdbqt
            n = 0
            filenames_pdbqt = []
            for i in filenames:
                ligand = [m for m in pybel.readfile(filename= str(i) ,format='mol2')][0]
                s = f"data/PDBQT_files/{ligs[n]}_H.pdbqt"
                filenames_pdbqt.append(s)
                ligand.write(filename = s, format='pdbqt', overwrite=True)
                n += 1
            
            # generate smiles string for each ligand
            ligand_smiles = []
            for i in ligs:
                mol = Chem.MolFromMol2File("data/MOL2_files/" + str(i) + "_H.mol2",sanitize=False)
                select_mol_smile = Chem.MolToSmiles(mol)
                ligand_smiles.append(select_mol_smile)
            # save ligand data to csv
            ligand_smiles_data = pd.DataFrame({"filename_hydrogens": filenames_H, "smiles": ligand_smiles})
            ligand_smiles_data.to_csv(f'data/ligand_smiles_data_id_{pdb_id}_{str(len(ligs))}.csv', index = False)
            
            st.session_state._ligs = ligs
            st.session_state._filenames = filenames
            st.session_state._filenames_H = filenames_H
            st.session_state._filenames_pdbqt = filenames_pdbqt
            status.update(label="Seperating and sanitizing of molecules completed!")
        
        dock_method = st.session_state.docking_method
        save_keys("current_dir")
        save_keys("protein_upload")
        save_keys("ligand_upload")
        save_keys("docking_engine")
        save_keys("exhaust_val")
        save_keys("poses_val")
        save_keys("ligs")
        save_keys("filenames")
        save_keys("filenames_H")
        save_keys("filenames_pdbqt")
        save_keys("pdb_id")
        if local:
            if dock_method == "Blind Docking":
                st.switch_page("pages/blind-docking.py")
            else:
                st.switch_page("pages/site-specific-docking.py")
if not local and "_pdb_id" in st.session_state:
    st.write("Make sure to download your receptor and ligand files before proceeding!")
    filenames = st.session_state._filenames
    filenames_H = st.session_state._filenames_H
    filenames_pdbqt = st.session_state._filenames_pdbqt
    pdb_id = st.session_state._pdb_id

    buf_mol2 = BytesIO()
    buf_mol2_H = BytesIO()
    buf_pdbqt = BytesIO()
    #f"data/PDB_files/{pdb_id}_protein.pdb"
    #f"data/PDB_files/{pdb_id}_protein_H.pdb"
    #f"data/PDB_files/{pdb_id}_protein_H.pdbqt"
    pdb_file = open(f"data/PDB_files/{pdb_id}_protein.pdb", "r", encoding="utf-8")
    pdb_H_file = open(f"data/PDB_files/{pdb_id}_protein_H.pdb", "r", encoding="utf-8")
    pdbqt_file = open(f"data/PDBQT_files/{pdb_id}_protein.pdbqt", "r", encoding="utf-8")

    with zipfile.ZipFile(buf_mol2, "x") as lig_mol_zip:
        for item in filenames:
            lig_mol_zip.write(item, os.path.basename(item))

    with zipfile.ZipFile(buf_mol2_H, "x") as lig_mol_zip_H:
        for item in filenames_H:
            lig_mol_zip_H.write(item, os.path.basename(item))

    with zipfile.ZipFile(buf_pdbqt, "x") as lig_mol_zip:
        for item in filenames_pdbqt:
            lig_mol_zip.write(item, os.path.basename(item))

    st.download_button(
        label="Download Sanitized Receptor (PDB)",
        data=pdb_file.read(),
        file_name=f"{pdb_id}_protein.pdb",)
    st.download_button(
        label="Download Sanitized Receptor (PDB, Protonated)",
        data=pdb_H_file.read(),
        file_name=f"{pdb_id}_protein_H.pdb",)
    st.download_button(
        label="Download Sanitized Receptor (PDBQT)",
        data=pdbqt_file.read(),
        file_name=f"{pdb_id}_protein.pdbqt",)
    st.download_button(
        label="Download Ligand Files (MOL2)",
        data=buf_mol2.getvalue(),
        file_name=f"{pdb_id}_ligands_mol2.zip",)
    st.download_button(
        label="Download Ligand Files (MOL2, Protonated)",
        data=buf_mol2_H.getvalue(),
        file_name=f"{pdb_id}_ligands_H.mol2.zip",)
    st.download_button(
        label="Download Ligand Files (PDBQT)",
        data=buf_pdbqt.getvalue(),
        file_name=f"{pdb_id}_ligands_pdbqt.zip",)
    st.download_button(
        label="Download Ligand SMILES Information (CSV)",
        data=ligand_smiles_data.to_csv().encode("utf-8"),
        file_name=f"ligand_smiles_data_id_{pdb_id}_{str(len(ligs))}.csv",)
    if st.button("Proceed to Docking Pages"):
        pdb_file.close()
        pdb_H_file.close()
        pdbqt_file.close()
        if dock_method == "Blind Docking":
            st.switch_page("pages/blind-docking.py")
        else:
            st.switch_page("pages/site-specific-docking.py")
        
            

