import streamlit as st
import streamlit.components.v1 as components
import sys, os
from io import StringIO
import glob
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, rdCoordGen
import py3Dmol

sys.path.insert(1, 'utilities/')
from ligandsplitter.ligandsplitter.basefunctions import create_folders
from ligandsplitter.ligandsplitter.ligandderive import get_func_groups, create_derivative_files

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def molecule_to_3d(molecule):
    mol = Chem.Mol(molecule)
    mol = AllChem.AddHs(mol, addCoords=True)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

def view_ligand(mol):
    view = py3Dmol.view(
    data=Chem.MolToMolBlock(mol),  # Convert the RDKit molecule for py3Dmol
    style={"stick": {}, "sphere": {"scale": 0.3}}
    )
    view.zoomTo()
    components.html(view._make_html(), height = 500,width=500)

try:
    load_keys("current_dir")
    current_dir = st.session_state.current_dir
except:
    current_dir = create_folders()
    st.session_state.current_dir = current_dir

if 'result_deriv_type_list' not in st.session_state:
    st.session_state.result_deriv_type_list = {}

if 'result_deriv_mol_list' not in st.session_state:
    st.session_state.result_deriv_mol_list = []

if 'result_deriv_smiles_list' not in st.session_state:
    st.session_state.result_deriv_smiles_list = []

if 'deriv_of_interest' not in st.session_state:
    st.session_state.deriv_of_interest = ""

st.title("Derive Ligands from Existing Ligands")
canon_ligand = st.file_uploader("Upload Ligand to Derive", accept_multiple_files= False, type="mol2", key = "canonical_ligand_upload")
if canon_ligand is not None:
    # check if canonical ligand is in MOL2 file folder
    file_location_data = os.path.join('data', 'MOL2_files', '*')
    ligands = glob.glob(file_location_data)
    if canon_ligand.name not in ligands:
        stringio = StringIO(canon_ligand.getvalue().decode("utf-8"))
        bytes_data = stringio.read()
        with open(f"data/MOL2_files/{canon_ligand.name}", "w+") as datafile:
            datafile.write(bytes_data)
    canon_lig_name = canon_ligand.name.split(".")[0]
    lig_list = [canon_lig_name]
    all_func_groups, all_derivs_mols, all_derivs_smiles, deriv_created = get_func_groups(lig_list)
    st.session_state.result_deriv_type_list = deriv_created
    st.session_state.result_deriv_mol_list = all_derivs_mols
    st.session_state.result_deriv_smiles_list = all_derivs_smiles

st.selectbox("Select derivative to view", st.session_state.result_deriv_type_list.values(), key = "deriv_of_interest")
if st.session_state.deriv_of_interest != None and st.session_state.deriv_of_interest != "":
    smiles = list(st.session_state.result_deriv_type_list.keys())[list(st.session_state.result_deriv_type_list.values()).index(st.session_state.deriv_of_interest)]
    st.write(f"Viewing derivative of {st.session_state.deriv_of_interest[2]}: {st.session_state.deriv_of_interest[0]} replaced with {st.session_state.deriv_of_interest[1]}")
    st.write(f"Derivative SMILES string: {smiles}")
    index = st.session_state.result_deriv_smiles_list.index(smiles)
    mol_of_interest = st.session_state.result_deriv_mol_list[index]
    # create personalized error message if mol_of_interest is None
    if mol_of_interest is None:
        st.write("Error: Could not create RDKit Mol object for this derivative. Cannot display 3D structure.")
    else:
        rdCoordGen.AddCoords(mol_of_interest)
        mol3d = molecule_to_3d(mol_of_interest)
        view_ligand(mol3d)
    

if st.button("Create Files for Selected Derivative of Ligand"):
    smiles = list(st.session_state.result_deriv_type_list.keys())[list(st.session_state.result_deriv_type_list.values()).index(st.session_state.deriv_of_interest)]
    index = st.session_state.result_deriv_smiles_list.index(smiles)
    selected_mol = [st.session_state.result_deriv_mol_list[index]]
    selected_smile = [st.session_state.result_deriv_smiles_list[index]]
    derivs, fxnal_groups_derivs = create_derivative_files(selected_mol, selected_smile, st.session_state.result_deriv_type_list)
    for value in derivs:
        filename = f"data/MOL2_files/{value}.mol2"
        filename_H = f"data/MOL2_files/{value}_H.mol2"
        filename_pdbqt = f"data/PDBQT_files/{value}.pdbqt"
    st.write(f"Files created for derivative: {st.session_state.deriv_of_interest}")
    st.write(f"MOL2 file: {filename}")
    st.write(f"Hydrogenated MOL2 file: {filename_H}")
    st.write(f"PDBQT file: {filename_pdbqt}")

if st.button("Create Files for All Derivatives of Ligand"):
    derivs, fxnal_groups_derivs = create_derivative_files(st.session_state.result_deriv_mol_list, st.session_state.result_deriv_smiles_list, st.session_state.result_deriv_type_list)
    filenames = []
    filenames_H = []
    filenames_pdbqt = []
    for value in derivs:
        filenames.append(f"data/MOL2_files/{value}.mol2")
        filenames_H.append(f"data/MOL2_files/{value}_H.mol2")
        filenames_pdbqt.append(f"data/PDBQT_files/{value}.pdbqt")
    st.write(f"Files created for all derivatives.")
    df = pd.DataFrame({"MOL2 files": filenames, "Hydrogenated MOL2 files": filenames_H, "PDBQT files": filenames_pdbqt})
    st.dataframe(df)