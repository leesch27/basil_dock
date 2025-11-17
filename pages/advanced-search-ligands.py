import streamlit as st
import streamlit.components.v1 as components
import sys, os
import glob

import py3Dmol
import requests
from openbabel import pybel
from rcsbapi.search import AttributeQuery, Attr, TextQuery, ChemSimilarityQuery
from rdkit import Chem
from rdkit.Chem import AllChem, rdCoordGen

sys.path.insert(1, 'utilities/')
from ligandsplitter.basefunctions import create_folders

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def molecule_to_3d(molecule):
    mol = Chem.Mol(molecule)
    mol = AllChem.AddHs(mol, addCoords=True)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

def view_ligands(ligand):
    lig = f"{current_dir}/data/test_files/{ligand}_ligand.sdf"
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    view.addModel(open(lig,'r').read(),format='mol2')
    ref_m = view.getModel()
    ref_m.setStyle({},{'stick':{'colorscheme':'greenCarbon','radius':0.2}})
    view.zoomTo()
    components.html(view._make_html(), height = 500,width=500)

chem_types = ("No Selection","D-beta-peptide, C-gamma linking",
                "D-gamma-peptide, C-delta linking",
                "D-peptide COOH carboxy terminus",
                "D-peptide NH3 amino terminus",
                "D-peptide linking",
                "D-saccharide",
                "D-saccharide, alpha linking",
                "D-saccharide, beta linking",
                "DNA OH 3 prime terminus",
                "DNA OH 5 prime terminus",
                "DNA linking",
                "L-DNA linking",
                "L-RNA linking",
                "L-beta-peptide, C-gamma linking",
                "L-gamma-peptide, C-delta linking",
                "L-peptide COOH carboxy terminus",
                "L-peptide NH3 amino terminus",
                "L-peptide linking",
                "L-saccharide",
                "L-saccharide, alpha linking",
                "L-saccharide, beta linking",
                "RNA OH 3 prime terminus",
                "RNA OH 5 prime terminus",
                "RNA linking",
                "non-polymer",
                "other",
                "peptide linking",
                "peptide-like",
                "saccharide")

result_lig_list = []
if 'result_lig_list' not in st.session_state:
    st.session_state.result_lig_list = []

load_keys("local")
local = st.session_state._local

try:
    load_keys("current_dir")
    current_dir = st.session_state._current_dir
except:
    current_dir = create_folders()
    st.session_state._current_dir = current_dir

title = st.columns([0.25, 0.75])
title[0].image("img/logo.png", width=200)
title[1].title("Advanced Ligand Search using RCSB PDB")

st.write("Select attributes to search for ligands in the RCSB PDB database. At least one attribute must be selected to perform a search.")
with st.container(border=True):
    st.text_input(label="Search by Chemical Name?", placeholder='Type Chemical Name Here (e.g. alanine)',key="chem_name")
    st.text_input(label="Search by Chemical Name Synonym?", placeholder='Type Synonym Here (e.g. acetylsalicylic acid)',key="chem_syn")
    st.text_input(label="Search by Chemical ID?", placeholder='Type RCSB Chemical ID Here (e.g. AIN)',key="chem_id")
    st.selectbox("Search by Chemical Type?", chem_types, placeholder="Select chemical type", accept_new_options=False, key = "chem_type")
    st.text_input(label="Search by Chemical Brand Name?", placeholder='Type DrugBank Brand Name Here (e.g. Aspirin)',key="chem_brand")
    st.text_input(label="Search by Formula Similarity?", placeholder='Type Ligand Formula Here (e.g. C9H8O4)',key="chem_formula")
    st.text_input(label="Search by Structure Similarity?", placeholder='Type Ligand SMILES Here',key="chem_struct")
    search = st.button("Search for ligands")

if search:
    attr = [st.session_state.chem_name, st.session_state.chem_syn, st.session_state.chem_id, st.session_state.chem_type, st.session_state.chem_brand, st.session_state.chem_formula, st.session_state.chem_struct]
    attr_bool = []
    values = []
    attr_bool_dict = {}
    attr_value_dict = {}
    for index, item in enumerate(attr):
        key = f"attr{index + 1}"
        values.append(item)
        attr_value_dict[key] = item
        if item == None or item == "" or item == "No Selection":
            attr_bool.append("No")
            attr_bool_dict[key] = "No"
        else:
            attr_bool.append("Yes")
            attr_bool_dict[key] = "Yes"
    q1 = AttributeQuery(attribute = "chem_comp.name", operator = "exact_match", value = values[0])
    q2 = AttributeQuery(attribute = "rcsb_chem_comp_synonyms.name", operator = "contains_phrase", value = values[1], service = "text_chem")
    q3 = AttributeQuery(attribute = "rcsb_id", operator = "exact_match", value = values[2], service = "text_chem")
    q4 = AttributeQuery(attribute = "chem_comp.type", operator = "exact_match", value = values[3])
    q5 = AttributeQuery(attribute = "drugbank_info.brand_names", operator = "contains_phrase", value = values[4])
    q6 = ChemSimilarityQuery(query_type = "formula", value = values[5])
    q7 = ChemSimilarityQuery(query_type = "descriptor", descriptor_type = "SMILES", match_type="fingerprint-similarity", value = values[6])
    attr_list = [q1, q2, q3, q4, q5, q6, q7]
    positives = []
    global query
    for number, value in enumerate(attr_bool):
        if value == "Yes":
            positives.append(attr_list[number])
    if len(positives) > 0:
        if len(positives) == 1:
            query = positives[0]
        else:
            query = ' & '.join(x for x in positives)
    else:
        print("Invalid.")
    result_lig = list(query())
    result_lig_list_temp = []
    for nonPoly in query(return_type="mol_definition"):
        result_lig_list_temp.append(nonPoly)
    st.session_state.result_lig_list = result_lig_list_temp

# view ligands that meet criteria
st.selectbox("Select ligand to view", st.session_state.result_lig_list, key = "lig_of_interest")
if st.session_state.lig_of_interest != None and st.session_state.lig_of_interest != "":
    try: # try getting ligand as sdf file first
        ligand = st.session_state.lig_of_interest
        lig_mol2 = requests.get(f'https://files.rcsb.org/ligands/download/{ligand}_ideal.sdf')
        with open(f"data/test_files/{ligand}_ligand.sdf", "w+") as file:
            file.write(lig_mol2.text)
        lig_filename = f"data/test_files/{ligand}_ligand.sdf"
    except: # if sdf doesn't work, get ligand as cif file
        lig_mol2 = requests.get(f'https://files.rcsb.org/ligands/download/{ligand}.cif')
        with open(f"data/test_files/{ligand}_ligand.cif", "w+") as file:
            file.write(lig_mol2.text)
        lig_filename = f"data/test_files/{ligand}_ligand.cif"
        pdb_mol2 = [m for m in pybel.readfile(filename = lig_filename, format='cif')][0]
        out_mol2 = pybel.Outputfile(filename = f"data/test_files/{ligand}_ligand.sdf", overwrite = True, format='sdf')
        out_mol2.write(pdb_mol2)
    view_ligands(ligand)

if not local:
    if st.button("Prepare download"):
        pdb_mol2 = [m for m in pybel.readfile(filename = lig_filename, format='sdf')][0]
        out_mol2 = pybel.Outputfile(filename = f"data/MOL2_files/{st.session_state.lig_of_interest}.mol2", overwrite = True, format='mol2')
        out_mol2.write(pdb_mol2)
        with open(out_mol2, "r") as pdb_file:
            st.download_button(
                label="Download selected ligand as MOL2",
                data=pdb_file.read().encode("utf-8"),
                file_name=f"data/MOL2_files/{st.session_state.lig_of_interest}.mol2",
                on_click="ignore",
                mime = "application/vnd.sybyl.mol2")
else:
    # save as MOL2
    if st.button("Download selected ligand as MOL2"):
        #add garbage collection for test files?
        pdb_mol2 = [m for m in pybel.readfile(filename = lig_filename, format='sdf')][0]
        out_mol2 = pybel.Outputfile(filename = f"data/MOL2_files/{st.session_state.lig_of_interest}.mol2", overwrite = True, format='mol2')
        out_mol2.write(pdb_mol2)

    # clean test files
    if st.button("Clean test folder"):
        testing_data = os.path.join('data', 'test_files', '*')
        testing_files = glob.glob(testing_data)
        for file in testing_files:
            os.remove(file)
        st.write("Test folder cleaned.")