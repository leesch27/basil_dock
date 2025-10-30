import streamlit as st
import streamlit.components.v1 as components
import sys, os

from rdkit import Chem
from openbabel import pybel
import py3Dmol

sys.path.insert(1, '../utilities/')
from ligandsplitter.ligandsplitter.basefunctions import create_folders

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def view_ligands(ligand):
    lig = f"{current_dir}/data/MOL2_files/{ligand}.mol2"
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    view.addModel(open(lig,'r').read(),format='mol2')
    ref_m = view.getModel()
    ref_m.setStyle({},{'stick':{'colorscheme':'greenCarbon','radius':0.2}})
    view.zoomTo()
    components.html(view._make_html(), height = 500,width=500)

try:
    load_keys("current_dir")
    current_dir = st.session_state.current_dir
except:
    current_dir = create_folders()
    st.session_state.current_dir = current_dir

st.title("Create Ligands from SMILES Strings")
row1 = st.columns([1,1])
selected_lig = row1[0].text_input(label="Ligand Name", placeholder='Type Ligand Name Here', key="lig_name")
selected_pocket = row1[1].text_input(label="Ligand SMILES String", placeholder='Type Ligand SMILES String Here', key="lig_smiles")
if st.button("Generate Ligand from SMILES"):
    if st.session_state.lig_name and st.session_state.lig_smiles:
        lig_name = st.session_state.lig_name
        lig_scratch = st.session_state.lig_smiles
        lig_test = Chem.MolFromSmiles(lig_scratch)
        if (len(lig_scratch) < 2000) & (lig_test is not None):
            out=pybel.Outputfile(filename = f"data/MOL2_files/{lig_name}.mol2",format='mol2',overwrite=True)
            mol = pybel.readstring(string=lig_scratch,format='smiles')
            mol.title= str(lig_name)
            mol.make3D('mmff94s')
            mol.localopt(forcefield = 'mmff94s', steps = 500)
            out.write(mol)
            out.close()
            st.success(f'Ligand {st.session_state.lig_name} created successfully!')
            view_ligands(st.session_state.lig_name)
        elif lig_test is None:
            st.error("Invalid SMILES string.")
    else:
        st.error("Please provide both a ligand name and a SMILES string.")