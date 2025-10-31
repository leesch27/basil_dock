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
import py3Dmol
import prolif as plf
from prolif.plotting.complex3d import Complex3D
from IPython.display import display

from ligandsplitter.ligandanalysis import get_ligand_properties, oral_bioactive_classifier, interaction_regressor

def docking_data_comparison(select_type, pose_mode, pdb_id, ligand_number, select_dock, form_items1, form_items2):
    if select_type == "Blind docking":
        if pose_mode == "Compare to original pose":
            pocket_selection_list = [int(form_items1[0])]
            pose_selection_list = [int(form_items2[0])]
        else:
            pocket_selection_list = [int(form_items1[0]), int(form_items1[1])]
            pose_selection_list = [int(form_items2[0]), int(form_items2[1])]

        #initialize py3dmol viewer
        view = py3Dmol.view(height = 800, width = 900)
        view.removeAllModels()
        view.setViewStyle({'style':'outline','color':'black','width':0.1})

        # add receptor model to all py3dmol viewers
        view.addModel(open(f"data/PDB_files/{pdb_id}_protein_H.pdb",'r').read(),format='pdb')
        Prot=view.getModel()
        Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
        view.addSurface(py3Dmol.VDW,{'opacity':0.6,'color':'white'})
    
        # add reference model of ligand to py3dmol viewer
        if pose_mode == "Compare to original pose":
            view.addModel(open(f"data/MOL2_files/{ligand_number}_H.mol2",'r').read(),format='mol2')
            ref_m = view.getModel()
            ref_m.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
        else:
            selected_ref = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_pocket_{pocket_selection_list[1]}_{select_dock}_out_2.sdf")
            p_ref=Chem.MolToMolBlock(selected_ref[pose_selection_list[1]],False)
            view.addModel(p_ref,'reference mol')
            x_ref = view.getModel()
            x_ref.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
            if select_dock == 'smina':
                print('Score: {}'.format(selected_ref[pose_selection_list[1]].GetProp('minimizedAffinity')))
            else:
                print('Pose: {} | Score: {}'.format(selected_ref[pose_selection_list[1]].GetProp('Pose'), selected_ref[pose_selection_list[1]].GetProp('Score')))
    
        # add experimental docking data of a desired pose in a pocket to py3dmol viewer
        selected = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_pocket_{pocket_selection_list[0]}_{select_dock}_out_2.sdf")
        p=Chem.MolToMolBlock(selected[pose_selection_list[0]],False)
        if select_dock == 'smina':
            print('Smina Pose (' + str(ligand_number) + '): Cyan')
            print('Score: {}'.format(selected[pose_selection_list[0]].GetProp('minimizedAffinity')))
        else:
            print('Vina Pose (' + str(ligand_number) + '): Cyan')
            print('Pose: {} | Score: {}'.format(selected[pose_selection_list[0]].GetProp('Pose'), selected[pose_selection_list[0]].GetProp('Score')))
        view.addModel(p,'mol')
        x = view.getModel()
        x.setStyle({},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})
    else:
        if pose_mode == "Compare to original pose":
            pose_selection_list = [int(form_items2[0])]
        else:
            pose_selection_list = [int(form_items2[0]), int(form_items2[1])]
    
        #initialize py3dmol viewer
        view = py3Dmol.view(height = 800, width = 900)
        view.removeAllModels()
        view.setViewStyle({'style':'outline','color':'black','width':0.1})

        # add receptor model to all py3dmol viewers
        view.addModel(open(f"data/PDB_files/{pdb_id}_protein_H.pdb",'r').read(),format='pdb')
        Prot=view.getModel()
        Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
        view.addSurface(py3Dmol.VDW,{'opacity':0.6,'color':'white'})
    
        # add reference model of ligand to py3dmol viewer
        if pose_mode == "Compare to original pose":
            view.addModel(open(f"data/MOL2_files/{ligand_number}_H.mol2",'r').read(),format='mol2')
            ref_m = view.getModel()
            ref_m.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
        else:
            selected_ref = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_{select_dock}_out_2.sdf")
            p_ref=Chem.MolToMolBlock(selected_ref[pose_selection_list[1]],False)
            view.addModel(p_ref,'reference mol')
            x_ref = view.getModel()
            x_ref.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
            if select_dock == 'smina':
                print('Score: {}'.format(selected_ref[pose_selection_list[1]].GetProp('minimizedAffinity')))
            else:
                print('Pose: {} | Score: {}'.format(selected_ref[pose_selection_list[1]].GetProp('Pose'), selected_ref[pose_selection_list[1]].GetProp('Score')))
    
        # add experimental docking data of a desired pose in a pocket to py3dmol viewer
        selected = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_{select_dock}_out_2.sdf")
        p=Chem.MolToMolBlock(selected[pose_selection_list[0]],False)
        if select_dock == 'smina':
            print('Smina Pose (' + str(ligand_number) + '): Cyan')
            print ('Score: {}'.format(selected[pose_selection_list[0]].GetProp('minimizedAffinity')))
        else:
            print('Vina Pose (' + str(ligand_number) + '): Cyan')
            print ('Pose: {} | Score: {}'.format(selected[pose_selection_list[0]].GetProp('Pose'), selected[pose_selection_list[0]].GetProp('Score')))
        view.addModel(p,'mol')
        x = view.getModel()
        x.setStyle({},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})
        
    view.zoomTo()
    view.show()
    components.html(view._make_html(), height = 500,width=1000)

def create_plf_viewer(ifp, plf_ligand, plf_prot):
    # view docking box in 3Dmol.js viewer prior to docking
    comp = Complex3D(ifp, plf_ligand, plf_prot)
    components.html(comp.display()._make_html(), height = 500,width=1000)

st.title("Docking Analysis")
st.write("Compare poses between ligand conformations")
data_retrieval = st.file_uploader("Upload Results", accept_multiple_files= False, type="csv", key = "docking_data_upload")
#initialize session state variables
if 'pdb_id' not in st.session_state:
    st.session_state.pdb_id = ""
if 'ligs' not in st.session_state:
    st.session_state.ligs = []
if "poses" not in st.session_state:
    st.session_state.poses = []
if "pockets" not in st.session_state:
    st.session_state.pockets = []
if "select_dock" not in st.session_state:
    st.session_state.select_dock = ""
if "select_type" not in st.session_state:
    st.session_state.select_type = ""
if "selected_ligand" not in st.session_state:
    st.session_state.selected_ligand = ""
if "selected_pocket" not in st.session_state:
    st.session_state.selected_pocket = ""
if "selected_pose" not in st.session_state:
    st.session_state.selected_pose = ""
if "selected_pocket_2" not in st.session_state:
    st.session_state.selected_pocket_2 = ""
if "selected_pose_2" not in st.session_state:
    st.session_state.selected_pose_2 = ""
if "form_items_1" not in st.session_state:
    st.session_state.form_items_1 = []
if "form_items_2" not in st.session_state:
    st.session_state.form_items_2 = []

if data_retrieval is not None:
    csv_name = data_retrieval.name
    pdb_id = csv_name.split("_")[0]
    st.session_state.pdb_id = pdb_id

    df = pd.read_csv(data_retrieval)
    st.dataframe(df)
    if "Pocket" in df.columns:
        st.session_state.select_type = "Blind docking"
        pockets = df['Pocket'].unique().tolist()
        st.session_state.pockets = pockets
    else:
        st.session_state.select_type = "Site-specific docking"
    if "vina" in csv_name.lower():
        st.session_state.select_dock = "vina"
    else:
        st.session_state.select_dock = "smina"
    ligs = df['Ligand'].unique().tolist()
    ligs_clean = []
    for lig in ligs:
        lig_name = lig.split(" ")[-1]
        ligs_clean.append(lig_name)
    st.session_state.ligs = ligs_clean
    poses = df['Frame'].unique().tolist()
    #poses = df.index.unique().tolist()
    st.session_state.poses = poses

    comp_method = st.selectbox("Select Comparison Method", ["Compare to original pose", "Compare to other docked poses"], key = "method_selection")
    if st.session_state.method_selection == "Compare to original pose":
        if st.session_state.select_type == "Blind docking":
            row1 = st.columns([1,1,1])
            selected_lig = row1[0].selectbox("Select Ligand to View", st.session_state.ligs, key = "selected_ligand")
            selected_pocket = row1[1].selectbox("Select Pocket to View", st.session_state.pockets, key = "selected_pocket")
            selected_pose = row1[2].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose")
            st.session_state.form_items_1 = [selected_pocket]
            st.session_state.form_items_2 = [selected_pose]

        else:
            row1 = st.columns([1,1])
            selected_lig = row1[0].selectbox("Select Ligand to View", st.session_state.ligs, key = "selected_ligand")
            selected_pose = row1[1].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose")
            st.session_state.form_items_2 = [selected_pose]
    else:
        if st.session_state.select_type == "Blind docking":
            row1 = st.columns([1,1,1])
            selected_lig = row1[0].selectbox("Select Ligand to View", st.session_state.ligs, key = "selected_ligand")
            selected_pocket = row1[1].selectbox("Select Pocket to View", st.session_state.pockets, key = "selected_pocket")
            selected_pose = row1[2].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose")

            row2 = st.columns([1,1,1])
            selected_pocket_2 = row2[1].selectbox("Select Pocket to View", st.session_state.pockets, key = "selected_pocket_2")
            selected_pose_2 = row2[2].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose_2")
            st.session_state.form_items_1 = [selected_pocket, selected_pocket_2]
            st.session_state.form_items_2 = [selected_pose, selected_pose_2]
        else:
            row1 = st.columns([1,1])
            selected_lig = row1[0].selectbox("Select Ligand to View", st.session_state.ligs, key = "selected_ligand")
            selected_pose = row1[1].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose")

            row2 = st.columns([1,1])
            selected_pose_2 = row2[1].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose_2")
            st.session_state.form_items_2 = [selected_pose, selected_pose_2]

    if st.button("Compare Docking Data"):
        form_items1 = st.session_state.form_items_1
        form_items2 = st.session_state.form_items_2
        docking_data_comparison(st.session_state.select_type, st.session_state.method_selection, pdb_id, st.session_state.selected_ligand, st.session_state.select_dock, form_items1, form_items2)

    st.write("View interactions between ligand and protein receptor")
    if st.session_state.select_type == "Blind docking":
        row1 = st.columns([1,1,1])
        selected_lig = row1[0].selectbox("Select Ligand to View", st.session_state.ligs, key = "selected_ligand_interaction")
        selected_pocket = row1[1].selectbox("Select Pocket to View", st.session_state.pockets, key = "selected_pocket_interaction")
        selected_pose = row1[2].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose_interaction")
        pose_num = int(st.session_state.selected_pose_interaction) + 1

    else:
        row1 = st.columns([1,1])
        selected_lig = row1[0].selectbox("Select Ligand to View", st.session_state.ligs, key = "selected_ligand_interaction")
        selected_pose = row1[1].selectbox("Select Pose to View", st.session_state.poses, key = "selected_pose_interaction")
        pose_num = int(st.session_state.selected_pose_interaction) + 1

    if st.button("View Interactions"):
        # find way to get ifps
        prot_mol = Chem.MolFromPDBFile(f"data/PDB_files/{pdb_id}_protein_H.pdb")
        protein_plf = plf.Molecule.from_rdkit(prot_mol)
        if st.session_state.select_type == "Blind docking":
            lig_suppl = plf.sdf_supplier(f"data/{st.session_state.select_dock}_out_2/{st.session_state.selected_ligand_interaction}_pocket_{st.session_state.selected_pocket_interaction}_{st.session_state.select_dock}_out_2.sdf")
            fp = plf.Fingerprint(count=True)
            fp.run_from_iterable(lig_suppl, protein_plf)
            ifp = fp.generate(lig_suppl[int(st.session_state.selected_pose_interaction)], protein_plf, metadata = True)
        else:
            lig_suppl = plf.sdf_supplier(f"data/{st.session_state.select_dock}_out_2/{st.session_state.selected_ligand_interaction}_{st.session_state.select_dock}_out_2.sdf")
            fp = plf.Fingerprint(count=True)
            fp.run_from_iterable(lig_suppl, protein_plf)
            ifp = fp.generate(lig_suppl[int(st.session_state.selected_pose_interaction)], protein_plf, metadata = True)
        lig_plf = lig_suppl[int(st.session_state.selected_pose_interaction)]
        create_plf_viewer(ifp, lig_plf, protein_plf)


    st.write("Determine which features affect ligand binding affinity the strongest")
    if st.button("Calculate Important Binding Features"):
        with st.status("Analyzing feature importances...") as status_regression:
            if df.shape[0] < 15:
                rf_affinity_importances = None
                xgb_affinity_importances = None
                status_regression.update(label="Not enough data points to perform regression analysis. Need at least 15 data points.")
            else:
                rf_affinity_importances, xgb_affinity_importances = interaction_regressor(df)
                status_regression.update(label="Analysis completed!")
        if(rf_affinity_importances is not None) and (xgb_affinity_importances is not None):
            st.write("Feature importances from Random Forest Regressor")
            st.dataframe(rf_affinity_importances, key = "rf_affinity_importances")
            st.write("Feature importances from XGBoost Regressor")
            st.dataframe(xgb_affinity_importances, key = "xgb_affinity_importances")

    st.write("Determine if features and properties suggest potential oral bioavailability of ligands")
    bioactive_df = pd.read_csv("data/test_train_bioactive_data.csv")
    ligand_smiles = pd.read_csv(f'data/ligand_smiles_data_id_{st.session_state.pdb_id}_{str(len(st.session_state.ligs))}.csv')
    try:
        bioactive_df = pd.concat([bioactive_df, ligand_smiles])
    except:
        bioactive_df = pd.read_csv("data/test_train_bioactive_data.csv")

    bioactive_method = st.selectbox("Select Method for Analysis:", ["Choose Method", "LRO5", "Ghose", "Veber"], key = "bioactive_method")
    if st.button("Calculate Important Bioavailability Features"):
        if bioactive_method != "Choose Method":
            with st.status("Analyzing feature importances...") as status_classification:
                lig_bioactive_df = get_ligand_properties(bioactive_df)
                bioactive_importances, class_labels = oral_bioactive_classifier(lig_bioactive_df, bioactive_method)
                for key, value in class_labels.items():
                    if value == 1:
                        st.write(f"Predicted orally bioactive value for {key}: Yes")
                    else:
                        st.write(f"Predicted orally bioactive value for {key}: No")
                st.write("Feature importances for oral bioavailability classification using " + bioactive_method + " method")
                status_classification.update(label="Analysis completed!")
        else:
            bioactive_importances = None
            st.write("Please select a method for analysis.")
        if bioactive_importances is not None:
            st.dataframe(bioactive_importances, key = "bioactive_importances")
            # write results to page