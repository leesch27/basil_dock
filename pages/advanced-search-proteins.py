#protein advanced search
import streamlit as st
import streamlit.components.v1 as components
import sys, os
import glob
import pandas as pd

from Bio.PDB import PDBList
from rcsbapi.search import AttributeQuery, Attr, TextQuery

import py3Dmol

from ligandsplitter.basefunctions import create_folders

def load_keys(key):
    st.session_state["_" + key] = st.session_state[key]

def view_prot(prot):
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    view.addModel(open(f'data/test_files/pdb{prot}.ent','r').read(),format='pdb')
    Prot=view.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
    view.zoomTo()
    view.show()
    components.html(view._make_html(), height = 500,width=500)

if 'result_prot_list' not in st.session_state:
    st.session_state.result_prot_list = []

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
title[1].title("Advanced Receptor Search using RCSB PDB")

# define advanced search attributes
comp_operators = ["==", ">", ">=", "<", "<="]

# get classification names and numbers
enzyme_class_df = pd.read_csv("data/enzyme_classification.csv")
enzyme_class_numbers = ["No Selection"] + enzyme_class_df.iloc[:, 0].tolist()
enzyme_class_names = ["No Selection"] + enzyme_class_df.iloc[:, 1].tolist()

#print(enzyme_class_names[:5]) TEST TEST

st.write("Select attributes to search for protein receptors in the RCSB PDB database. At least one attribute must be selected to perform a search.")
with st.container(border=True):
    row1 = st.columns([2,1,2])
    selected_lig = row1[0].write("Search for Proteins by Enzyme Classification Name")
    selected_pocket = row1[1].selectbox("Select Operator", ["is", "is not empty"], key = "class_name_operator")
    selected_pose = row1[2].selectbox("Search by Enzyme Classification Name?", enzyme_class_names, accept_new_options=False, key="class_name")

    row2 = st.columns([2,1,2])
    selected_lig = row2[0].write("Search for Proteins by Enzyme Classification Number")
    selected_pocket = row2[1].selectbox("Select Operator", ["is any of", "is not empty"], key = "class_number_operator")
    selected_pose = row2[2].selectbox("Search by Enzyme Classification Number?", enzyme_class_numbers, accept_new_options=False, key="class_number")

    row3 = st.columns([2,1,2])
    selected_lig = row3[0].write("Search for Proteins by Number of Chains")
    selected_pocket = row3[1].selectbox("Select Operator", comp_operators, key = "num_chain_operator")
    selected_pose = row3[2].text_input(label="Search by Number of Chains?", placeholder='Type Number of Protein Chains Here', label_visibility = "hidden", key="num_chain")

    row4 = st.columns([2,1,2])
    selected_lig = row4[0].write("Search for Proteins by Length of Sequence")
    selected_pocket = row4[1].selectbox("Select Operator", comp_operators, key = "length_seq_operator")
    selected_pose = row4[2].text_input(label="Search by Length of Sequence?", placeholder='Type Length of Sequence Here', label_visibility = "hidden", key="length_seq")

    row5 = st.columns([2,1,2])
    selected_lig = row5[0].write("Search for Proteins by Molecular Weight")
    selected_pocket = row5[1].selectbox("Select Operator", comp_operators, key = "molec_weight_operator")
    selected_pose = row5[2].text_input(label="Search by Molecular Weight?", placeholder='Type Molecular Weight Here', label_visibility = "hidden", key="molec_weight")

    search = st.button("Search for proteins")

if search:
    try:
        attr = [st.session_state.class_name, st.session_state.class_number, st.session_state.num_chain, st.session_state.length_seq, st.session_state.molec_weight]
        attr_comp = [st.session_state.class_name_operator, st.session_state.class_number_operator, st.session_state.num_chain_operator, st.session_state.length_seq_operator, st.session_state.molec_weight_operator]
        attr_bool = []
        values = []
        attr_bool_dict = {}
        attr_value_dict = {}
        comp_vals = []
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
    
        #translate initial comparison symbols into terms needed for search
        for value in attr_comp:
            if value == "==":
                comp_vals.append("equals")
            elif value == ">":
                comp_vals.append("greater")
            elif value == ">=":
                comp_vals.append("greater_or_equal")
            elif value == "<":
                comp_vals.append("less")
            elif value == "<=":
                comp_vals.append("less_or_equal")
            elif value == "is":
                comp_vals.append("exact_match")
            elif value == "is not empty":
                comp_vals.append("exists")
            elif value == "is any of":
                comp_vals.append("in")
            else:
                comp_vals.append("equals")
    
        q0 = AttributeQuery(attribute = "rcsb_entry_info.selected_polymer_entity_types", operator = "exact_match", value = "Protein (only)")
        # Search for Proteins by Enzyme Classification Name
        # LEE NOTE TO LEE: appears to be controlled vocabulary. possible to make dropdown menu instead of text input?
        # LEE UPDATE: using dropdown menu now with BRENDA database names for both ec name and number. need to make sure names/numbers match RCSB PDB database
        if attr_bool[0] == "Yes":
            if comp_vals[0] == "exists":
                q1 = AttributeQuery(attribute = "rcsb_polymer_entity.rcsb_ec_lineage.name", operator = "exists")
            else:
                q1 = AttributeQuery(attribute = "rcsb_polymer_entity.rcsb_ec_lineage.name", operator = comp_vals[0], value = values[0])
        else:
            q1 = AttributeQuery(attribute = "rcsb_polymer_entity.rcsb_ec_lineage.name", operator = "equals", value = "")
        
        # Search for Proteins by Enzyme Classification Number
        if attr_bool[1] == "Yes":
            if comp_vals[1] == "exists":
                q2 = AttributeQuery(attribute = "rcsb_polymer_entity.rcsb_ec_lineage.id", operator = "exists")
            else:
                value_list = [x.strip() for x in values[1].split(",")]
                q2 = AttributeQuery(attribute = "rcsb_polymer_entity.rcsb_ec_lineage.id", operator = comp_vals[1], value = value_list)
        else:
            q2 = AttributeQuery(attribute = "rcsb_polymer_entity.rcsb_ec_lineage.id", operator = "equals", value = "")
        
        # Search for Proteins by Number of Chains
        if attr_bool[2] == "Yes":
            q3 = AttributeQuery(attribute = "rcsb_entry_info.polymer_entity_count_protein", operator = comp_vals[2], value = int(values[2]))
        else:
            q3 = AttributeQuery(attribute = "rcsb_entry_info.polymer_entity_count_protein", operator = "equals", value = 0)
        
        # Search for Proteins by Length of Sequence
        if attr_bool[3] == "Yes":
            q4 = AttributeQuery(attribute = "entity_poly.rcsb_sample_sequence_length", operator = comp_vals[3], value = int(values[3]))
        else:
            q4 = AttributeQuery(attribute = "entity_poly.rcsb_sample_sequence_length", operator = "equals", value = 0)
        
        # Search for Proteins by Molecular Weight
        if attr_bool[4] == "Yes":
            q5 = AttributeQuery(attribute = "rcsb_entry_info.molecular_weight", operator = comp_vals[4], value = float(values[4]))
        else:
            q5 = AttributeQuery(attribute = "rcsb_entry_info.molecular_weight", operator = "equals", value = 0)
        
        attr_list = [q1, q2, q3, q4, q5]
        positives = [q0]
        global query
        for number, value in enumerate(attr_bool):
            if value == "Yes":
                positives.append(attr_list[number])
        if len(positives) > 0:
            if len(positives) == 1:
                query = positives[0]
            else:
                current_len = 1
                query = positives[0]
                while current_len < len(positives):
                    query = query & positives[current_len]
                    current_len += 1
        else:
            print("Invalid.")
        result_prot = list(query())
        st.session_state.result_prot_list = result_prot
    except:
        # need to add personalized error message
        pass

# view proteins that meet criteria
st.selectbox("Select protein to view", st.session_state.result_prot_list, key = "prot_of_interest")
if st.session_state.prot_of_interest != None and st.session_state.prot_of_interest != "":
    pdb_list = PDBList()
    pdb_prot_filename = pdb_list.retrieve_pdb_file(st.session_state.prot_of_interest, pdir="data/test_files", file_format="pdb")
    lower_id = st.session_state.prot_of_interest.lower()
    view_prot(lower_id)

if not local:
    if st.button("Prepare download"):
        pdb_list = PDBList()
        pdb_lower = st.session_state.prot_of_interest.lower()
        pdb_prot_filename = pdb_list.retrieve_pdb_file(st.session_state.prot_of_interest, pdir="data/PDB_files", file_format="pdb")
        cif_prot_filename = pdb_list.retrieve_pdb_file(st.session_state.prot_of_interest, pdir="data/PDB_files", file_format="mmCif")
        with open(pdb_prot_filename, "r") as pdb_file:
            st.download_button(
                label="Download selected receptor as PDB/ENT",
                data=pdb_file.read().encode("utf-8"),
                file_name=f"{pdb_lower}.ent",
                on_click="ignore",
                mime="text/plain",)
        with open(cif_prot_filename, "r") as cif_file:
            st.download_button(
                label="Download selected receptor as CIF",
                data=cif_file.read().encode("utf-8"),
                file_name=f"{pdb_lower}.cif",
                on_click="ignore",
                mime="text/plain",)
else:
    # save as PDB/ENT
    if st.button("Download selected receptor as PDB/ENT"):
        #add garbage collection for test files?
        pdb_lower = st.session_state.prot_of_interest.lower()
        file_location_data = os.path.join('data', 'PDB_files', '*.ent')
        receptors = glob.glob(file_location_data)
        if f"pdb{pdb_lower}.ent" in receptors:
            st.write("PDB file already exists in PDB_files folder.")
        else:
            pdb_list = PDBList()
            pdb_prot_filename = pdb_list.retrieve_pdb_file(st.session_state.prot_of_interest, pdir="data/PDB_files", file_format="pdb")
            st.write("PDB file saved to data/PDB_files folder.")

    # save as mmCif
    if st.button("Download selected receptor as CIF"):
        pdb_lower = st.session_state.prot_of_interest.lower()
        file_location_data = os.path.join('data', 'PDB_files', '*.cif')
        receptors = glob.glob(file_location_data)
        if f"{pdb_lower}.cif" in receptors:
            st.write("CIF file already exists in PDB_files folder.")
        else:
            pdb_list = PDBList()
            pdb_prot_filename = pdb_list.retrieve_pdb_file(st.session_state.prot_of_interest, pdir="data/PDB_files", file_format="mmCif")
            st.write("CIF file saved to data/PDB_files folder.")

    # clean test files
    if st.button("Clean test folder"):
        testing_data = os.path.join('data', 'test_files', '*')
        testing_files = glob.glob(testing_data)
        for file in testing_files:
            os.remove(file)
        st.write("Test folder cleaned.")