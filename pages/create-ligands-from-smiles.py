import streamlit as st
import sys, os
sys.path.insert(1, '../utilities/')
from ligandsplitter.ligandsplitter.ligandgenerate import create_ligands_from_smiles, display_smiles_form, create_mols_from_smiles

attr_bool, attr_val, form_items1, form_items2 = create_search_for_expo()
display_expo_form(form_items1, form_items2)