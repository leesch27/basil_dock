import streamlit as st
import sys, os
sys.path.insert(1, 'utilities/')
from ligandsplitter.ligandsplitter.ligandgenerate import create_ligands_from_smiles, display_smiles_form, create_mols_from_smiles