import streamlit as st
import sys, os
sys.path.insert(1, 'utilities/')
from ligandsplitter.ligandsplitter.ligandanalysis import group_idxes_from_mol, get_ligand_properties, rf_classifier, rf_regressor
from ligandsplitter.ligandsplitter.ligandderive import get_func_groups, create_derivative_files