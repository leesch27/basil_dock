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
from IPython.display import display

sys.path.insert(1, '../utilities/')
from ligandsplitter.ligandsplitter.ligandanalysis import group_idxes_from_mol, get_ligand_properties, rf_classifier, rf_regressor
from basil_utils import fetch_data_files, fetch_docking_data, get_ifps, get_scores, get_largest_array_column, expand_df, fill_df, compare_poses_form, compare_poses