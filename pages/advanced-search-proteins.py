#protein advanced search
import streamlit as st
import sys, os

import rcsbapi
from rcsbapi.search import AttributeQuery, Attr, TextQuery, ChemSimilarityQuery
sys.path.insert(1, 'utilities/')
from ligandsplitter.ligandsplitter.ligandgenerate import create_search_for_protein, display_expo_form, create_proteins