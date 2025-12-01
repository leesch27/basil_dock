import streamlit as st

title = st.columns([0.25, 0.75])
title[0].image("img/logo.png", width=200)
title[1].title("Basic Concepts in Molecular Docking")
st.header("Overview", divider="gray")

st.header("Steps involved in Molecular Docking", divider="gray")
st.subheader("Obtaining the Receptor and Ligand Structure Files")
st.markdown('''
    When performing molecular docking simulations, the first step is to obtain the structure files for both the receptor (in this case, a protein) and the ligand (a small molecule that has the potential to bind to the receptor).
    These structure files can be obtained from various sources, such as the RCSB Protein Data Bank (PDB) for protein structures and either the RCSB PDB or the PubChem database for small molecule structures.
    Receptors are usually stored in PDB or CIF format, and ligands are typically stored in MOL2 or SDF format. Alternatively, users can create ligand representations using SMILES strings, which encode ligands as text strings.''')

st.subheader("Sanitizing Structure Files")
st.markdown('''
    Before performing molecular docking simulations, it is important to sanitize the receptor and ligand structure files to ensure they 
    contain all necessary information. This involves adding missing hydrogen atoms, assigning proper bond orders, and removing any unwanted molecules such as water or ions from the receptor structure.
    Structures may also need to be converted to a different format depending on the docking engine being used.''')

st.subheader("Defining the Binding Location")
st.markdown(''' 
    The binding location is the specific region on the receptor where the ligand is expected to bind. 
    There are two main types of docking: blind docking, where the entire receptor surface is searched for potential binding sites, and site-specific docking, 
    where the binding site is predefined by the user.''')

st.subheader("Docking Ligands to the Receptor")
st.markdown('''
    There are two main steps involved in docking ligands to the receptor: sampling (also known as searching) and scoring. Sampling involves generating multiple possible binding poses for the ligand within the binding site of the receptor.
    Common sampling algorithms are stochastic methods, which include genetic algorithms, Monte Carlo simulations, and molecular dynamics. Other methods such as the incremental adjustment 
    algorithm breaks ligands into smaller fragments and docks them sequentially to the binding site, reassembling them into a complete ligand structure.
    Scoring involves evaluating each binding pose using a scoring function that estimates the binding affinity between the ligand and receptor. The binding pose with the highest score is typically selected as the predicted binding mode. 
    Common scoring functions include empirical scoring functions (which assigns weights to different energy terms such as electrostatic energy, hydrophobicity, and entropy), knowledge-based scoring functions (which take structural information 
    from experimental data and calculate distance-dependent pairwise potentials), and force-field-based scoring functions (which use molecular force fields to evaluate binding interactions). After docking, the best binding poses are saved to structural files (e.g., PDBQT, SDF) for further analysis.''')

st.subheader("Analyzing Docking Results")
st.markdown('''
    After docking is complete, it is important to analyze the results to gain insights into the binding interactions between the ligand and receptor.
    This involves visualizing the predicted binding poses, examining the binding affinities, and identifying key interactions that may contribute to binding.
    Various tools and techniques can be used for this analysis, including molecular visualization software and interaction fingerprinting methods.
    Structures may also need to be converted to a different format depending on the docking engine being used.''')
