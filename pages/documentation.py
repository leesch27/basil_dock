import streamlit as st

st.title("basil_dock Documentation")

st.header("Overview", divider="gray")

st.header("Pages", divider="gray")

st.subheader("General Docking Parameters")
st.write("This page allows users to set various parameters that will be used throughout the docking process, including selecting docking software, selecting docking method (blind or site-specific), and setting variables.")
st.write("*Select Docking Engine to Use*: Select docking software to be used (AutoDock Vina or Smina)")
st.write("*Select Docking Method to Use*: Select whether to perform blind docking (binding site not predefined) or site-specific docking (binding site predefined by user)")
st.write("*Select Exhaustiveness Value*: Set exhaustiveness parameter for docking runs (higher values increase accuracy but also increase computation time)")
st.write("*Set Number of Docking Modes to Generate*: Set number of binding poses to generate for each ligand during docking")
st.write("*Upload Protein Receptor To Use*: Upload the receptor structure file to use in docking in CIF, ENT, or PDB format")
st.write("*Upload Ligand/s To Use*: Upload one or more ligand structure files to use in docking in MOL2 or SDF format")

st.subheader("Blind Docking")
st.write("This page provides options for performing blind docking studies, where the binding site is not predefined.")


st.subheader("Site-Specific Docking")
st.write("This page provides options for performing site-specific docking studies, where the binding site is predefined by the user.")


st.subheader("Docking Analysis")
st.write("This page allows users to analyze the results of their docking studies, including viewing binding poses and calculating interaction fingerprints.")


st.subheader("Get Ligands via RCSB PDB Advanced Search")
st.write("This page allows users to search for ligands in the RCSB PDB database using advanced search options.")
st.write("*Search by Chemical Name*: Enter the chemical name to find relevant ligand structures.")
st.write("*Search by Chemical Name Synonym*: Enter the chemical name synonym to find relevant ligand structures.")
st.write("*Search by Chemical ID*: Enter the RCSB Chemical ID to find relevant ligand structures.")
st.write("*Search by Chemical Type*: Select the chemical type to find relevant ligand structures.")
st.write("*Search by Chemical Brand Name*: Enter the DrugBank brand name to find relevant ligand structures.")
st.write("*Search by Formula Similarity*: Enter the ligand formula to find relevant ligand structures.")
st.write("*Search by Structure Similarity*: Enter the ligand SMILES string to find relevant ligand structures.")

st.subheader("Generate Ligands from SMILES Strings")
st.write("This page provides options for generating 3D ligand structures from SMILES strings.")


st.subheader("Derive Ligands from Existing Ligands")
st.write("This page provides options for deriving new ligands from existing ligand structures.")


st.subheader("Get Receptors via RCSB PDB Advanced Search")
st.write("This page allows users to search for receptors in the RCSB PDB database using advanced search options.")
st.write("*Search for Proteins by Enzyme Classification Name*: Enter the enzyme classification name to find relevant protein structures.")
st.write("*Search for Proteins by Enzyme Classification Number*: Enter the enzyme classification number to find relevant protein structures.")
st.write("*Search for Proteins by Number of Chains*: Enter the number of chains to find relevant protein structures.")
st.write("*Search for Proteins by Length of Sequence*: Enter the length of the sequence to find relevant protein structures.")
st.write("*Search for Proteins by Molecular Weight*: Enter the molecular weight to find relevant protein structures.")