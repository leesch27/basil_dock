"""Additional functions used within basil_dock notebooks."""
import sys, os, glob, re
from ligandsplitter.ligandsplitter.ligandanalysis import group_idxes_from_mol
from ligandsplitter.ligandsplitter.basefunctions import LigandVariables

vars = LigandVariables()

import numpy as np
import pandas as pd
import prolif as plf

from rdkit import Chem
from rdkit.Chem import AllChem,rdFMCS, rdMolAlign
import py3Dmol
from ipywidgets import Layout, Label, Dropdown, Box, HBox
import MDAnalysis as mda 

def get_prot_pockets_data(current_dir, pdb_id, prot_pockets):
    fpocket_out = "data/PDB_files/" + str(pdb_id)+ "_protein_out/"
    f_pocket_dir = os.path.join(current_dir, fpocket_out)
    u_prot = mda.Universe(f"data/PDB_files/{pdb_id}_protein.pdb")
    for file in os.listdir(f_pocket_dir):
        if 'env_atm' in file:
            atoms = []
            res_and_atoms = []
            pocket_num = int(file.split('_')[0].replace('pocket',''))
            out_dir = os.path.join(f_pocket_dir, file)
            with open(out_dir, 'r') as outfile:
                data = outfile.readlines()
            for line in data:
                split_line = line.split()
                if len(split_line) > 1:
                    select_atom_num = split_line[1]
                    select_atom = split_line[2]
                    select_residue_num = split_line[5]
                    # if the residue number for a protein has four digits (greater than 999), split_line[5] will be 
                    # equal to the x coordinate of the atom as the whitespace between the chain identifier and the 
                    # residue number will disappear. the following if statement addresses this
                    if "." in select_residue_num: 
                        temp_residue_num = re.findall(r'\d+', split_line[4])
                        temp2_residue_num = ''.join(str(x) for x in temp_residue_num)
                        select_residue_num = int(temp2_residue_num)
                    atoms.append(select_atom_num)
                    md_input1 = "(resid " + str(select_residue_num) + " and name " + str(select_atom) + ")"
                    res_and_atoms.append(md_input1)

            # get center of docking box
            atom_string = ', '.join(str(x) for x in atoms)
            res_and_atom_string = ' or '.join(str(x) for x in res_and_atoms)
            pocket_mda = u_prot.select_atoms(res_and_atom_string)
            pocket_center = pocket_mda.center_of_geometry()
            pocket_center_list = np.ndarray.tolist(pocket_center)

            # get size of docking box
            ligand_box = pocket_mda.positions.max(axis=0) - pocket_mda.positions.min(axis=0)
            ligand_box_list = np.ndarray.tolist(ligand_box)
            ligand_box_list2 = []
            for value in ligand_box_list:
                if value < 0:
                    ligand_box_list2.append(float(value - 5))
                elif value > 0:
                    ligand_box_list2.append(float(value + 5))
                else:
                    ligand_box_list2.append(float(0))
        
            prot_pockets.loc[pocket_num,'center_x'] = pocket_center_list[0]
            prot_pockets.loc[pocket_num,'center_y'] = pocket_center_list[1]
            prot_pockets.loc[pocket_num,'center_z'] = pocket_center_list[2]
            prot_pockets.loc[pocket_num,'size_x'] = abs(ligand_box_list2[0])
            prot_pockets.loc[pocket_num,'size_y'] = abs(ligand_box_list2[1])
            prot_pockets.loc[pocket_num,'size_z'] = abs(ligand_box_list2[2])

def fetch_data_files():
    receptor_files = []
    prot_pocket_csvs = []
    ligand_info_csvs = []

    file_location_data = os.path.join('data', 'PDB_files', '*')
    receptors = glob.glob(file_location_data)

    for file in receptors:
        if ".pdb" in file:
            file_short = file.split("/")[-1]
            pdb_id_initial = file_short.split("_")[:-1]
            # do not add "clean" or "protein" to the list of receptor files
            if "clean" not in pdb_id_initial and "protein" not in pdb_id_initial and "ligand" not in pdb_id_initial:
                if len(pdb_id_initial) > 1:
                    pdb_id = '_'.join(str(x) for x in pdb_id_initial)
                elif len(pdb_id_initial) == 1:
                    pdb_id = pdb_id_initial[0]
                if pdb_id not in receptor_files:
                    receptor_files.append(pdb_id)
        elif ".cif" in file:
            file_short = file.split("/")[-1]
            pdb_id = file_short.split(".")[0]
            if pdb_id not in receptor_files:
                receptor_files.append(pdb_id)
        elif ".ent" in file:
            file_short = file.split("/")[-1]
            pdb_id_long = file_short.split(".")[0]
            if "pdb" in pdb_id_long[:3]:
                pdb_id = pdb_id_long[3:]
                if pdb_id not in receptor_files:
                    receptor_files.append(pdb_id)
            else:
                if pdb_id_long not in receptor_files:
                    receptor_files.append(pdb_id_long)
    
    file_location_data = os.path.join('data','*.csv')
    csvs = glob.glob(file_location_data)

    for file in csvs:
        if "protein_pockets" in file:
            prot_pocket_csvs.append(file)
        elif "ligand_information" in file:
            ligand_info_csvs.append(file)
    
    return receptor_files, prot_pocket_csvs, ligand_info_csvs

def fetch_docking_data():
    docking_files = []

    file_location_data = os.path.join('data','*.csv')
    csvs = glob.glob(file_location_data)

    for file in csvs:
        if "docking_information" in file:
            docking_files.append(file)
        elif("protein_pockets" not in file) and ("ligand_information" not in file) and ("smiles_data" not in file):
            if "test_train_bioactive_data" in file:
                continue
            data = pd.read_csv(file)
            columns = data.columns
            if("Frame" in columns) and ("Score" in columns) and ("Ligand" in columns):
                docking_files.append(file)
    
    return docking_files

def get_ifps(dock_value, dock_engine, select_ligs, protein_plf, prot_pockets = pd.DataFrame()):
    all_ligand_plf = []
    ligand_plf_descriptors = []
    all_df = []
    all_ifps = []
    for i in select_ligs:
        if dock_value == "Blind docking":
            for pocket in prot_pockets.index:
                lig_suppl = plf.sdf_supplier(f"data/{dock_engine}_out_2/{i}_pocket_{pocket}_{dock_engine}_out_2.sdf")
                fp = plf.Fingerprint(count=True)
                fp.run_from_iterable(lig_suppl, protein_plf)
                results_df = fp.to_dataframe()
                all_df.append(results_df)
                for ind, lig in enumerate(lig_suppl):
                    all_ligand_plf.append(lig)
                    ligand_plf_descriptors.append(f"Ligand {i}, Pocket {pocket}, Pose {ind + 1}")
                    ifp = fp.generate(lig, protein_plf, metadata = True)
                    all_ifps.append(ifp)
        else:
            lig_suppl = plf.sdf_supplier(f"data/{dock_engine}_out_2/{i}_{dock_engine}_out_2.sdf")
            fp = plf.Fingerprint(count=True)
            fp.run_from_iterable(lig_suppl, protein_plf)
            results_df = fp.to_dataframe()
            all_df.append(results_df)
            for ind, lig in enumerate(lig_suppl):
                all_ligand_plf.append(lig)
                ligand_plf_descriptors.append(f"Ligand {i}, Pose {ind + 1}")
                ifp = fp.generate(lig, protein_plf, metadata = True)
                all_ifps.append(ifp)
    return all_df, all_ifps, all_ligand_plf, ligand_plf_descriptors

def get_scores(dock_value, dock_engine, select_ligs, prot_pockets = pd.DataFrame()):
    all_results = []
    scores = [] # get list of scores for each pose
    for h, i in enumerate(select_ligs):
        # initialize list that contains values all poses in all pockets (for 1 ligand at a time)
        nested_results = []
        # append pose data for each pocket to nested_results
        if dock_value == "Blind docking":
            for pocket in prot_pockets.index:
                results = Chem.SDMolSupplier(f"data/{dock_engine}_out_2/{i}_pocket_{pocket}_{dock_engine}_out_2.sdf")
                nested_results.append(results)
        else:
            results = Chem.SDMolSupplier(f"data/{dock_engine}_out_2/{i}_{dock_engine}_out_2.sdf")
            nested_results.append(results) 
        # add all values in nested_results to allResults list
        all_results.append(nested_results)
    
    # get score values for every pose in allResults
    if dock_value == "Blind docking":
        for linenum, i in enumerate(all_results):
            for num, pocket in enumerate(i):
                for num2, pose in enumerate(pocket):
                    if dock_engine == "smina":
                        scores.append(float(all_results[linenum][num][num2].GetProp('minimizedAffinity')))
                    else:
                        scores.append(float(all_results[linenum][num][num2].GetProp('Score')))
    else:
        for linenum, i in enumerate(all_results):
            for num, pose in enumerate(i):
                for num2, item in enumerate(pose):
                    if dock_engine == "smina":
                        scores.append(float(all_results[linenum][num][num2].GetProp('minimizedAffinity')))
                    else:
                        scores.append(float(all_results[linenum][num][num2].GetProp('Score')))
    return scores

def get_largest_array_column(df, docking_type):
    largest_array_column = {}
    if docking_type == "Blind docking":
        start_col_num = 2
    else:
        start_col_num = 1
    for col_num, column in enumerate(df):
        largest_array = 0
        if col_num > start_col_num:
            for row in df[column]:
                if int(row) > largest_array:
                    largest_array = int(row)
            largest_array_column[column] = largest_array
    return largest_array_column

def expand_df(all_ifps, df, df2, largest_array_column):
    col_names_list = []
    residues = []
    interactions = []
    total_counter = 0
    for key in all_ifps:
        for key_new in key:
            for key_2 in key[key_new]:
                residues.append(str(key_new[1]))
                interactions.append(str(key_2))
                lig_name = str(key_new[0])
                res_name = str(key_new[1])
                column_name = (lig_name, res_name, key_2)
                if column_name not in col_names_list:
                    df2[column_name] = df[column_name]
                    number_of_ints = int(largest_array_column[column_name])
                    counter_ind = 0
                    while counter_ind < number_of_ints:
                        df2[(lig_name, res_name, f"Functional group involved ({key_2}){counter_ind}")] = pd.Series([0] * df.shape[0])
                        #df2[(lig_name, res_name, f"Residue type({key_2}){counter_ind}")] = pd.Series([0] * df.shape[0])
                        df2[(lig_name, res_name, f"Distance ({key_2}){counter_ind}")] = pd.Series([0] * df.shape[0])
                        df2[(lig_name, res_name, f"Index 1 (Ligand) ({key_2}){counter_ind}")] = pd.Series([0] * df.shape[0])
                        df2[(lig_name, res_name, f"Index 2 (Ligand) ({key_2}){counter_ind}")] = pd.Series([0] * df.shape[0])
                        df2[(lig_name, res_name, f"Index 3 (Protein) ({key_2}){counter_ind}")] = pd.Series([0] * df.shape[0])
                        df2[(lig_name, res_name, f"Index 4 (Protein) ({key_2}){counter_ind}")] = pd.Series([0] * df.shape[0])
                        counter_ind += 1
                total_counter += 1

def fill_df(df2, all_ifps, all_ligand_plf, largest_array_column):
    total_counter = 0
    for number, key in enumerate(all_ifps):
        for key_new in key:
            for key_2 in key[key_new]:
                lig_name = str(key_new[0])
                res_name = str(key_new[1])
                column_name = (lig_name, res_name, key_2)
                get_pose = df2["Frame"][number]
                x = key[key_new]
                y = x[key_2]
                df_groups = [0] * largest_array_column[column_name]
                #df_residue = [0] * largest_array_column[column_name]
                df_distance = [0] * largest_array_column[column_name]
                df_ind_1 = [0] * largest_array_column[column_name]
                df_ind_2 = [0] * largest_array_column[column_name]
                df_ind_3 = [0] * largest_array_column[column_name]
                df_ind_4 = [0] * largest_array_column[column_name]
                for inst_num, instance in enumerate(y):
                    distance = instance["distance"]
                    df_distance[inst_num] = distance
                    found_res = res_name[:3]
                    #df_residue[inst_num] = (vars.type_dict[found_res])
                    parent_index = instance["parent_indices"]
                    if len(parent_index["ligand"]) == 2:
                        df_ind_1[inst_num] = parent_index["ligand"][0]
                        df_ind_2[inst_num] = parent_index["ligand"][1]
                    else:
                        df_ind_1[inst_num] = parent_index["ligand"][0]
                        df_ind_2[inst_num] = 0
                    if len(parent_index["protein"]) == 2:
                        df_ind_3[inst_num] = parent_index["protein"][0]
                        df_ind_4[inst_num] = parent_index["protein"][1]
                    else:
                        df_ind_3[inst_num] = parent_index["protein"][0]
                        df_ind_4[inst_num] = 0
                    current = all_ligand_plf[number]
                    group_ints = group_idxes_from_mol(current, renumber = False)
                    for value in group_ints.keys():
                        if len(parent_index["ligand"]) == 2:
                            if value == str(parent_index["ligand"][0]) or value == str(parent_index["ligand"][1]):
                                df_groups[inst_num] = int(vars.groups_to_numbers[group_ints[value][0]])
                        else:
                            if value == str(parent_index["ligand"][0]):
                                df_groups[inst_num] = int(vars.groups_to_numbers[group_ints[value][0]])
                number_of_ints = int(largest_array_column[column_name])
                counter_ind = 0
                while counter_ind < number_of_ints:
                    df2.at[number, (lig_name, res_name, f"Functional group involved ({key_2}){counter_ind}")] = df_groups[counter_ind]
                    #df2.at[number, (lig_name, res_name, f"Residue type({key_2}){counter_ind}")] = df_residue[counter_ind]
                    df2.at[number, (lig_name, res_name, f"Distance ({key_2}){counter_ind}")] = df_distance[counter_ind]
                    df2.at[number, (lig_name, res_name, f"Index 1 (Ligand) ({key_2}){counter_ind}")] = df_ind_1[counter_ind]
                    df2.at[number, (lig_name, res_name, f"Index 2 (Ligand) ({key_2}){counter_ind}")] = df_ind_2[counter_ind]
                    df2.at[number, (lig_name, res_name, f"Index 3 (Protein) ({key_2}){counter_ind}")] = df_ind_3[counter_ind]
                    df2.at[number, (lig_name, res_name, f"Index 4 (Protein) ({key_2}){counter_ind}")] = df_ind_4[counter_ind]
                    counter_ind += 1
                total_counter += 1

def save_dataframe(df2, dock_engine, pdb_id, lig_value, csv_name = ""):
    col_names_prelim = df2.keys()
    col_names = {}
    for column in col_names_prelim:
        index1 = column[0]
        index2 = column[1]
        index3 = column[2]
        if "Frame" in index1:
            col_names[str(column)] = "Frame"
        elif "Score" in index1:
            col_names[str(column)] = "Score"
        elif "Ligand" in index1:
            col_names[str(column)] = "Ligand"
        elif "Pocket" in index1:
            col_names[str(column)] = "Pocket"
        else:
            new_name = f"{index2} {index3}"
            col_names[str(column)] = new_name

    if csv_name == "":
        csv_title = f'data/{pdb_id}_{str(len(lig_value))}_ligands_docking_information_{dock_engine}.csv'
    else:
        csv_title = f'data/{csv_name}.csv'
    df2 = df2.convert_dtypes()
    df2.to_csv(csv_title, index = False, header = col_names_prelim)
    df3 = pd.read_csv(csv_title)
    df3 = df3.rename(columns = col_names)
    df3.to_csv(csv_title, index = False, header = col_names)

def compare_poses_form(select_type, select_ligs, prot_pockets, pose_mode):
    all_pock_nums = {}
    all_pose_nums = {}
    form_item_layout = Layout(
        display='flex',
        flex_flow='row',
        justify_content='space-between')

    if select_type == "Blind docking":
        pocket_list = []
        for pocket in prot_pockets.index:
            pocket_list.append((str(pocket), int(pocket)))

        ligand_number = Dropdown(options = select_ligs)

        pock_number1 = Dropdown(options = pocket_list)
        pose_number1 = Dropdown(options = [('1', 1), ('2', 2), ('3', 3), ('4', 4), ('5', 5)])
        if pose_mode == 2:
            pock_number2 = Dropdown(options = pocket_list)
            pose_number2 = Dropdown(options = [('1', 1), ('2', 2), ('3', 3), ('4', 4), ('5', 5)])

        form_items1 = [Box([Label(value='Ligand'), ligand_number], layout=form_item_layout)]

        if pose_mode == 2:
            form_items2 = [Box([Label(value='Pocket Number (View 1)'), pock_number1], layout=form_item_layout),
                        Box([Label(value='Pose Number (View 1)'), pose_number1], layout=form_item_layout),
                        Box([Label(value='Pocket Number (View 2)'), pock_number2], layout=form_item_layout),
                        Box([Label(value='Pose Number (View 2)'), pose_number2], layout=form_item_layout)]
            all_pock_nums = {1: pock_number1, 2: pock_number2}
            all_pose_nums = {1: pose_number1, 2: pose_number2}
        else:
            form_items2 = [Box([Label(value='Pocket Number (View 1)'), pock_number1], layout=form_item_layout),
                        Box([Label(value='Pose Number (View 1)'), pose_number1], layout=form_item_layout)]
            all_pock_nums = {1: pock_number1}
            all_pose_nums = {1: pose_number1}

        form1 = Box(form_items1, layout=Layout(
            display='flex',
            flex_flow='column',
            border='solid 2px',
            align_items='stretch',
            width='50%'
        ))
        form2 = Box(form_items2, layout=Layout(
            display='flex',
            flex_flow='column',
            border='solid 2px',
            align_items='stretch',
            width='50%'
        ))

        form = HBox([form1, form2])
    else:
        ligand_number = Dropdown(options = select_ligs)
        pose_number1 = Dropdown(options = [('1', 1), ('2', 2), ('3', 3), ('4', 4), ('5', 5)])
        if pose_mode == 2:
            pose_number2 = Dropdown(options = [('1', 1), ('2', 2), ('3', 3), ('4', 4), ('5', 5)])
    
        form_items1 = [Box([Label(value='Ligand'), ligand_number], layout=form_item_layout)]
    
        if pose_mode == 2:
            form_items2 = [Box([Label(value='Pose Number (View 1)'), pose_number1], layout=form_item_layout),
                        Box([Label(value='Pose Number (View 2)'), pose_number2], layout=form_item_layout)]
            all_pose_nums = {1: pose_number1, 2: pose_number2}
        else:
            form_items2 = [Box([Label(value='Pose Number (View 1)'), pose_number1], layout=form_item_layout)]
            all_pose_nums = {1: pose_number1}
    
        form1 = Box(form_items1, layout=Layout(
            display='flex',
            flex_flow='column',
            border='solid 2px',
            align_items='stretch',
            width='50%'
        ))
        form2 = Box(form_items2, layout=Layout(
            display='flex',
            flex_flow='column',
            border='solid 2px',
            align_items='stretch',
            width='50%'
        ))
    
        form = HBox([form1, form2])
    return form, ligand_number, all_pock_nums, all_pose_nums

def compare_poses(pdb_id, select_type, pose_mode, ligand_number, select_dock, form_items1, form_items2):
    if select_type == "Blind docking":
        if pose_mode == 1:
            pocket_selection_list = [int(form_items1[1])]
            pose_selection_list = [int(form_items2[1])]
        else:
            pocket_selection_list = [int(form_items1[1]), int(form_items1[2])]
            pose_selection_list = [int(form_items2[1]), int(form_items2[2])]

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
        if pose_mode == 1:
            view.addModel(open(f"data/MOL2_files/{ligand_number}.mol2",'r').read(),format='mol2')
            ref_m = view.getModel()
            ref_m.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
        else:
            selected_ref = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_pocket_{pocket_selection_list[1]}_{select_dock}_out_2.sdf")
            p_ref=Chem.MolToMolBlock(selected_ref[pose_selection_list[1] - 1],False)
            view.addModel(p_ref,'reference mol')
            x_ref = view.getModel()
            x_ref.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
            if select_dock == 'smina':
                print('Score: {}'.format(selected_ref[pose_selection_list[1] - 1].GetProp('minimizedAffinity')))
            else:
                print('Pose: {} | Score: {}'.format(selected_ref[pose_selection_list[1] - 1].GetProp('Pose'), selected_ref[pose_selection_list[1] - 1].GetProp('Score')))
    
        # add experimental docking data of a desired pose in a pocket to py3dmol viewer
        selected = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_pocket_{pocket_selection_list[0]}_{select_dock}_out_2.sdf")
        p=Chem.MolToMolBlock(selected[pose_selection_list[0] - 1],False)
        if select_dock == 'smina':
            print('Smina Pose (' + str(ligand_number) + '): Cyan')
            print('Score: {}'.format(selected[pose_selection_list[0] - 1].GetProp('minimizedAffinity')))
        else:
            print('Vina Pose (' + str(ligand_number) + '): Cyan')
            print('Pose: {} | Score: {}'.format(selected[pose_selection_list[0] - 1].GetProp('Pose'), selected[pose_selection_list[0] - 1].GetProp('Score')))
        view.addModel(p,'mol')
        x = view.getModel()
        x.setStyle({},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})
    else:
        if pose_mode == 1:
            pose_selection_list = [int(form_items2[1])]
        else:
            pose_selection_list = [int(form_items2[1]), int(form_items2[2])]
    
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
        if pose_mode == 1:
            view.addModel(open(f"data/MOL2_files/{ligand_number}.mol2",'r').read(),format='mol2')
            ref_m = view.getModel()
            ref_m.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
        else:
            selected_ref = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_{select_dock}_out_2.sdf")
            p_ref=Chem.MolToMolBlock(selected_ref[pose_selection_list[1] - 1],False)
            view.addModel(p_ref,'reference mol')
            x_ref = view.getModel()
            x_ref.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            print('Reference (' + str(ligand_number) + '): Magenta ')
            if select_dock == 'smina':
                print('Score: {}'.format(selected_ref[pose_selection_list[1] - 1].GetProp('minimizedAffinity')))
            else:
                print('Pose: {} | Score: {}'.format(selected_ref[pose_selection_list[1] - 1].GetProp('Pose'), selected_ref[pose_selection_list[1] - 1].GetProp('Score')))
    
        # add experimental docking data of a desired pose in a pocket to py3dmol viewer
        selected = Chem.SDMolSupplier(f"data/{select_dock}_out_2/{ligand_number}_{select_dock}_out_2.sdf")
        p=Chem.MolToMolBlock(selected[pose_selection_list[0] - 1],False)
        if select_dock == 'smina':
            print('Smina Pose (' + str(ligand_number) + '): Cyan')
            print ('Score: {}'.format(selected[pose_selection_list[0] - 1].GetProp('minimizedAffinity')))
        else:
            print('Vina Pose (' + str(ligand_number) + '): Cyan')
            print ('Pose: {} | Score: {}'.format(selected[pose_selection_list[0] - 1].GetProp('Pose'), selected[pose_selection_list[0] - 1].GetProp('Score')))
        view.addModel(p,'mol')
        x = view.getModel()
        x.setStyle({},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})
        
    view.zoomTo()
    view.show()