from rdkit import Chem
from rdkit.Chem import AllChem, Draw, DataStructs, RDConfig, rdBase

def get_fxnal_groups(ligand):
    mol_groups = []
    mol = Chem.MolFromMol2File(f"data/MOL2_files/{ligand}_H.mol2",sanitize=False)
    #common functional groups
    ester = '[CX3H1,CX3](=O)[OX2H0][C;!$([C]=[O])]'
    ether = '[OD2]([C;!$([C]=[O])])[C;!$([C]=[O])]'
    hydroxy = '[OX2H;$([O]([C])[H]);!$([O](C=O)[H])][H]'
    carbox_acid = '[CX3;!$([CX3][H])](=O)[OX2H][H]'
    aldehyde = '[CX3;$(C([H])(=O)[C])](=[O;!$([O][O])])[H]'
    anhydr = 'C-[CX3](=O)[O][CX3](=O)-C'

    #nitrogen containing
    amine = '[C][NX3;H2;!$(NC=O)]([H])[H]'
    amine_2 = '[C][NX3;H;!$(NC=O)]([C])[H]'
    amine_3 = '[C][NX3;H0;!$(NC=O);!$(N=O)]([CX4])[CX4]'
    amide = '[CX3](=O)[NX3;!$(N=O)]([H])[H]' 
    amide_2 = '[CX3](=O)[NX3;!$(N=O)]([C])[H]'
    amide_3 = '[CX3](=O)[NX3;!$(N=O)]([C])[C]'
    nitro = 'C-N(=O)-O'
    imine = '[C;$(C([C,H])([C,H])=[N])]=[N][C,H]'

    #halogens
    f_hal = 'F'
    cl_hal = 'Cl'
    br_hal = 'Br'
    i_hal = 'I'

    #multiply bonded and rings
    alkene = 'C=C'
    alkyne = 'C#C-[H]'
    phenyl = '[CX4;$([C][c;$(c1cc[c]cc1)]);!$([C]([H])([H])([C])[c;$(c1cc[c]cc1)])][c]1:[c]:[c]:[c]:[c]:[c]1'
    benzyl = '[CX4;$([C](c1ccccc1)([H])([H])[C])][c]1:[c]:[c]:[c]:[c]:[c]1'
    pyrrole = '[$([NH]1C=CC=C1)]'
    imidiz = '[$([NH]1C=NC=C1)]'

    functional_groups = [ester, ether, hydroxy, carbox_acid, aldehyde, anhydr, amine, amine_2, amine_3, amide, amide_2, amide_3, nitro, imine, f_hal, cl_hal, br_hal, i_hal, alkene, alkyne, phenyl, benzyl, pyrrole, imidiz]

    functional_groups_dict = {
        '[CX3H1,CX3](=O)[OX2H0][C;!$([C]=[O])]': 'ester',
        '[OD2]([C;!$([C]=[O])])[C;!$([C]=[O])]': 'ether',
        '[OX2H;$([O]([C])[H]);!$([O](C=O)[H])][H]': 'hydroxy',
        '[CX3;!$([CX3][H])](=O)[OX2H][H]': 'carbox_acid',
    	'[CX3;$(C([H])(=O)[C])](=[O;!$([O][O])])[H]':'aldehyde',
    	'C-[CX3](=O)[O][CX3](=O)-C)': 'anhydr',
    	'[C][NX3;H2;!$(NC=O)]([H])[H]': 'amine',
        '[C][NX3;H;!$(NC=O)]([C])[H]': 'amine_2',
    	'[C][NX3;H0;!$(NC=O);!$(N=O)]([CX4])[CX4]': 'amine_3',
    	'[CX3](=O)[NX3;!$(N=O)]([H])[H]': 'amide',
    	'[CX3](=O)[NX3;!$(N=O)]([C])[H]': 'amide_2',
    	'[CX3](=O)[NX3;!$(N=O)]([C])[C]': 'amide_3',
    	'C-N(=O)-O': 'nitro',
    	'[C;$(C([C,H])([C,H])=[N])]=[N][C,H]': 'imine',
    	'F':'f_hal',
    	'Cl': 'cl_hal',
    	'Br': 'br_hal',
    	'I': 'i_hal',
    	'C=C': 'alkene',
    	'C#C-[H]': 'alkyne',
    	'[CX4;$([C][c;$(c1cc[c]cc1)]);!$([C]([H])([H])([C])[c;$(c1cc[c]cc1)])][c]1:[c]:[c]:[c]:[c]:[c]1': 'phenyl',
    	'[CX4;$([C]([H])([H])([C])[c;$(c1cc[c]cc1)])][c]1:[c]:[c]:[c]:[c]:[c]1': 'benzyl',
    	'[$([NH]1C=CC=C1)]': 'pyrrole',
        '[$([NH]1C=NC=C1)]': 'imidiz'
    }
    for j in functional_groups:
        k = Chem.MolFromSmarts(j)
        if mol.HasSubstructMatch(k):
            mol_groups.append(1) # 1 corresponds to a functional group being present
            print(f'Atom indices for substructure match (type: {functional_groups_dict[j]}) in ligand {i}: ')
            print(mol.GetSubstructMatches(k))
        else:
            mol_groups.append(0) # 0 corresponds to a functional group being absent
    return mol_groups
