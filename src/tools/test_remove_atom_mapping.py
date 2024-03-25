
from rdkit import Chem

def clear_atom_mapping(str_smiles):
    m = Chem.MolFromSmiles(str_smiles)
    for a in m.GetAtoms():
        a.ClearProp('molAtomMapNumber')
    new_smi = Chem.MolToSmiles(m)
    return new_smi

smi = "[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]"
print(clear_atom_mapping(smi))