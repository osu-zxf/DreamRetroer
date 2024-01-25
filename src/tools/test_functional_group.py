from rdkit import Chem

smis = ['BrCCc1ccccc1', 'ICCc1ccccc1', 'ClCCc1ccccc1']

smart = Chem.MolFromSmarts('[F,Cl,Br,I][CH2;D2;+0:1][C:2]')

for smi in smis:
    mol = Chem.MolFromSmiles(smi)
    print(mol.GetSubstructMatches(smart))