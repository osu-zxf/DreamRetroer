from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def draw_molecule_to_png(smi, file_name):
    '''
    use RDKit to draw molecule & save as [file_name](should be png)
    input:  1. smi: str, SMILES
            2. file_name : str, target file path
    '''
    mol = Chem.MolFromSmiles(smi)
    Draw.MolToFile(mol, file_name, size=(350, 300))

def draw_reaction_to_png(reaction_smarts, file_name):
    '''
    use RDKit to draw reaction SMARTS & save as [file_name](should be png)
    input:  1. reaction_smarts: str
            2. file_name: str, target file path
    '''
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)  
    img = Draw.ReactionToImage(rxn)
    img.save(file_name)

if __name__ == "__main__":
    # smis=['[C:2]-[NH;D2;+0:1]-[C:3]']
    # cnt = 0
    # for smi in smis:
    #     cnt += 1
    #     file_name = 'mol%d.png' % cnt
    #     draw_molecule_to_png(smi, file_name)

    # smarts = '[F,Cl,Br,I][CH2;D2;+0:1][C:2]'
    # smarts = "O=[C;H0;D3;+0:1](-[C:2])-[C;D1;H3:3]"
    # mol = Chem.MolFromSmarts(smarts)
    # Draw.MolToFile(mol, 'tmp.png', size=(350, 300))

    # react = '([C:3]-[N;H0;D3;+0:4](-[C:5])-[CH2;D2;+0:1]-[C:2])>>[F,Cl,Br,I]-[CH2;D2;+0:1]-[C:2].[C:3]-[NH;D2;+0:4]-[C:5]'

    # reacts = ['([C:3]-[N;H0;D3;+0:4](-[C:5])-[CH2;D2;+0:1]-[C:2])>>Br-[CH2;D2;+0:1]-[C:2].[C:3]-[NH;D2;+0:4]-[C:5]',
    #             '([C:3]-[N;H0;D3;+0:4](-[C:5])-[CH2;D2;+0:1]-[C:2])>>O=[CH;D2;+0:1]-[C:2].[C:3]-[NH;D2;+0:4]-[C:5]',
    #             '([C:2]-[CH2;D2;+0:1]-[N;H0;D3;+0:4](-[C:3])-[C:5])>>C-S(=O)(=O)-O-[CH2;D2;+0:1]-[C:2].[C:3]-[NH;D2;+0:4]-[C:5]',
    #             '([C:2]-[CH2;D2;+0:1]-[N;H0;D3;+0:4](-[C:3])-[C:5])>>I-[CH2;D2;+0:1]-[C:2].[C:3]-[NH;D2;+0:4]-[C:5]',
    #             '([C:3]-[N;H0;D3;+0:4](-[C:5])-[CH2;D2;+0:1]-[C:2])>>Cl-[CH2;D2;+0:1]-[C:2].[C:3]-[NH;D2;+0:4]-[C:5]',
    #             '([C:2]-[CH2;D2;+0:1]-[N;H0;D3;+0:4](-[C:3])-[C:5])>>C-c1:c:c:c(-S(=O)(=O)-O-[CH2;D2;+0:1]-[C:2]):c:c:1.[C:3]-[NH;D2;+0:4]-[C:5]',
    #             '([C:2]-[CH2;D2;+0:1]-[N;H0;D3;+0:4](-[C:3])-[C:5])>>O-[CH2;D2;+0:1]-[C:2].[C:3]-[NH;D2;+0:4]-[C:5]']
    
    reacts = ["[Br:1][CH2:2][CH2:3][OH:4].[CH2:5]([S:7](Cl)(=[O:9])=[O:8])[CH3:6].CCOCC>>[CH2:5]([S:7]([O:4][CH2:3][CH2:2][Br:1])(=[O:9])=[O:8])[CH3:6]",
              "[C:5]-[O;H0;D2;+0:6]-[S;H0;D4;+0:1](-[C:2])(=[O;D1;H0:3])=[O;D1;H0:4]>>Cl-[S;H0;D4;+0:1](-[C:2])(=[O;D1;H0:3])=[O;D1;H0:4].[C:5]-[OH;D1;+0:6]"]
    cnt = 0
    reacts = ["CCOC(=O)c1cccc(N2CC[C@H](NC(=O)c3nc(C(F)(F)F)c(CC)[nH]3)[C@H](OC)C2)c1>>CCOC(=O)c1cccc(N2CC[C@H](N)[C@H](OC)C2)c1.CCc1[nH]c(C(=O)O)nc1C(F)(F)F",
              "CCc1[nH]c(C(=O)N[C@H]2CCN(C(=O)OC(C)(C)C)C[C@H]2OC)nc1C(F)(F)F>>CO[C@@H]1CN(C(=O)OC(C)(C)C)CC[C@@H]1N.CCc1[nH]c(C(=O)O)nc1C(F)(F)F"]
    for re in reacts:
        cnt += 1
        draw_reaction_to_png(re, 'tmp%d.png' % cnt)