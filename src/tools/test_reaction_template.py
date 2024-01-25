from rdkit import Chem

from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def getrxns(rxn, product_smi):
    """
    获取反应规则的引擎对象
    product_smi 反应物
    """
    product_mol = Chem.MolFromSmiles(product_smi)
    reactions = rxn.RunReactants([product_mol])
    rxns = []
    for reaction in reactions:
        smis = []
        for compound in reaction:
            smi = Chem.MolToSmiles(compound)
            smis.append(smi)

        rxnstr = '.'.join(smis) + '>>' + product_smi
        rxns.append(rxnstr)
    return rxns


tem = '([#8:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[N;H0;D3;+0:5](-[C:4])-[C:6])>>Cl-[C;H0;D3;+0:1](-[#8:2])=[O;D1;H0:3].[C:4]-[NH;D2;+0:5]-[C:6]'
rxn = AllChem.ReactionFromSmarts(tem)
product_smi = 'O=C(OCc1ccccc1)N1CC[C@H]2CCCN(CCc3ccccc3)C[C@H]21'
reactions = getrxns(rxn, product_smi)
cnt = 0
for reaction in reactions:
    cnt += 1
    rxn_tmp = AllChem.ReactionFromSmarts(reaction)
    img = Draw.ReactionToImage(
        rxn_tmp
    )
    file_name = 'reaction_template%d.jpg' % cnt
    img.save(
        file_name
    )