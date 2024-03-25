import json
from rdkit import Chem
from time import time
from rdkit.Chem import Draw
from rdchiral import template_extractor

reaction_json = {}
reaction_json['_id'] = 7
reaction_json['reactants'] = 'c1ccc(CCN2CCC[C@@H]3CCN[C@@H]3C2)cc1.O=C(Cl)OCc1ccccc1'
reaction_json['products'] = 'O=C(OCc1ccccc1)N1CC[C@H]2CCCN(CCc3ccccc3)C[C@H]21'

mol1 = Chem.MolFromSmiles('c1ccc(CCN2CCC[C@@H]3CCN[C@@H]3C2)cc1')
mol2 = Chem.MolFromSmiles('O=C(Cl)OCc1ccccc1')
print(template_extractor.extract_from_reaction(reaction_json))