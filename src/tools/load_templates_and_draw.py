import os
import random
import numpy as np
from tqdm import tqdm, trange

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

def draw_reaction(reaction_smarts, reaction_id):
    '''
    reaction_smarts: str
    reaction_id: int
    '''
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)  
    img = Draw.ReactionToImage(rxn)
    file_name = './reactions/route-4-step-1-reaction-%d.jpg' % (reaction_id)
    img.save(file_name)

template_rule_path = '../dataset/template_rules_1.dat'
template_rules = []

# Halogen_Reacions = []

with open(template_rule_path, 'r') as f:
    for i, l in tqdm(enumerate(f), desc='template rules'):
        rule= l.strip()
        template_rules.append(rule)

#         if not (rule.find('F') == -1 and rule.find('Cl') == -1 and rule.find('Br') == -1 and rule.find('I') == -1):
#             Halogen_Reacions.append(rule)

# print('Halogen_Reacions:', len(Halogen_Reacions))
# template_Halogen_reactions = []

# for reaction in Halogen_Reacions:
#     rule = reaction.replace('F', 'Cl')
#     rule = rule.replace('Br', 'Cl')
#     rule = rule.replace('I', 'Cl')
#     template_Halogen_reactions.append(rule.split('>>')[1])

# print('replace filter len:', len(set(template_Halogen_reactions)))

target_template_id_list = [214796, 359734, 15004, 604, 11656, 115605]
for target_template_id in target_template_id_list:
    print(target_template_id, ':', template_rules[target_template_id])
    draw_reaction(template_rules[target_template_id], target_template_id)


# with open('../dataset/template_library.txt', 'w') as f:
#     for rule in template_rules:
#         f.write('%s\n' % rule)

# route_reaction_list = [
#     "NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](c3nnco3)C[C@H]2CS1>>O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](c3nnco3)C[C@H]2CS1)c1ccccc1",
#     "O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](c3nnco3)C[C@H]2CS1)c1ccccc1>>Cc1ccc(S(=O)(=O)O)cc1.NNC(=O)[C@H]1C[C@H]2CSC(NC(=O)c3ccccc3)=N[C@@]2(c2ccc(F)cc2F)CO1",
#     "NNC(=O)[C@H]1C[C@H]2CSC(NC(=O)c3ccccc3)=N[C@@]2(c2ccc(F)cc2F)CO1>>CC(C)(C)OC(=O)NNC(=O)[C@H]1C[C@H]2CSC(NC(=O)c3ccccc3)=N[C@@]2(c2ccc(F)cc2F)CO1",
#     "CC(C)(C)OC(=O)NNC(=O)[C@H]1C[C@H]2CSC(NC(=O)c3ccccc3)=N[C@@]2(c2ccc(F)cc2F)CO1>>CC(C)(C)OC(=O)NN.COC(CNC(=O)[C@H]1C[C@H]2CSC(NC(=O)c3ccccc3)=N[C@@]2(c2ccc(F)cc2F)CO1)OC",
#     "COC(CNC(=O)[C@H]1C[C@H]2CSC(NC(=O)c3ccccc3)=N[C@@]2(c2ccc(F)cc2F)CO1)OC>>COC(CN)OC.O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](C(=O)O)C[C@H]2CS1)c1ccccc1",
#     "O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](C(=O)O)C[C@H]2CS1)c1ccccc1>>C[N+]1([O-])CCOCC1.O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](CO)C[C@H]2CS1)c1ccccc1",
#     "O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](CO)C[C@H]2CS1)c1ccccc1>>O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](COCc3ccccc3)C[C@H]2CS1)c1ccccc1",
#     "O=C(NC1=N[C@@]2(c3ccc(F)cc3F)CO[C@@H](COCc3ccccc3)C[C@H]2CS1)c1ccccc1>>O=C(NC(=S)N[C@@]1(c2ccc(F)cc2F)CO[C@@H](COCc2ccccc2)C[C@H]1CO)c1ccccc1",
#     "O=C(NC(=S)N[C@@]1(c2ccc(F)cc2F)CO[C@@H](COCc2ccccc2)C[C@H]1CO)c1ccccc1>>N[C@@]1(c2ccc(F)cc2F)CO[C@@H](COCc2ccccc2)C[C@H]1CO.O=C(N=C=S)c1ccccc1",
#     "N[C@@]1(c2ccc(F)cc2F)CO[C@@H](COCc2ccccc2)C[C@H]1CO>>Fc1ccc([C@]23CO[C@@H](COCc4ccccc4)C[C@H]2CON3)c(F)c1",
#     "Fc1ccc([C@]23CO[C@@H](COCc4ccccc4)C[C@H]2CON3)c(F)c1>>Fc1ccc(I)c(F)c1.c1ccc(COC[C@H]2C[C@H]3CON=C3CO2)cc1",
#     "c1ccc(COC[C@H]2C[C@H]3CON=C3CO2)cc1>>C=CC[C@H](COCc1ccccc1)OCC=NO",
#     "C=CC[C@H](COCc1ccccc1)OCC=NO>>NO.C=CC[C@H](COCc1ccccc1)OCC(OCC)OCC",
#     "C=CC[C@H](COCc1ccccc1)OCC(OCC)OCC>>C=CC[C@@H](O)COCc1ccccc1.CCOC(CBr)OCC"
# ]

# reaction_list = [
#     "([#16:3]-[C:2](=[#7:4])-[NH2;D1;+0:1])>>O=C(-[NH;D2;+0:1]-[C:2](-[#16:3])=[#7:4])-c1:c:c:c:c:c:1",
#     "([#8:2]-[C:3]-[c;H0;D3;+0:4]1:[n;H0;D2;+0:6]:[n;H0;D2;+0:7]:[cH;D2;+0:1]:[o;H0;D2;+0:5]:1)>>C-[c;H0;D3;+0:1]1:c:c:c(-S(-O)(=O)=O):c:c:1.[#8:2]-[C:3]-[C;H0;D3;+0:4](=[O;H0;D1;+0:5])-[NH;D2;+0:6]-[NH2;D1;+0:7]",
#     "([C:2]-[NH2;D1;+0:1])>>O=C(-O-C-c1:c:c:c:c:c:1)-[NH;D2;+0:1]-[C:2]",
#     "([C:2]-[NH2;D1;+0:1])>>C-C(-C)(-C)-O-C(=O)-[NH;D2;+0:1]-[C:2]",
#     "([C:1]-[NH2;D1;+0:2])>>[C:1]-[N;H0;D2;+0:2]=[N+]=[N-]"
# ]

# route4_step2_reaction_list = [
#     "([#8:2]-[C:3]-[c;H0;D3;+0:4]1:[n;H0;D2;+0:6]:[n;H0;D2;+0:7]:[cH;D2;+0:1]:[o;H0;D2;+0:5]:1)>>C-[c;H0;D3;+0:1]1:c:c:c(-S(-O)(=O)=O):c:c:1.[#8:2]-[C:3]-[C;H0;D3;+0:4](=[O;H0;D1;+0:5])-[NH;D2;+0:6]-[NH2;D1;+0:7]",
#     "([C:2]-[NH;D2;+0:1]-[C:3])>>C-C(-C)(-C)-O-C(=O)-[N;H0;D3;+0:1](-[C:2])-[C:3]",
#     "([O;D1;H0:9]=[C:8]-[#7:7]-[C;H0;D3;+0:5]1=[N;H0;D2;+0:4]-[C:3]-[C:2]-[CH2;D2;+0:1]-[S;H0;D2;+0:6]-1)>>O-[CH2;D2;+0:1]-[C:2]-[C:3]-[NH;D2;+0:4]-[C;H0;D3;+0:5](=[S;H0;D1;+0:6])-[#7:7]-[C:8]=[O;D1;H0:9]",
#     "([#8:1]-[C@@H;D3;+0:2](-[c:3])-[CH2;D2;+0:4]-[C:5])>>[#8:1]-[C;H0;D3;+0:2](-[c:3])=[CH;D2;+0:4]-[C:5]"
# ]

# route4_step3_reaction_list = [
#     "([#16:3]-[C:2](=[#7:4])-[NH2;D1;+0:1])>>O=C(-[NH;D2;+0:1]-[C:2](-[#16:3])=[#7:4])-c1:c:c:c:c:c:1",
#     "([C:2]-[C;H0;D3;+0:1](=[O;H0;D1;+0:3])-[NH;D2;+0:5]-[N;D1;H2:4])>>O=[C;H0;D3;+0:1](-[C:2])-[OH;D1;+0:3].[N;D1;H2:4]-[NH2;D1;+0:5]",
#     "([#8:6]-[C:5]-[C;H0;D3;+0:3](=[O;H0;D1;+0:4])-[NH;D2;+0:2]-[NH2;D1;+0:1])>>O-n1:[n;H0;D2;+0:1]:[n;H0;D2;+0:2]:c2:c:c:c:c:c:2:1.O=[C;H0;D3;+0:3](-[OH;D1;+0:4])-[C:5]-[#8:6]",
#     "([#8:8]-[C@@H;D3;+0:7](-[C:9])-[C:6]-[C@H;D3;+0:5]1-[C:10]-[#16:11]-[C:2](-[NH2;D1;+0:1])=[#7:3]-[C:4]-1)>>C-C(-C)(-C)-O-C(=O)-[NH;D2;+0:1]-[C:2]1=[#7:3]-[C:4]-[C@@H;D3;+0:5](-[C:6]-[C@@H;D3;+0:7](-[#8:8])-[C:9])-[C:10]-[#16:11]-1",
#     "([C:3]-[C;H0;D3;+0:2](=[O;H0;D1;+0:1])-[NH;D2;+0:5]-[N;D1;H2:4])>>C-[O;H0;D2;+0:1]-[C;H0;D3;+0:2](=O)-[C:3].[N;D1;H2:4]-[NH2;D1;+0:5]",
#     "([C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[NH;D2;+0:5]-[N;D1;H2:4])>>C-C-O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[N;D1;H2:4]-[NH2;D1;+0:5]",
#     "([c:2]:[cH;D2;+0:1]:[c:3])>>F-[c;H0;D3;+0:1](:[c:2]):[c:3]"
# ]

# i=0
# for reaction in route4_step3_reaction_list:
#     draw_reaction(reaction, i)
#     i+=1