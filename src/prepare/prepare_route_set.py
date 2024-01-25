import pandas as pd
import json
from mol_tree import MolTree
from rdkit import Chem

def load_json_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data

def prepare_starting_molecules(filename):
    '''
    load starting molecules from csv file
    return: list of mol
    '''
    print('Loading starting molecules from %s' % filename)

    if filename[-3:] == 'csv':
        mols_list = list(pd.read_csv(filename)['mol'])
        starting_mols = set(mols_list)

    print('%d starting molecules loaded' % len(starting_mols))
    print('start mol example:', list(starting_mols)[0])
    return starting_mols

def clear_atom_mapping(str_smiles):
    m = Chem.MolFromSmiles(str_smiles)
    for a in m.GetAtoms():
        a.ClearProp('molAtomMapNumber')
    new_smi = Chem.MolToSmiles(m)
    return new_smi

def extract_reactions_and_templates(json_data_reaction, json_data_template):
    '''
    return a dict:
    key: product smiles
    value: a list [(reaction_id, a list of reactants[]), ...]
    '''
    reaction_dict = {}
    template_dict = {}
    for item in json_data_reaction:
        if "_id" in item and "reactants" in item and "products" in item and "spectators" in item:
            reaction = ([clear_atom_mapping(x) for x in item["reactants"].split('.')], item["products"])
            reaction_dict[item["_id"]] = reaction   # ([reactants, ...], product)
    print('process reaction file finished')
    
    for item in json_data_template:
        if "reaction_id" in item and "reaction_smarts" in item and "products" in item and "reactants" in item:
            # just solve only one product
            if len(item["products"].split('.')) >= 2:
                continue
            retro_template = item["reaction_smarts"]
            template_dict[item["reaction_id"]] = retro_template
    print('process template file finished')

    # all routes must have template
    all_single_steps_set = {}
    for reaction_id, rxn_template in template_dict.items():
        if reaction_dict.__contains__(reaction_id):
            product = reaction_dict[reaction_id][1]
            if not all_single_steps_set.__contains__(product):
                all_single_steps_set[product] = []
            all_single_steps_set[product].append((reaction_id, reaction_dict[reaction_id][0]))
    return all_single_steps_set

def prepare_single_steps_set(reaction_filename, template_filename):
    json_data_template = load_json_data(template_filename)
    json_data_reaction = load_json_data(reaction_filename)
    all_single_steps_set = extract_reactions_and_templates(json_data_reaction, json_data_template)
    return all_single_steps_set

def get_test_starting_molecules():
    return ['a', 'b', 'c', 'd', 'e', 'j']

def get_test_single_steps_set():
    single_steps_set = {}
    single_steps_set['i'] = [(1, ['f','g']), (2, ['d', 'h'])]
    single_steps_set['f'] = [(3, ['a']), (7, ['j'])]
    single_steps_set['g'] = [(4, ['b']), (5, ['c'])]
    single_steps_set['h'] = [(6, ['e'])]
    return single_steps_set

def test_extract_route():
    starting_mols = get_test_starting_molecules()
    single_steps_set = get_test_single_steps_set()

    products = ['i']

    for product in products:
        mol_tree = MolTree(
            target_mol=product,
            known_mols=starting_mols
        )
        iter = 0
        while iter < 1000:
            m_next = mol_tree.get_fist_open_node()
            print(m_next.mol)
            if single_steps_set.__contains__(m_next.mol):
                mol_tree.expand(m_next, single_steps_set[m_next.mol])
            else:
                mol_tree.expand(m_next, None)
                print('Expansion fails on %s!' % m_next.mol)
            iter += 1

            if mol_tree.is_all_succ():
                break

        print(mol_tree.save_all_route())


def extract_route():
    starting_mols = prepare_starting_molecules('../dataset/origin_dict.csv')
    single_steps_set = prepare_single_steps_set('../dataset/uspto.reactions.json', '../dataset/uspto.templates.json')
    cnt = 0
    with open('../dataset/routes.txt', 'w') as f:
        for product in single_steps_set.keys():
            mol_tree = MolTree(
                target_mol=product,
                known_mols=starting_mols
            )
            iter = 0
            while iter < 1000:
                m_next = mol_tree.get_fist_open_node()
                if m_next is None:
                    break
                
                if single_steps_set.__contains__(m_next.mol):
                    mol_tree.expand(m_next, single_steps_set[m_next.mol])
                else:
                    mol_tree.expand(m_next, None)
                    print('Expansion fails on %s!' % m_next.mol)
                iter += 1

                if mol_tree.is_all_succ():
                    break
            
            if mol_tree.root.succ:
                routes = mol_tree.save_all_route()
                if len(routes) == 0:
                    continue
                f.write('------------------------- product %d-------------------------\n' % cnt)
                for i in range(len(routes)):
                    f.write('************************* route %d *************************\n' % i)
                    for reaction in routes[i]:
                        f.write('%s \n' % reaction)
                
                cnt += 1            
                if cnt >= 20000:
                    break
        
extract_route()