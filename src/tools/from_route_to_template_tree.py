import pickle
import pandas as pd
import json
from rdkit import Chem

def load_json_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data


def save_data_to_pkl(data, saveto):
    with open(saveto, 'wb') as f:  # open a text file
        pickle.dump(data, f) # serialize the list

## target product
def prepare_target_product_list():
    ### load route
    routes = pickle.load(open('../dataset/routes_possible_test_hard.pkl', 'rb'))
    print('total load %d routes' % len(routes))
    ### extract product
    final_product_list = []
    all_product_list = []
    for route in routes:
        target_product = route[0].split('>')[0]
        final_product_list.append(target_product)
        for reaction in route:
            inter_product = reaction.split('>')[0]
            all_product_list.append(inter_product)
    print('target product count: %d' % (len(final_product_list)))
    print('inter product count: %d' % (len(all_product_list)))
    return final_product_list, all_product_list, routes

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
            reaction = ([clear_atom_mapping(x) for x in item["reactants"].split('.')], clear_atom_mapping(item["products"]))
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
            all_single_steps_set[product].append((reaction_id, reaction_dict[reaction_id][0], rxn_template))
    return all_single_steps_set

def prepare_single_steps_set(reaction_filename, template_filename):
    json_data_template = load_json_data(template_filename)
    json_data_reaction = load_json_data(reaction_filename)
    all_single_steps_set = extract_reactions_and_templates(json_data_reaction, json_data_template)
    return all_single_steps_set

## source dataset
def prepare_source_product_list_from_reaction_set():
    ### extract product from source dataset
    # single_steps_set = prepare_single_steps_set('../dataset/uspto.reactions.json', '../dataset/uspto.templates.json')
    # source_product_list = list(single_steps_set.keys())
    source_product_list = []
    with open('../dataset/uspto_solved_reaction.csv', 'r') as f:
        while True:
            # Get next line from file
            line = f.readline()
            # if line is empty
            # or end of file is reached
            if not line:
                break
            source_product_list.append(line.strip())
    return source_product_list


# Python program to illustrate the intersection
# of two lists using set() and intersection()
def Intersection(lst1, lst2):
    return set(lst1).intersection(lst2)

source_products = prepare_source_product_list_from_reaction_set()
target_products, all_inter_products, routes = prepare_target_product_list()
print("target overlap with inter product count: %d" % (len(Intersection(target_products, all_inter_products))))
print('target products intersection count: %d' % (len(Intersection(target_products, source_products))))
print('all inter products intersection count: %d' % (len(Intersection(all_inter_products, source_products))))

# routes_start_with_valid_product = []
# source_products_set = set(source_products)
# for i in range(len(target_products)):
#     if i % 10000 == 0:
#         print(i)
#     if target_products[i] in source_products_set:
#         routes_start_with_valid_product.append(routes[i])
# save_data_to_pkl(routes_start_with_valid_product, "../dataset/valid_routes.pkl")
