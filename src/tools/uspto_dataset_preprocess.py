import json

def load_json_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
        print(len(data))
        print(data[0])
    return data

def process_template_file(filename, saveto="", nums=10):
    json_data = load_json_data(filename)
    with open(saveto, 'w') as fw:
        json.dump(json_data[0:nums], fw, indent=2)

def process_reaction_file(filename, saveto="", nums=10):
    json_data = load_json_data(filename)

    with open(saveto, 'w') as fw:
        json.dump(json_data[0:nums], fw, indent=2)

def preprocess_to_subset():
    filename = '../dataset/uspto.templates.json'
    saveto = '../dataset/uspto.templates_subset.json'
    process_template_file(filename, saveto, 100)
    filename = '../dataset/uspto.reactions.json'
    saveto = '../dataset/uspto.reactions_subset.json'
    process_reaction_file(filename, saveto)

def extract_templates(json_data):
    templates = []
    for item in json_data:
        if "reaction_smarts" in item:
            templates.append(item["reaction_smarts"])
    return templates

def extract_target_molecules(json_data):
    products = []
    for item in json_data:
        if "products" in item:
            if len(item["products"].split('.')) == 1:
                products.append(item["products"])
    return products

def extract_reactions_and_templates(json_data_reaction, json_data_template):
    results = []
    reaction_dict = {}
    template_dict = {}
    for item in json_data_reaction:
        if "_id" in item and "reactants" in item and "products" in item and "spectators" in item:
            reaction = item["reactants"] + ">" + item["spectators"] + ">" + item["products"]
            reaction_dict[item["_id"]] = reaction

    for item in json_data_template:
        if "reaction_id" in item and "reaction_smarts" in item and "products" in item and "reactants" in item:
            retro_template = item["reaction_smarts"]
            # reaction = item["reactants"] + ">>" + item["products"]
            template_dict[item["reaction_id"]] = retro_template

    for reaction_id, rxn_smiles in reaction_dict.items():
        if template_dict.__contains__(reaction_id):
            results.append((reaction_id, rxn_smiles, template_dict[reaction_id]))
    return results

def extract_only_templates(filename, saveto=""):
    json_data = load_json_data(filename)
    templates = extract_templates(json_data)
    templates_set = set(templates)
    print('extract done, total %d items' % len(templates))
    print('after filter, start to save, total %d items' % len(templates_set))
    # with open(saveto, 'w') as f:
    #     for item in templates_set:
    #         f.write('%s\n' % item)

def extract_only_target_molcules(filename, saveto=""):
    json_data = load_json_data(filename)
    products = extract_target_molecules(json_data)
    products_set = set(products)
    print('extract done, total %d items' % len(products))
    print('after filter, start to save, total %d items' % len(products_set))
    # with open(saveto, 'w') as f:
    #     for item in products_set:
    #         f.write('%s\n' % item)

def extract_details_to_file(reaction_filename, template_filename, saveto=''):
    json_data_template = load_json_data(template_filename)
    json_data_reaction = load_json_data(reaction_filename)
    results = extract_reactions_and_templates(json_data_reaction, json_data_template)
    with open(saveto, 'w') as f:
        f.write("id,rxn_smiles,retro_templates\n")
        id = 0
        for reaction_pair in results:
            reaction_id, reaction, retro_template = reaction_pair[0], reaction_pair[1], reaction_pair[2]
            f.write('%d,%s,%s\n' % (reaction_id, reaction, retro_template))
            id += 1

# extract_only_templates('../dataset/uspto.templates.json', '../dataset/uspto.templates.txt')
# extract_only_target_molcules('../dataset/uspto.reactions.json', '../dataset/uspto.reactions.txt')
extract_details_to_file('../dataset/uspto.reactions.json', '../dataset/uspto.templates.json', '../dataset/uspto_reactions.csv')
# preprocess_to_subset()

