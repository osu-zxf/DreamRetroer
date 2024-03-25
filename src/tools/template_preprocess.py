import json
from tqdm import tqdm, trange

def load_json_data(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
        print(len(data))
        print(data[0])
    return data

def extract_templates(json_data):
    templates = []
    for item in json_data:
        if "reaction_smarts" in item:
            templates.append(item["reaction_smarts"])
    return templates


## uspto templates
json_data = load_json_data('../dataset/uspto.templates.json')
uspto_templates = extract_templates(json_data)
uspto_templates = list(set(uspto_templates))
print('load %d uspto templates' % len(uspto_templates))

## retro_star templates
template_rule_path = '../dataset/template_rules_1.dat'
template_rules = []
with open(template_rule_path, 'r') as f:
    for i, l in tqdm(enumerate(f), desc='template rules'):
        rule = l.strip()
        ## the template in retro star always format as (xxxx)>>a.b
        ## however, in uspto json, it look like xxxx>>a.b
        rule_str_list = rule.split(">>")
        rule = rule_str_list[0][1:-1]+">>"+rule_str_list[1]
        template_rules.append(rule)
template_rules = list(set(template_rules))
print('load %d retrostar templates' % len(template_rules))
print(len(list(set(uspto_templates).intersection(set(template_rules)))))
