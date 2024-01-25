import os
from graphviz import Digraph
from queue import Queue

def solve_route_to_pathList(file_name):
    reaction_list = []
    node_dict = {}
    with open(file_name, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            node = line.split('>')[0]
            if node.isdigit():
                reaction_list.append(node)
            reactants = (line.split('>')[1][:-1]).split('...')
            children = []
            for reactant in reactants:
                if (reactant is not '...') and (reactant is not ''):
                    children.append(reactant)
            node_dict[node] = children
            
    pathList, stack = [], [(reaction_list[0], "")]
    while stack:
        node, pathStr = stack.pop(0)
        if not node_dict.__contains__(node):
            pathList.append(pathStr + node)
        else:
            for child in node_dict[node]:
                stack.append((child, pathStr + node + "->"))
    return pathList

def extract_reaction_chain(pathList):
    reaction_chain_list = []
    for path in pathList:
        nodes = path.split('->')
        new_path = ""
        count = 0
        for node in nodes:
            if node.isdigit():
                new_path += node + "->"
                count += 1
        if count >= 2:
            reaction_chain_list.append(new_path[:-2])
    # remove repeat
    reaction_chain_list = list(set(reaction_chain_list))
    # extract all segment
    segment_list = []
    for reaction_chain in reaction_chain_list:
        reaction_ordered_list = reaction_chain.split('->')
        if len(reaction_ordered_list) > 2:
            for i in range(len(reaction_ordered_list) - 2):
                for j in range(i+1, len(reaction_ordered_list)):
                    segment = '->'.join(reaction_ordered_list[i:j+1])
                    segment_list.append(segment)
        else:
            segment_list.append(reaction_chain)
    return list(set(segment_list)), reaction_chain_list
# for i in range(190):
#     file_name = 'syn_route_%d.txt' % (i)
#     if os.path.exists(file_name):

# # print(pathList)
# print('-----------------------')
# print(reaction_chain)

reaction_chains_list = []
file_id_list = []
for i in range(190):
    file_name = 'syn_route_%d.txt' % i
    save_file_name = 'new_route/syn_route_%d.txt' % i
    if os.path.exists(file_name):
        pathList = solve_route_to_pathList(file_name)
        segment_list, reaction_chains = extract_reaction_chain(pathList)
        reaction_chains_list.append(segment_list)
        file_id_list.append(i)

        # with open(save_file_name, 'w') as f:
        #     f.write(str(pathList))
        #     f.write('\n')
        #     f.write('----------------------------------\n')
        #     f.write(str(reaction_chains))

## write to file: which two mol syn tree have common sub-tre
# with open('inter_list.txt', 'w') as f:
#     for i in range(len(reaction_chains_list)-1):
#         for j in range(i+1, len(reaction_chains_list)):
#             inter_list = list(set(reaction_chains_list[i]).intersection(set(reaction_chains_list[j])))
#             if len(inter_list) > 0:
#                 file_str = 'route %d and route %d has repeat segment:' % (file_id_list[i], file_id_list[j])
#                 print(file_str)
#                 print(inter_list)
#                 f.write(file_str)
#                 f.write('\n')
#                 f.write(str(inter_list))
#                 f.write('\n')

inter_dict = {}
for i in range(len(reaction_chains_list)-1):
    for j in range(i+1, len(reaction_chains_list)):
        inter_list = list(set(reaction_chains_list[i]).intersection(set(reaction_chains_list[j])))
        if len(inter_list) > 0:
            if not inter_dict.__contains__(file_id_list[i]):
                inter_dict[file_id_list[i]] = []
            inter_dict[file_id_list[i]].append(file_id_list[j])

# no_common_subtree_id_list = []
# for file_id in file_id_list:
#     if not inter_dict.__contains__(file_id):
#         no_common_subtree_id_list.append(file_id)
# print('no common:', len(no_common_subtree_id_list))
# print('total:', len(file_id_list))

# with open('inter_list_summize.txt', 'w') as f:
#     for key, value in inter_dict.items():
#         f.write(str(key))
#         f.write(':')
#         f.write(str(list(set(value))))
#         f.write('\n')

import networkx as nx
import matplotlib.pyplot as plt

def my_grap():
    G = nx.Graph(my_seq='test_mol_set_graph')
    pair_list = []
    id_list = []
    for key, value in inter_dict.items():
        for item in value:
            pair_list.append((key, item))
            id_list.append(item)
    id_list = list(set(id_list))
    G.add_nodes_from(id_list)
    G.add_edges_from(pair_list)

    pos = nx.spring_layout(G)
    nx.draw(G, pos,
            # node_color='red',
            node_size=500,
            with_labels=True)
 
    # weights = nx.get_edge_attributes(G, 'weight')
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=weights)
    plt.show()

my_grap()