import numpy as np
from queue import Queue
import logging
import networkx as nx
from graphviz import Digraph
from dreamretroer.algorithm.mol_node import MolNode
from dreamretroer.algorithm.reaction_node import ReactionNode
from dreamretroer.algorithm.syn_route import SynRoute
# from queue import Queue

from dreamretroer.algorithm.priority_queue import PriorityQueue

class SearchTree:
    def __init__(self, target_mol, building_blocks, value_fn, zero_known_value=True):
        self.target_mol = target_mol
        # building_blocks
        self.building_blocks = building_blocks
        # Experience Guidance Network
        self.value_fn = value_fn
        self.zero_known_value = zero_known_value
        self.mol_nodes = []
        self.mol_smiles = {}
        self.overlap_count = 0
        self.reaction_nodes = []
        self.root = self._add_mol_node(target_mol, None)
        self.succ = target_mol in building_blocks

        self.expand_model_call = 0
        self.value_model_call = 0
        if self.succ:
            logging.info('Synthesis route found: target in building blocks')
        

    def _add_mol_node(self, mol, parent):
        # already in building blocks
        is_known = mol in self.building_blocks

        mol_node = MolNode(
            mol=mol,
            parent=parent,
            is_known=is_known,
            zero_known_value=self.zero_known_value
        )
        self.mol_nodes.append(mol_node)
        mol_node.id = len(self.mol_nodes)
        self.mol_smiles[mol] = len(self.mol_nodes)-1

        return mol_node

    def _add_reaction_and_mol_nodes(self, pro, mols, parent, template, init_value, ancestors=[]):
        assert pro >= 0
        # repeat : ignored
        for mol in mols:
            if mol in ancestors:
                return None

        reaction_node = ReactionNode(parent, pro, template, init_value)
        succe = True
        for mol in mols:
            if self.mol_smiles.__contains__(mol):
                mol_id = self.mol_smiles[mol]
                reaction_node.add_mol_node(self.mol_nodes[mol_id])
                self.overlap_count+=1
                succe = succe & self.mol_nodes[mol_id].succ
            else:
                node = self._add_mol_node(mol, reaction_node)
                succe = succe & node.succ

            
        # all mol are successful
        if succe:
            reaction_node.set_success()
        self.reaction_nodes.append(reaction_node)
        reaction_node.id = len(self.reaction_nodes)

        return reaction_node

    # expand phase
    def expand(self, mol_node, reactant_lists, pros, templates):
        assert not mol_node.is_known #and not mol_node.children
        if pros is None:
            # there is no decomposition action: mol node failed -> -inf
            mol_node.value = np.NINF
            # already expanded
            mol_node.open = False
            return
        else :
            #  mol node has some decomposition actions
            # assert mol_node.open
            ancestors = mol_node.get_ancestors()
            succe = False
            for i in range(len(pros)):
                # experience guidance network
                init_value = self.value_fn(mol_node.mol, templates[i])
                self.value_model_call += 1
                node = self._add_reaction_and_mol_nodes(pros[i], reactant_lists[i], mol_node, templates[i], init_value, ancestors)

                if node != None:
                    succe = succe|node.succ
                
            # print('now total mol nodes: %d, overlap mol nodes: %d, ratio: %.2f' % (len(self.mol_nodes), self.overlap_count, (float(self.overlap_count)/(self.overlap_count+len(self.mol_nodes)))))
            max_depth = -1
            max_parent_num = -1
            for mol in self.mol_nodes:
                max_depth = max(max_depth, mol.depth)
                max_parent_num = max(max_parent_num, len(mol.parents))
            # print('max depth:', max_depth, ', max_parent_num:', max_parent_num)
            # all are invalid actions : mol node fails -> -inf and open = false: has been expanded
            if len(mol_node.children) == 0:
                mol_node.value = np.NINF
                mol_node.open = False
                return
            # there is a successful action
            if succe:
                # then  mol node is successful  (or node)
                mol_node.succ = True
                mol_node.value = 10.0
            else:
                # mol node value = max child.Q
                valid_children = mol_node.children_list()
                valid_children_Q = [valid_child.Q_value for valid_child in valid_children]
                max_child_Q = np.max(valid_children_Q)

                mol_node.value = max_child_Q
            # has been expanded
            mol_node.open = False
            return

    # update phase
    def update(self, mol_node):
        mol_queue = PriorityQueue()
        mol_queue.push(mol_node.mol, 0)
        update_count = 0
        while not mol_queue.is_empty():
            node_depth, mol_node_smile = mol_queue.pop()
            mol_node_tmp = self.mol_nodes[self.mol_smiles[mol_node_smile]]
            mol_node_tmp.backup()
            push_count = 0
            modify_count = 0
            grand_parents = []
            for parent in mol_node_tmp.parents:
                if parent is not None:
                    if not mol_queue.contains(parent.parent.mol):
                        mol_queue.push(parent.parent.mol, node_depth+1)
                        grand_parents.append(parent.parent.mol)
                        push_count+=1
                    else:
                        mol_queue.update_priority(parent.parent.mol, node_depth+1)
                        modify_count+=1
                        grand_parents.append(parent.parent.mol)
            update_count += 1
            debug_str = 'mol_smile:%s, depth:%d, update phase count:%d, push count:%d, modify count:%d, total:%d' % (mol_node_smile, node_depth, update_count, push_count, modify_count, mol_queue.size())
            # logging.info(debug_str)
            # logging.info(','.join(grand_parents))
        # print('*****update phase final count:%d' % update_count)
        if self.root.succ:
            self.succ = True
        return self.succ

    # get the synthetic experience
    def get_tree_Qvalue(self):
        experience = []
        for reaction_node in self.reaction_nodes:
            # valid but not succ and updata count<3 and Q>0  not collected
            if reaction_node.count < 3 and reaction_node.valid == True and reaction_node.succ == False and reaction_node.Q_value > 0.0:
                continue
            mol = reaction_node.parent.mol
            template = reaction_node.template
            Q_value = reaction_node.Q_value
            # succ = reaction_node.succ
            # childs = reaction_node.children
            # child_list = [(child.value, child.succ, child.open) for child in childs]
            # experience( (m,T),Q_value)
            experience.append((mol, template, Q_value))
        return experience

    def get_solution_route(self,iterations):
        if not self.succ:
            return None

        syn_route = SynRoute(
            target_mol=self.root.mol,
        )

        mol_queue = Queue()
        mol_queue.put(self.root)
        while not mol_queue.empty():
            mol = mol_queue.get()
            if mol.is_known:
                syn_route.set_value(mol.mol, mol.value)
                continue

            reaction_nodes = mol.children_list()

            reaction_succ = [reaction_node.succ for reaction_node in reaction_nodes]
            best_reaction = reaction_nodes[np.argmax(reaction_succ)]

            assert best_reaction.succ == True

            reactants = []
            for reactant in best_reaction.children:
                mol_queue.put(reactant)
                reactants.append(reactant.mol)

            syn_route.add_reaction(
                mol=mol.mol,
                value=mol.value,
                template=best_reaction.template,
                prob = best_reaction.pro,
                Q0 = best_reaction.value,
                reactants=reactants,
            )

        return syn_route

    # from <https://github.com/binghong-ml/retro_star/>
    def viz_search_tree(self, viz_file):
        G = Digraph('G', filename=viz_file)
        G.attr(rankdir='LR')
        G.attr('node', shape='box')
        G.format = 'pdf'

        node_queue = Queue()
        node_queue.put((self.root, None))
        while not node_queue.empty():
            node, parent = node_queue.get()

            color = 'lightgrey'
            if hasattr(node, 'mol'):
                shape = 'box'
            else:
                shape = 'rarrow'

            if node.succ:
                color = 'red'
                if hasattr(node, 'mol') and node.is_known:
                    color = 'lightyellow'
                print(color)
            G.node(node.serialize(), shape=shape, color=color, style='filled')

            label = ''
            if hasattr(parent, 'mol'):
                label = '%.3f' % node.Q_value
                # label = 'pro:%.3f  Q:%,3f'%(node.pro,node.Q_value)

            if parent is not None:
                G.edge(parent.serialize(), node.serialize(), label=label)

            if node.children is not None:
                for c in node.children:
                    node_queue.put((c, node))

        G.render()

    def save_to_file(self, file_name):
        with open(file_name, 'w') as f:
            for mol_nd in self.mol_nodes:
                f.write(mol_nd.mol)
                f.write('>>>')
                for reaction_nd in mol_nd.children:
                    f.write(str(reaction_nd.template))
                    f.write(',')
                f.write('\n')
            for reaction_nd in self.reaction_nodes:
                f.write(str(reaction_nd.template))
                f.write('>>>')
                for mol_nd in reaction_nd.children:
                    f.write(mol_nd.mol)
                    f.write('...')
                f.write('\n')

