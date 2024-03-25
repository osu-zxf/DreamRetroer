import numpy as np
from queue import Queue
import logging
import networkx as nx
from graphviz import Digraph
from mol_node import MolNode
from reaction_node import ReactionNode
from syn_route import SynRoute

class MolTree:
    def __init__(self, target_mol, known_mols, zero_known_value=True):
        self.target_mol = target_mol
        self.known_mols = known_mols
        self.mol_nodes = []
        self.reaction_nodes = []

        self.root = self._add_mol_node(target_mol, None)
        self.succ = target_mol in known_mols
        self.search_status = 0

        if self.succ:
            print('Synthesis route found: target in starting molecules')

    def _add_mol_node(self, mol, parent):
        is_known = mol in self.known_mols

        mol_node = MolNode(
            mol=mol,
            parent=parent,
            is_known=is_known
        )
        self.mol_nodes.append(mol_node)
        mol_node.id = len(self.mol_nodes)

        return mol_node

    def _add_reaction_and_mol_nodes(self, mols, parent, ancestors, reaction_id):

        for mol in mols:
            if mol in ancestors:
                return

        reaction_node = ReactionNode(parent, reaction_id)
        for mol in mols:
            self._add_mol_node(mol, reaction_node)
        reaction_node.init_values()
        self.reaction_nodes.append(reaction_node)

        return reaction_node

    def expand(self, mol_node, reaction_lists):
        assert not mol_node.is_known and not mol_node.children
        if reaction_lists is None:
            mol_node.init_values(no_child=True)
            return False

        assert mol_node.open
        ancestors = mol_node.get_ancestors()
        for i in range(len(reaction_lists)):
            self._add_reaction_and_mol_nodes(reaction_lists[i][1],
                                             mol_node, ancestors, reaction_lists[i][0])

        mol_node.init_values()
        if mol_node.parent:
            mol_node.parent.backup()
        if not self.succ and self.root.succ:
            print('Synthesis route found!')
            self.succ = True

        return self.succ
    
    def get_fist_open_node(self):
        for node in self.mol_nodes:
            if node.open:
                return node
        return None
    
    def is_all_succ(self):
        all_succ = True
        for reaction_node in self.reaction_nodes:
            all_succ &= reaction_node.succ
        return all_succ

    def save_all_route(self):
        if not self.succ:
            return None
        
        syn_routes = []
        
        for route in self.root.serialize():
            syn_routes.append([x for x in route.split(',') if x])

        return syn_routes

    def viz_search_tree(self, viz_file):
        G = Digraph('G', filename=viz_file)
        G.attr(rankdir='LR')
        G.attr('node', shape='box')
        G.format = 'pdf'

        node_queue = Queue()
        node_queue.put((self.root, None))
        while not node_queue.empty():
            node, parent = node_queue.get()

            if node.open:
                color = 'lightgrey'
            else:
                color = 'aquamarine'

            if hasattr(node, 'mol'):
                shape = 'box'
            else:
                shape = 'rarrow'

            if node.succ:
                color = 'lightblue'
                if hasattr(node, 'mol') and node.is_known:
                    color = 'lightyellow'

            G.node(node.serialize(), shape=shape, color=color, style='filled')

            label = ''
            if hasattr(parent, 'mol'):
                label = '%.3f' % node.cost
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
                f.write('>')
                for reaction_nd in mol_nd.children:
                    f.write(str(reaction_nd.template_id))
                    f.write(',')
                f.write('\n')
            for reaction_nd in self.reaction_nodes:
                f.write(str(reaction_nd.template_id))
                f.write('>')
                for mol_nd in reaction_nd.children:
                    f.write(mol_nd.mol)
                    f.write('...')
                f.write('\n')
