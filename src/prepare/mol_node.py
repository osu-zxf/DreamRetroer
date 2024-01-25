import numpy as np
import logging


class MolNode:
    def __init__(self, mol, parent=None, is_known=False):
        self.mol = mol
        self.parent = parent

        self.id = -1
        if self.parent is None:
            self.depth = 0
        else:
            self.depth = self.parent.depth

        self.is_known = is_known
        self.children = []
        self.succ = is_known
        self.open = True    # before expansion: True, after expansion: False
        if is_known:
            self.open = False
        
        if parent is not None:
            parent.children.append(self)

    def init_values(self, no_child=False):
        assert self.open and (no_child or self.children)

        self.succ = False
        for reaction in self.children:
            self.succ |= reaction.succ
        
        self.open = False

    def backup(self, succ):
        assert not self.is_known

        new_succ = self.succ | succ
        self.succ = new_succ

        if self.parent:
            return self.parent.backup()

    def __str__(self):
        text = '%s' % (self.mol)
        return text
    
    def serialize(self):
        '''
        return a list of str
        '''
        str_list = []
        if len(self.children) == 0:
            return [""]
        for child_reaction_node in self.children:
            if child_reaction_node.succ:
                str_list.extend(child_reaction_node.serialize())
        if len(str_list) == 0:
            return [""]
        return str_list

    def get_ancestors(self):
        if self.parent is None:
            return {self.mol}

        ancestors = self.parent.parent.get_ancestors()
        ancestors.add(self.mol)
        return ancestors