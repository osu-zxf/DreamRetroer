import numpy as np
import logging


class ReactionNode:
    def __init__(self, parent, idx):
        self.parent = parent
        self.depth = self.parent.depth + 1
        self.id = idx
        self.children = []
        self.succ = None    # successfully found a valid synthesis route
        self.open = True    # before expansion: True, after expansion: False
        parent.children.append(self)

    def init_values(self):
        assert self.open
        self.succ = True
        for mol in self.children:
            self.succ &= mol.succ
        self.open = False

    def backup(self):
        self.succ = True
        for mol in self.children:
            self.succ &= mol.succ

        return self.parent.backup(self.succ)

    def __str__(self):
        children_str = [child.mol for child in self.children]
        return '%s>%d>%s' % (self.parent.mol, self.id, '.'.join(children_str))
    
    def serialize(self):
        '''
        return a list of str
        '''
        if not self.succ:
            return [""]
        str_list = []
        self_str = self.__str__()
        str_list.append(self_str)
        for child_mol_node in self.children:
            child_str_list = child_mol_node.serialize()
            tmp_str_list = []
            for str1 in str_list:
                for str2 in child_str_list:
                    tmp_str_list.append(str1+","+str2)
            str_list = tmp_str_list
        return str_list
        
        