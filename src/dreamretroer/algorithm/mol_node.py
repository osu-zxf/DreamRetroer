import numpy as np
import logging


class MolNode:
    def __init__(self, mol, init_value=np.inf, parent=None, is_known = False, zero_known_value = True):
        self.mol = mol
        # Vm = 1 if m is known
        # Vm = -inf if m can not be expanded
        # Vm = argmax Q
        self.value = init_value
        self.parents = [parent]

        self.id = -1
        if parent is None:
            self.depth = 0
        else: 
            self.depth = parent.depth

        self.is_known = is_known
        self.children = []
        self.succ = is_known
        # before expansion: open = True and VM = inf, after expansion: False and has Vm
        self.open = True
        # already in B
        if is_known:
            self.open = False
            if zero_known_value:
                self.value = 10.0

        self.ancestors = {self.mol}
        if parent is not None:
            parent.children.append(self)
            for ances in parent.parent.get_ancestors():
                self.ancestors.add(ances)

    def add_parent_reaction_node(self, reaction_node):
        self.parents.append(reaction_node)
        for ances in reaction_node.parent.get_ancestors():
            self.ancestors.add(ances)

    def children_list(self):
        if len(self.children) != 0:
            # all_valid reactions
            valid_child = [child.valid for child in self.children]
            valid_index = np.where(np.array(valid_child)==1)
            if len(valid_index[0]) == 0:
                return None
            return [self.children[i] for i in valid_index[0]]
        else :
            return None
    
    def unsucc_children_list(self):
        if len(self.children) != 0:
            # all_valid reactions
            valid_child = [child.valid for child in self.children]
            valid_index = np.where(np.array(valid_child)==1)
            if len(valid_index[0]) == 0:
                return None
            child_list = []
            for i in valid_index[0]:
                if not self.children[i].succ:
                    child_list.append(self.children[i])
            if len(child_list)>0:
                return child_list
            else:
                return None
        else :
            return None
        
    def check_invalid(self):
        valid_children = self.children_list()
        if valid_children is None:
            self.value = np.NINF
            return True
        return False

    def v_self(self):
        return self.value

    def serialize(self):
        text = '%d | %s' % (self.id, self.mol)
        return text
    
    def backup(self):
        assert not self.is_known
        
        if self.value == np.NINF:
            for parent in self.parents:
                if parent is not None:
                    parent.set_invaild()
                    parent.parent.check_invalid()
                    # parent.backup()
            return

        updated = True
        old_Vm = self.value

        valid_children = self.children_list()
        if valid_children is None:
            self.value = np.NINF
        else:
            succe = False
            for child in valid_children: 
                if child.succ == True:
                    succe = True
                    break
            if succe :
                self.value = 10.0
                self.succ = True
            else:
                valid_children_Q = [valid_child.Q_value for valid_child in valid_children]
                max_child_Q = np.max(valid_children_Q)
                if abs(max_child_Q - old_Vm) < 0.01:
                    # no more update
                    updated = False
                else:
                    self.value = max_child_Q

        if updated and len(self.parents) != 0:
            for parent in self.parents:
                if parent is not None:
                    parent.backup()

    def get_ancestors(self):
        return self.ancestors
        if self.parent is None:
            return {self.mol}
        ancestors = self.parent.parent.get_ancestors()
        ancestors.add(self.mol)
        return ancestors