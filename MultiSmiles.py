from CGRtools import ReactionContainer
from collections import OrderedDict
from operator import itemgetter

class Get_MS():
    def __init__(self, reactions):
        self.reactions = reactions

    def get_rules(self, reaction):
        rules = []
        for n, partial_reaction in enumerate(reaction.enumerate_centers()):
            if n == 5:
                break
            reaction_center = set(partial_reaction.extended_centers_list[0])
            bare_reaction_center = set(partial_reaction.compose().center_atoms)
            reactants = []
            products = []
            for mol in partial_reaction.reactants:
                group_atoms = reaction_center.intersection(mol)
                if group_atoms:
                    group = mol.substructure(group_atoms, as_query=True)
                    for i in group_atoms.difference(bare_reaction_center):
                        group._neighbors[i] = ()
                    reactants.append(group)
            for mol in partial_reaction.products:
                group_atoms = reaction_center.intersection(mol)
                if group_atoms:
                    group = mol.substructure(group_atoms, as_query=True)
                    for i in group_atoms.difference(bare_reaction_center):
                        group._neighbors[i] = ()
                    products.append(group)
            rule = ReactionContainer(reactants=reactants, products=products,
                                     meta=reaction.meta)
            rules.append(rule)
        return rules

    def get_nodes(self):
        nodes = dict()
        for reaction in self.reactions:
            if len(reaction.reactants) == 1:
                nodes[str(reaction.reactants[0])] = [None, str(reaction.products[0]), self.get_rules(reaction)[0]]
            elif len(reaction.reactants) > 1:
                nodes[str(reaction.reactants[0])] = [str(reaction.reactants[1]), str(reaction.products[0]), self.get_rules(reaction)[0]]
                nodes[str(reaction.reactants[1])] = [str(reaction.reactants[0]), str(reaction.products[0]), self.get_rules(reaction)[0]]
        return nodes

    def get_target(self):
        nodes = self.get_nodes()
        for p in [x for y in nodes.values() for x in y]:
            if p not in [i for i in nodes.keys()] and p is not None:
                target_mol = p
                return target_mol

    def _get_child(self, node):
        childs = list()
        for k,l_v in self.get_nodes().items():
            if node == l_v[1]:
                childs.append(k)
                childs.extend(self._get_child(k))
        return childs

    def _get_depth(self, node_start, count=0):
        nodes = self.get_nodes()
        _, v1 = nodes[node_start]
        if self.get_target() != v1:
            if v1 in nodes:
                count += 1
                count = self._get_depth(v1, count)
        return count

    def fit(self):
        smiles_all = []
        nodes = self.get_nodes()
        nodes_copy, target_mol = nodes.copy(), self.get_target()
        for target in [k for k, l_v in nodes_copy.items() if target_mol in l_v]:
            count = 0
            all_child = []
            smile_branch = []
            for p in self._get_child(target):
                all_child.append(p)
            branch = {child: self._get_depth(child) for child in all_child}
            st_branch = OrderedDict(sorted(branch.items(), key=itemgetter(1), reverse=True))
            for mol in st_branch:
                if mol in nodes_copy:
                    rd = str(nodes_copy[mol][2])
                    if count == 0:
                        if nodes_copy[mol][0] is None:
                            smile = str(mol) + str({rd}) + '^'
                            del nodes_copy[mol]
                        else:
                            smile = str(mol) + '.' + str(nodes_copy[mol][0]) + str({rd}) + '+'
                            neigh = nodes_copy[mol][0]
                            del nodes_copy[mol], nodes_copy[neigh]
                        count += 1

                    else:
                        if nodes_copy[mol][0] is None:
                            smile = str({rd}) + '^'
                            del nodes_copy[mol]
                        else:
                            if mol in [v[1] for v in nodes_copy.values()]:
                                smile = str(mol) + str({rd}) + '+'
                            if mol not in [v[1] for v in nodes_copy.values()]:
                                smile = str(nodes_copy[mol][0]) + str({rd}) + '+'
                            if nodes_copy[mol][1] != target_mol:
                                neigh = nodes_copy[mol][0]
                                del nodes_copy[mol], nodes_copy[neigh]
                            else:
                                del nodes_copy[mol]
                    smile_branch.append(smile)

            rd = str(nodes_copy[mol][2])
            smiles_all.extend(smile_branch)
            if nodes[target][0] == None:  # one branch
                my_smile = ''.join(smiles_all) + str({rd}) + '^'
            elif len(self._get_child(target)) == 0:
                my_smile = ''.join(smiles_all) + str(target) + str({rd}) + '+'
            else:
                my_smile = ''.join(smiles_all) + str({rd}) + '+'

        return my_smile