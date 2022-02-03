from CGRtools import ReactionContainer
from collections import OrderedDict
from operator import itemgetter
import string
import random

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
            reactants, products = list(), list()
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
        for p in [x for y in nodes.values() for x in y[:2]]:
            if p not in [i for i in nodes.keys()] and p is not None:
                target_mol = p
                return target_mol

    def _get_child(self, node):
        childs = list()
        for k, l_v in self.get_nodes().items():
            if node == l_v[1]:
                childs.append(k)
                childs.extend(self._get_child(k))
        return childs

    def _get_depth(self, node_start, count=0):
        nodes = self.get_nodes()
        _, p, _ = nodes[node_start]
        if self.get_target() != p:
            if p in nodes:
                count += 1
                count = self._get_depth(p, count)
        return count

    def fit(self):
        smiles_all = []
        nodes = self.get_nodes()
        nodes_copy, target_mol = nodes.copy(), self.get_target()
        for target in [k for k, l_v in nodes_copy.items() if target_mol in l_v]:
            count = 0
            all_child, smile_branch = list(), list()
            for p in self._get_child(target):
                all_child.append(p)
            branch = {child: self._get_depth(child) for child in all_child}
            st_branch = OrderedDict(sorted(branch.items(), key=itemgetter(1), reverse=True))
            for mol in st_branch:
                if mol in nodes_copy:
                    rule = str(nodes_copy[mol][2])
                    if count == 0:
                        if nodes_copy[mol][0] is None:
                            smile_branch.append(str(mol) + str({rule}) + '^')
                            del nodes_copy[mol]
                        else:
                            smile_branch.append(str(mol) + '.' + str(nodes_copy[mol][0]) + str({rule}) + '+')
                            neigh = nodes_copy[mol][0]
                            del nodes_copy[mol], nodes_copy[neigh]
                        count += 1
                    else:
                        if nodes_copy[mol][0] is None:
                            smile_branch.append(str({rule}) + '^')
                            del nodes_copy[mol]
                        else:
                            if mol in [v[1] for v in nodes_copy.values()]:
                                smile_branch.append(str(mol) + str({rule}) + '+')
                            if mol not in [v[1] for v in nodes_copy.values()]:
                                smile_branch.append(str(nodes_copy[mol][0]) + str({rule}) + '+')
                            if nodes_copy[mol][1] != target_mol:
                                neigh = nodes_copy[mol][0]
                                del nodes_copy[mol], nodes_copy[neigh]
                            else:
                                del nodes_copy[mol]

            smiles_all.extend(smile_branch)
            rule = str([l_v[2] for k, l_v in nodes.items() if target_mol in l_v][0])
            if nodes[target][0] is None:  # one branch
                return ''.join(smiles_all) + str({rule}) + '^'
            elif len(self._get_child(target)) == 0:
                return ''.join(smiles_all) + str(target) + str({rule}) + '+'
            else:
                return ''.join(smiles_all) + str({rule}) + '+'

class parser_MS():
    def __init__(self, smile):
        self.smile = smile

    def _split_smile(self, sm_t):
        for i in sm_t[1:]:
            if i == '+' or i == '^':
                l = sm_t[1:][:sm_t[1:].index(i)].split('{')
                if len(l) == 2:
                    return (l[0], l[1][:-1], sm_t[sm_t[1:].index(i) + 1:])
                if len(l) == 3:
                    return (l[0], l[1][:-1], sm_t[sm_t[1:].index(i) + 1:])

    def fit(self):
        d = OrderedDict()
        flag = True
        for num, part in enumerate(self.smile.partition('}'), 1):
            if num == 1:
                mol_rule = part.split('{')
                par = ''.join([random.choice(string.ascii_letters + string.digits) for n in range(2)])
                if self.smile.partition('}')[2].startswith('^'):
                    d[mol_rule[0]] = (None, par, mol_rule[1])
                elif self.smile.partition('}')[2].startswith('+'):
                    mols = mol_rule[0].split('.')
                    d[mols[0]] = (mols[1], par, mol_rule[1])
                    d[mols[1]] = (mols[0], par, mol_rule[1])
            elif part != '}':
                next_part = part
                while flag:
                    out = self._split_smile(next_part)
                    reag, rule, next_part = out
                    par = ''.join([random.choice(string.ascii_letters + string.digits) for n in range(2)])
                    if len(reag) == 0:  # no reagent
                        if next_part.startswith('^'):
                            d[[_ for _ in d.values()][-1][1]] = (None, par, rule)
                        if next_part.startswith('+'):
                            point = par
                            d[[_ for _ in d.values()][-1][1]] = (point, par, rule)  # для соединения веток
                            d[point] = ([_ for _ in d.values()][-1][1], par, rule)  # для соединения веток
                    else:
                        if len(reag) == 1:
                            if next_part.startswith('^'):  # new b
                                d[reag] = (None, par, rule)
                            if next_part.startswith('+'):
                                d[reag] = ([_ for _ in d.values()][-1][1], par, rule)
                                d[[_ for _ in d.values()][-1][0]] = (reag, par, rule)
                        if len(reag) == 2:  # new b
                            d[reag[0]] = (reag[1], par, rule)
                            d[reag[1]] = (reag[0], par, rule)

                    if next_part == '+' or next_part == '^':
                        flag = False
        return d