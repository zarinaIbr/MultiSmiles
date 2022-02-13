from CGRtools import ReactionContainer, SMILESRead
from collections import OrderedDict
from operator import itemgetter
from random import sample
from string import digits
from CGRtools.reactor import Reactor
from tempfile import mkdtemp
from shutil import rmtree
from pathlib import Path

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
            rule = self.get_rules(reaction)[0]
            if len(reaction.reactants) == 1:
                nodes[str(reaction.reactants[0])] = [None, str(reaction.products[0]), rule]
            elif len(reaction.reactants) > 1:
                nodes[str(reaction.reactants[0])] = [str(reaction.reactants[1]), str(reaction.products[0]), rule]
                nodes[str(reaction.reactants[1])] = [str(reaction.reactants[0]), str(reaction.products[0]), rule]
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
        smiles_all, rules = [], []
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
                    rule = nodes_copy[mol][2]
                    if count == 0:
                        if nodes_copy[mol][0] is None:
                            smile_branch.append(str(mol) + str({str(rule)}) + '^')
                        else:
                            smile_branch.append(str(mol) + '.' + str(nodes_copy[mol][0]) + str({str(rule)}) + '+')
                            neigh = nodes_copy[mol][0]
                            del nodes_copy[neigh]
                        rules.append(rule)
                        del nodes_copy[mol]
                        count += 1
                    else:
                        if nodes_copy[mol][0] is None:
                            smile_branch.append(str({str(rule)}) + '^')
                            rules.append(rule)
                            del nodes_copy[mol]
                        else:
                            if mol in [v[1] for v in nodes_copy.values()]:
                                smile_branch.append(str(mol) + str({str(rule)}) + '+')
                            if mol not in [v[1] for v in nodes_copy.values()]:
                                smile_branch.append(str(nodes_copy[mol][0]) + str({str(rule)}) + '+')
                            rules.append(rule)
                            if nodes_copy[mol][1] != target_mol:
                                neigh = nodes_copy[mol][0]
                                del nodes_copy[mol], nodes_copy[neigh]
                            else:
                                del nodes_copy[mol]

            smiles_all.extend(smile_branch)
        rule = [l_v[2] for k, l_v in nodes.items() if target_mol in l_v][0]
        rules.append(rule)
        if len([k for k, l_v in self.get_nodes().items() if target_mol in l_v]) == 1:  # one branch
            return ''.join(smiles_all) + str({str(rule)}) + '^', rules
        elif any(len(self._get_child(t_child)) == 0 for t_child in nodes_copy):  # one molecule in the second branch
            return ''.join(smiles_all) + str([t_child for t_child in nodes_copy if len(self._get_child(t_child)) == 0][0]) + str({str(rule)}) + '+', rules
        else:
            return ''.join(smiles_all) + str({str(rule)}) + '+', rules


class parser_MS():
    def __init__(self, smile_s):
        self.smile = smile_s[0]
        self.rules = smile_s[1].copy()

    def generate(self, rules, reactants):
        reactors = [Reactor(x, delete_atoms=True) for x in rules]
        seen = set()
        queue = [(r(reactants, False), [0] * len(rules)) for r in reactors]
        while queue:
            r, n = queue.pop(0)
            for i in r:
                if i in seen:
                    continue
                seen.add(i)
                yield i

    def _split_smile(self, sm_t):
        for i in sm_t[1:]:
            if i == '+' or i == '^':
                l = sm_t[1:][:sm_t[1:].index(i)].split('{')
                if len(l) == 2:
                    return l[0], sm_t[sm_t[1:].index(i) + 1:]
                if len(l) == 3:
                    return l[0], sm_t[sm_t[1:].index(i) + 1:]

    def gen_graph(self):
        d = OrderedDict()
        flag = True
        for num, part in enumerate(self.smile.partition('}'), 1):
            if num == 1:
                rule = self.rules.pop(0)
                par = ''.join(sample(digits, 4))
                if self.smile.partition('}')[2].startswith('^'):
                    d[part.split('{')[0]] = [None, par, rule]
                elif self.smile.partition('}')[2].startswith('+'):
                    mols = part.split('{')[0].split('.')
                    d[mols[0]] = [mols[1], par, rule]
                    d[mols[1]] = [mols[0], par, rule]
            elif part != '}':
                next_part = part
                while flag:
                    reag, next_part = self._split_smile(next_part)
                    rule = self.rules.pop(0)
                    par = ''.join(sample(digits, 4))
                    if len(reag) == 0:  # no reagent
                        if next_part.startswith('^'):
                            d[list(d.values())[-1][1]] = [None, par, rule]
                        elif next_part.startswith('+'):
                            d[list(d.values())[-2][1]] = [list(d.values())[-1][1], par, rule]  # для соединения веток
                            d[list(d.values())[-1][0]] = [list(d.keys())[-1], par, rule]  # для соединения веток
                    else:  # reag
                        if '.' not in reag:
                            if next_part.startswith('^'):  # new b
                                d[reag] = [None, par, rule]
                            elif next_part.startswith('+'):
                                d[reag] = [list(d.values())[-1][1], par, rule]
                                d[list(d.values())[-1][0]] = [reag, par, rule]
                        if '.' in reag:  # new b
                            d[reag[0]] = [reag[1], par, rule]
                            d[reag[1]] = [reag[0], par, rule]
                    if next_part == '+' or next_part == '^':
                        flag = False
        return d

    def get_p(self, sm_mol, d):
        work_dir = Path(mkdtemp(prefix='sm_'))
        idx = d[sm_mol][1]
        file = work_dir / 'sm_mol'
        with open(file, 'w') as f:
            f.write(sm_mol)
        mol = SMILESRead(Path(file)).read()[0]
        mol.clean2d()
        reaction = [_ for _ in self.generate([d[sm_mol][2]], [mol])]
        rmtree(work_dir)
        if reaction:
            return idx, reaction
        return

    def fit(self):
        d = self.gen_graph()
        mols = [m for m in d if not m.isdigit()]
        for num, mol in enumerate(mols):
            args = self.get_p(mol, d)
            if args:
                idx, product = args[0], args[1][0].products[0]
                d[mol][1] = product
                if idx in d.keys():
                    d[product] = d[idx]
                    del d[idx]
                for mol in d:
                    if idx in d[mol]:
                        d[mol][d[mol].index(idx)] = product
        return d