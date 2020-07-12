# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell Lim
# University of California, San Francisco

from warnings import warn

from .. import syslayer
from ..functions import Selector, Constant


class BioLayer:

    def __init__(self):
        self.compartments = []  # unique compartment names
        self.compartment_linkages = []  # tuples of compartment names

        self.tx_cells = []
        # tx celltypes represented as dictionaries with:
        # 'name' : name
        # 'states' : ['state1', 'state2', ...]
        # 'state_linkages' : [statelink1, statelink2, ...]
        # 'daughter_state' : [('state1', 0), ('state2', 0), ...]
        # 'cytokine_linkages' : [cytokineLink1, cytokineLink2, ...]
        # 'killing_linkages' : [killLink1, killLink2, ...]

        self.cells = []
        # normal celltypes are representated as dictionaries with:
        # 'name' : name
        # 'compartment' : compartment_name

        self.cytokines = []
        # cytokines represented as dictionaries with:
        # 'name' : name

        self.interactions = []
        # custom interactions may be defined using a dictionary representation
        # 'a' : SingleSelector
        # 'b' : SingleSelector
        # 'func' : sympy expression of constants and single selectors

    # Adding Elements - Add instance to relevant dict
    def add_compartment(self, name):
        self.compartments.append(name)

    def link_compartments(self, a, b):
        self.compartment_linkages.append((a, b))

    def add_tx_cells(self, name, states, state_linkages=[], daughter_state=[], cytokine_linkages=[], killing_linkages=[]):
        c = {}
        c['name'] = name
        c['states'] = states
        c['state_linkages'] = state_linkages
        c['daughter_state'] = daughter_state
        c['cytokine_linkages'] = cytokine_linkages
        c['killing_linkages'] = killing_linkages
        self.tx_cells.append(c)

    def add_cells(self, name, compartment):
        c = {}
        c['name'] = name
        c['compartment'] = compartment
        self.cells.append(c)

    def add_cytokine(self, name):
        c = {}
        c['name'] = name
        self.cytokines.append(c)

    def add_interaction(self, a, b, func):
        i = {}
        i['a'] = a
        i['b'] = b
        i['func'] = func
        self.interactions.append(i)

    def validate_selector(self, selector):
        # TODO: Implement selector validation
        return True

    # SELECTOR Getters - get a selector for something
    def get_tx_celltype(self, name):
        kind = 'tx_celltype'
        cs = {}
        cs['type'] = kind
        cs['name'] = name

        selector_name = '[%s].[%s]' % (kind, name)

        return Selector(selector_name, kind, cs)

    def get_tx_cellstate(self, name, state):
        kind = 'tx_cellstate'
        cs = {}
        cs['type'] = kind
        cs['name'] = name
        cs['state'] = state
        selector_name = '[%s].[%s].[%s]' % (kind, name, state)

        return Selector(selector_name, kind, cs)

    def get_cells(self, name, compartment=None):
        kind = 'cell'
        cs = {}
        cs['type'] = kind
        cs['name'] = name
        if compartment:
            cs['compartment'] = compartment

        selector_name = '[%s].[%s].[%s]' % (kind, name, compartment)

        return Selector(selector_name, kind, cs)

    def get_cytokine(self, name, compartment=None):
        kind = 'cytokine'
        cs = {}
        cs['type'] = kind
        cs['name'] = name
        if compartment:
            cs['compartment'] = compartment

        selector_name = '[%s].[%s].[%s]' % (kind, name, compartment)

        return Selector(selector_name, kind, cs)

    # SPECIES Getters - get something from one of self.x based on name and type
    def get_species(self, type, name):
        map = {}
        map['tx_cell'] = self.tx_cells
        map['cell'] = self.cells
        map['cytokine'] = self.cytokines

        source = map[type]

        for x in source:
            if x['name'] == name:
                return x
        return None

    # Linkage Specifiers
    def constant(self, name, value):
        return Constant(name, value)

    def add_tx_state_linkage(self, tx_cell_name, a, b, func):
        target = self.get_species('tx_cell', tx_cell_name)
        # stateLink:
        # 'a' : selector of type 'tx_cellstate'
        # 'b' : selector of type 'tx_cellstate'
        # 'c' : sympy expression of constants and selectors

        stateLink = {}
        stateLink['a'] = a
        stateLink['b'] = b
        stateLink['func'] = func

        if stateLink in target['state_linkages']:
            warn('tried to add duplicate tx_state_linkage.')
        else:
            target['state_linkages'].append(stateLink)

    def add_tx_cell_killtarget(self, tx_cell_name, target_cells, killer_states):
        target = self.get_species('tx_cell', tx_cell_name)
        # killLink:
        # 'target' : selector of type 'cells'
        # 'killer_states' = array of selectors of type 'tx_cellstate'

        killLink = {}
        killLink['target'] = target_cells
        killLink['killer_states'] = killer_states

        target['killing_linkages'].append(killLink)

    def set_tx_cell_daughter_state(self, tx_cell_name, daughter_state):
        target = self.get_species('tx_cell', tx_cell_name)
        # daughter_state : selector of type 'tx_cellstate'
        target['daughter_state'] = daughter_state

    def add_tx_cytokine_interaction(self, tx_cell_name, cytokine, states, action):
        target = self.get_species('tx_cell', tx_cell_name)

        # cytokineLink:
        # 'target_cytokine' : selector of type 'cytokine'
        # 'states' : array of selectors of type 'tx_cellstate'
        # 'action' : either 'sink' or 'secrete'

        cytokineLink = {}
        cytokineLink['target_cytokine'] = cytokine
        cytokineLink['states'] = states
        cytokineLink['action'] = action

        target['cytokine_linkages'].append(cytokineLink)

    # HELPERS
    def get_adjacent_compartments(self, compartment):
        output = []
        for linkage in self.compartment_linkages:
            if compartment in linkage:
                target = None
                if linkage[0] == compartment:
                    target = linkage[1]
                elif linkage[1] == compartment:
                    target = linkage[0]
                if target not in output:
                    output.append(target)
        return output

    def gen_states_for_tx_cell(self, tx_cell):
        state_names = tx_cell['states']
        # vars is an array of string names
        # desired output is:
        # [[(activated, 0), (primed, 0)], [(activated, 1), (primed, 1)], ...]
        # for each state, iterate through the other states
        n = len(state_names)
        out = []
        # Iterate through all possible states of n binary vars
        for i in range(1 << n):
            s = bin(i)[2:]
            s = '0' * (n - len(s)) + s
            state = list(map(int, list(s)))
            state_description = [(state_names[j], state[j]) for j in range(n)]
            out.append(state_description)
        return out

    def convert_bio_sel_to_sys(self, sys, selector, compartment_override=None):
        # map every biolayer selector type onto a syslayer relationship
        cs = selector.selector
        biolayer_sel_type = cs['type']

        if compartment_override != None:
            cs['compartment'] = compartment_override

        # biolayer selectors types: 'tx_celltype', 'tx_cellstate', 'cell', 'cytokine'
        # syslayer selector types: 'element' and 'element_state'
        # both have 'type', 'target_type', 'target_name', 'target_compartment'
        # 'element_state' has 'target_state' specifier

        nonstate_types = ['tx_cellstate', 'cell', 'cytokine']

        if biolayer_sel_type in nonstate_types:
            # get the corresponding element selector from the syslayer
            sel = sys.get_element(kind=biolayer_sel_type, name=cs['name'], compartment=cs['compartment'])
            return sel

        elif biolayer_sel_type == 'tx_cellstate':
            # get the corresponding elements_state selector from the syslayer
            sel = sys.get_element_state(kind=biolayer_sel_type, name=cs['name'], compartment=cs['compartment'], state=cs['state'])
            return sel
        else:
            warn('BioLayer was was unable to convert selector type %s to sys' % biolayer_sel_type)

    def convert_bio_func_to_sys(self, sys, func):
        # map a biolayer expression (func) onto a syslayer compatible one
        # i.e., convert each selector within the expression to a syslayer compatible one
        # and leave the constants as is
        # ALGORITHM
        # For each term in expr.args:
        # if it's a constant, leave it
        # if it's a selector, convert it via above method and substitute
        # if it's an expression, pass rescursively to this function and substitute

        # TODO: I have no clue if this recursion is going to work.

        for arg in func.args:
            if isinstance(arg, Selector):
                # term is a biolayer selector, convert to syslayer selector
                sys_selector = self.convert_bio_sel_to_sys(sys, arg)
                func.subs(arg, sys_selector)
            elif isinstance(arg, Constant):
                # term is a constant, don't change it
                sys_constant = arg
                func.subs(arg, sys_constant)
            else:
                # term is a composite, recurse
                result = self.convert_bio_func_to_sys(sys, arg)
                func.subs(arg, result)

        return func

    # CORE
    def compose(self):
        # Generate a system layer from this biology layer
        sys = syslayer.SysLayer()

        # Passthrough the compartments
        sys.compartments = self.compartments
        sys.compartment_linkages = self.compartment_linkages

        # CREATE ALL THE ELEMENTS
        # Elements exist for each species for each compartment.
        # Create elements corresponding to the tx cells in each compartment
        for tx in self.tx_cells:
            for compartment in self.compartments:
                sys.add_element(kind='tx_cell', name=tx['name'], compartment=compartment, states=tx['states'])

        # Create elements corresponding to regular cells
        for cell in self.cells:
            sys.add_element(kind='cell', name=cell['name'], compartment=cell['compartment'])

        # Create elements corresponding to cytokines in each compartment
        for cytokine in self.cytokines:
            for compartment in self.compartments:
                sys.add_element(kind='cytokine', name=cytokine['name'], compartment=compartment)

        # Add the Tx cell linkages
        # * migration (inter-compartment, full-autogen)
        # * death (intra-compartment, full-autogen)
        # * proliferation (intra-compartment, semi-autogen: daughter_cellstate)
        # * state changes (intra-compartment, semi-autogen: state_linkages)
        # * killing (intra-compartmental, semi-autogen: killer_state, kill_target)
        # * cytokine modulation (intra-compartment, semi-autogen: target, secretion_state, action)

        for tx_cell in self.tx_cells:
            # migration - for each compartment, join each state to equiv in each adjacent compartment
            for compartment in self.compartments:
                adjacent = self.get_adjacent_compartments(compartment)
                for state in self.gen_states_for_tx_cell(tx_cell):
                    a = sys.get_element_state('tx_cell', tx_cell['name'], compartment, state)
                    for adj_compartment in adjacent:
                        b = sys.get_element_state('tx_cell', tx_cell['name'], adj_compartment, state)
                        k_name = 'k_mig_%s_to_%s' % (compartment, adj_compartment)
                        func = Constant(k_name, 10)*a
                        sys.add_relationship('migration', a, b, func)

                        # TODO: This algorithm will probably add some redundant linkages.

            # death - for each compartment, join each state to self with death coefficient
            for compartment in self.compartments:
                for state in self.gen_states_for_tx_cell(tx_cell):
                    a = sys.get_element_state('tx_cell', tx_cell['name'], compartment, state)
                    sys.add_relationship('death', a, a, -Constant('k_death',5)*a)

            # proliferation - for each compartment, join each state to daughter_state with birth coefficient
            for compartment in self.compartments:
                for state in self.gen_states_for_tx_cell(tx_cell):
                    a = sys.get_element_state('tx_cell', tx_cell['name'], compartment, state)
                    b = sys.get_element_state('tx_cell', tx_cell['name'], compartment, tx_cell['daughter_state'])
                    func = a*Constant('k_proliferation', 10)
                    sys.add_relationship('prolif', a, b, func)

            # state changes - for each compartment, join each state to relevant other states based on state_linkages
            for compartment in self.compartments:
                for state in self.gen_states_for_tx_cell(tx_cell):
                    a = sys.get_element_state('tx_cell', tx_cell['name'], compartment, state)
                    for linkage in tx_cell['state_linkages']:
                        # if the linkage stems from a, add the linkage
                        if linkage['a'] == a:
                            b = self.convert_bio_sel_to_sys(sys, linkage['b'])
                            func = self.convert_bio_func_to_sys(sys, linkage['func'])
                            sys.add_relationship('state_link', a, b, func)

            # killing - for each compartment, join the killer state to the killed state for the tx cell
            for compartment in self.compartments:
                for kill_link in tx_cell['killing_linkages']:
                    # Get the compartmental target
                    b = self.convert_bio_sel_to_sys(sys, kill_link['target'], compartment_override=compartment)
                    for killer_state in kill_link['killer_states']:
                        a = sys.get_element_state('tx_cellstate', tx_cell['name'], compartment, killer_state)
                        func = a*b*Constant('k_kill', 10)
                        sys.add_relationship('killing', a, b, func)

            # cytokine modulation - for each compartment, join the action state to the cytokine
            for compartment in self.compartments:
                for cytokine_link in tx_cell['cytokine_linkages']:
                    # Get the compartmental cytokine
                    # TODO: THIS IS WRONG THE LINKAGE NEEDS TO BE A -> A NOT TX -> CYTOKINE
                    b = self.convert_bio_sel_to_sys(sys, cytokine_link['target_cytokine'], compartment_override=compartment)
                    for action_state in cytokine_link['states']:
                            a = sys.get_element_state('tx_cellstate', tx_cell['name'], compartment=compartment, state=action_state)
                            func = 1
                            if cytokine_link['action'] == 'secrete':
                                func = a*Constant('k_secrete', 10)
                            elif cytokine_link['action'] == 'sink':
                                func = -a*b*constant('k_sink', 10)
                            sys.add_relationship('cytokine_modulation', a, b, func)

        # Add the cell linkages
        # * proliferation (intra-compartment, full-autogen)
        # * death (intra-compartment, full-autogen)

        for cell in self.cells:
            a = sys.get_element('cell', cell['name'], cell['compartment'])
            pro_func = a*Constant('k_proliferation', 10)
            sys.add_relationship('proliferation', a, a, pro_func)

            death_func = -a*self.constant('k_death', 5)
            sys.add_relationship('death', a, a, death_func)

        # Add the cytokine linkages
        # * diffusion (inter-compartment, full-autogen)
        # * degradation (intra-compartment, full-autogen)

        for cytokine in self.cytokines:
            # for each compartment, add migration linkage to adjacent compartments
            for compartment in self.compartments:
                a = sys.get_element('cytokine', cytokine['name'], compartment)
                for adj in self.get_adjacent_compartments(compartment):
                    b = sys.get_element('cytokine', cytokine['name'], adj)
                    func = a*Constant('k_diffuse', 10)
                    sys.add_relationship('diffusion', a, b, func)

            for compartment in self.compartments:
                a = sys.get_element('cytokine', cytokine['name'], compartment)
                deg_func = -Constant('k_degrade', 10)*a
                sys.add_relationship('degradation', a, a, deg_func)

        return sys
