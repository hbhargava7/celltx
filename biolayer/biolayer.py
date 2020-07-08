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
        kind = 'cells'
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
            pass

        # Add the cell linkages
        # * proliferation (intra-compartment, full-autogen)
        # * death (intra-compartment, full-autogen)

        for cell in self.cells:
            pass

        # Add the cytokine linkages
        # * diffusion (inter-compartment, full-autogen)
        # * degradation (intra-compartment, full-autogen)

        for cytokine in self.cytokines:
            pass

        return sys
