# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell A. Lim
# University of California, San Francisco

"""Specify the architecture of a biological system at molecular, cellular, and tissue-level scales."""

from warnings import warn

from .. import syslayer
from ..functions import Selector, Constant


class BioLayer:
    """
    Class which enables specification and storage of multiscale models of cell therapies in a biologically intuitive fashion.

    Helper Pseudo-Classes
    ---------------------
    celltx currently uses a variety of pseudo-classes which are currently implemented as ``dict``s. Their specifications
    are provided here.

    StateLink : dict
        Dictionary describing a linkage between two states of a tx_celltype (e.g. [(primed, 0)] and [(primed, 1)].
        StateLink has the following structure:
            * ``'a'`` : :class:`Selector` of type ``tx_cellstate`` addressing the origin state.
            * ``'b'`` : :class:`Selector` of type ``tx_cellstate`` addressing the destination state.
            * ``'func'`` : sympy.core.expr.Expr expression of constants and selectors describing the linkage.

    CytokineLink : dict
        Dictionary describing a linkage between a tx_cell and a cytokine (e.g. for secretion or sinking)
        CytokineLink has the following structure:
            * ``'target_cytokine'`` : :class:`Selector` of type ``cytokine`` addressing the cytokine being linked.
            * ``'states'`` : list[:class:`Selector`] of type ``tx_cellstate`` addressing the tx_cellstates being linked.
            * ``'action'`` : str either 'sink' or 'secrete' describing the action of the linkage.

    KillLink : dict
        Dictionary describing a linkage between some tx_cellstates and a cell killed by those tx_cellstates.
        KillLink has the following structure:
            * ``'target'`` : :class:`Selector` of type 'cells'
            * ``'killer_states'`` list[:class:`Selector`] of type 'tx_cellstate'

    BioLink : dict
        Dictionary describing a generic, custom linkage between two :class:`Selector` objects.
        BioLink has the following structure:
            * ``'a'`` : :class:`Selector` addressing the origin entity.
            * ``'b'`` : :class:`Selector` addressing the destination entity.
            * ``'func'`` : sympy.core.expr.Expr in terms of constants and selectors.

    Attributes
    __________
    compartments : list[str]
        List of unique compartment names. Linkages specified by :class:`.compartment_linkages`.
    compartment_linkages : list[tuple[str]]
        List of 2-tuples of strings each describing a linkage between two compartments in :class:`.compartments`
    tx_cells : list[dict]
        List containing specifications of each therapeutic celltype in the model (e.g. CAR-T cells).
        Each dictionary contains the following keys and structure:
            * ``'name'`` (str) unique name.
            * ``'states'`` (list[str]) list of names of binary states (e.g. 'primed', 'activated')
            * ``'state_linkages'`` (list[:class:`StateLink`]) list describing linkages between states.
            * ``'daughter_state'`` (list[tuple[str]]) list of 2-tuples addressing the cellstate yielded by proliferation.
            * ``'cytokine_linkages'`` (list[:class:`CytokineLink`]) list describing linkages with cytokines.
            * ``'killing_linkages'`` (list[:class:`KillLink`]) list describing what the tx cell can kill.
    cells : list[dict]
        List describing the ordinary celltypes as ``dict`` objects with structure:
            * ``'name'`` : Unique name
            * ``'compartment'`` : Name of the compartment in which the cells reside.
    cytokines : list[str]
        List describing the cytokines in the model by unique names.
    interactions : list[:class:`SingleSelector`]
        List of custom interactions to add to the model.
    """

    def __init__(self):
        self.compartments = []
        self.compartment_linkages = []
        self.tx_cells = []
        self.cells = []
        self.cytokines = []
        self.interactions = []

    def add_compartment(self, name):
        """Add a compartment to the model."""
        self.compartments.append(name)

    def link_compartments(self, a, b):
        """Add a compartment linkage to the model."""
        self.compartment_linkages.append((a, b))

    def add_tx_cells(self, name, states, state_linkages=[], daughter_state=[], cytokine_linkages=[], killing_linkages=[]):
        """Add a type of therapeutic cells to the model."""
        c = {}
        c['name'] = name
        c['states'] = states
        c['state_linkages'] = state_linkages
        c['daughter_state'] = daughter_state
        c['cytokine_linkages'] = cytokine_linkages
        c['killing_linkages'] = killing_linkages
        self.tx_cells.append(c)

    def add_cells(self, name, compartment, growth_type='Logistic'):
        """Add a type of non-therapeutic cell to the model (e.g. tumor or normal cells)."""
        c = {}
        c['name'] = name
        c['compartment'] = compartment
        c['growth_type'] = growth_type
        self.cells.append(c)

    def add_cytokine(self, name):
        """Add a type of cytokine to the model."""
        c = {}
        c['name'] = name
        self.cytokines.append(c)

    def add_interaction(self, a, b, func):
        """Add a custom interaction to the model."""
        i = {}
        i['a'] = a
        i['b'] = b
        i['func'] = func
        self.interactions.append(i)

    def validate_selector(self, selector):
        """Validate that a selector actually addresses something in the BioLayer."""
        # TODO: Implement selector validation
        return True

    # SELECTOR Getters - get a selector for something
    def get_tx_celltype(self, name):
        """
        Get a :class:`Selector` addressing a therapeutic celltype.

        Parameters
        ----------
        name : str
            Name of the tx_celltype.

        Returns
        -------
        :class:`Selector` of type `tx_celltype`
            Selector addressing the desired therapeutic celltype.

        """
        kind = 'tx_celltype'
        cs = {}
        cs['type'] = kind
        cs['name'] = name

        selector_name = '[%s].[%s]' % (kind, name)
        selector_latex = '\\text{[%s].[%s]}' % (kind, name)
        return Selector(selector_name, kind, cs, latex=selector_latex)

    def get_tx_cellstate(self, name, state):
        """
        Get a :class:`Selector` addressing a therapeutic cellstate.

        Parameters
        ----------
        name : str
            Name of the tx_celltype
        state : list[tuple]
            List of 2-tuples addressing a unique state of the tx_celltype

        Returns
        -------
        :class:`Selector` of type 'tx_cellstate'
            Selector addressing the desired state of the desired therapeutic celltype.
        """

        kind = 'tx_cellstate'
        cs = {}
        cs['type'] = kind
        cs['name'] = name
        cs['state'] = state
        selector_name = '[%s].[%s].[%s]' % (kind, name, state)
        selector_latex = '\\text{[%s].[%s].[%s]}' % (kind, name, state)

        return Selector(selector_name, kind, cs, latex=selector_latex)

    def get_cells(self, name, compartment=None):
        kind = 'cell'
        cs = {}
        cs['type'] = kind
        cs['name'] = name
        if compartment:
            cs['compartment'] = compartment

        selector_name = '[%s].[%s].[%s]' % (kind, name, compartment)
        selector_latex = '\\text{[%s].[%s].[%s]}' % (kind, name, compartment)

        return Selector(selector_name, kind, cs, latex=selector_latex)

    def get_cytokine(self, name, compartment=None):
        kind = 'cytokine'
        cs = {}
        cs['type'] = kind
        cs['name'] = name
        if compartment:
            cs['compartment'] = compartment

        selector_name = '[%s].[%s].[%s]' % (kind, name, compartment)
        selector_latex = '\\text{[%s].[%s].[%s]}' % (kind, name, compartment)

        return Selector(selector_name, kind, cs, latex=selector_latex)

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
        # The states might be passed as a selector, e.g. from get_tx_cellstate.
        # If this is the case, need to extract the arrays.
        comp_states = []
        for state in killer_states:
            if not isinstance(state, list):
                comp_states.append(state.selector['state'])
            else:
                comp_states.append(state)

        killLink = {}
        killLink['target'] = target_cells
        killLink['killer_states'] = comp_states

        target['killing_linkages'].append(killLink)

    def set_tx_cell_daughter_state(self, tx_cell_name, daughter_state):
        target = self.get_species('tx_cell', tx_cell_name)
        # daughter_state : selector of type 'tx_cellstate'

        # The daughter_state might be passed as a selector, e.g. from get_tx_cellstate.
        # If this is the case, need to extract the array.
        if not isinstance(daughter_state, list):
            daughter_state = daughter_state.selector['state']
        target['daughter_state'] = daughter_state

    def add_tx_cytokine_interaction(self, tx_cell_name, cytokine, states, action):
        target = self.get_species('tx_cell', tx_cell_name)

        # cytokineLink:
        # 'target_cytokine' : selector of type 'cytokine'
        # 'states' : array of selectors of type 'tx_cellstate'
        # 'action' : either 'sink' or 'secrete'

        # The states might be passed as a selector, e.g. from get_tx_cellstate.
        # If this is the case, need to extract the arrays.
        comp_states = []
        for state in states:
            if not isinstance(state, list):
                comp_states.append(state.selector['state'])
            else:
                comp_states.append(state)

        cytokineLink = {}
        cytokineLink['target_cytokine'] = cytokine
        cytokineLink['states'] = comp_states
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
        """
        Convert a BioLayer selector into a SysLayer compartible one
        """

        # map every biolayer selector type onto a syslayer relationship
        cs = selector.selector
        biolayer_sel_type = cs['type']

        if compartment_override != None:
            cs['compartment'] = compartment_override

        # biolayer selectors types: 'tx_celltype', 'tx_cellstate', 'cell', 'cytokine'
        # syslayer selector types: 'element' and 'element_state'
        # both have 'type', 'target_type', 'target_name', 'target_compartment'
        # 'element_state' has 'target_state' specifier

        nonstate_types = ['tx_celltype', 'cell', 'cytokine']

        if biolayer_sel_type in nonstate_types:
            # get the corresponding element selector from the syslayer
            sel = sys.get_element(kind=biolayer_sel_type, name=cs['name'], compartment=cs['compartment'])
            return sel

        elif biolayer_sel_type == 'tx_cellstate':
            # get the corresponding elements_state selector from the syslayer
            sel = sys.get_element_state(kind='tx_cell', name=cs['name'], compartment=cs['compartment'], state=cs['state'])
            return sel
        else:
            warn('BioLayer was was unable to convert selector type %s to sys' % biolayer_sel_type)

    def convert_bio_func_to_sys(self, sys, func, compartment_context=None):
        """
        Convert a :class:`BioLayer` level selector into a :class:`SysLayer` level selector.

        Parameters
        ----------
        sys : :class:`SysLayer`
            SysLayer that has been initialized and populated with elements that will be used to validate selectors.
        func : sympy.core.expr.Expr
            Expr in terms of BioLayer Selectors and Constants
        compartment_context : bool
            Override the compartment context while converting the function

        Returns
        -------
        Sympy.core.expr.Expr
            The converted expression

        Algorithm
        ---------
        For each term in func.args:
            * if it's a constant, leave it alone
            * If it's a selector, convert it via :meth:`convert_bio_sel_to_sys` and substitute
            * If it's an expr, recurse and substitute.
        """

        for arg in func.args:
            if isinstance(arg, Selector):
                # term is a biolayer selector, convert to syslayer selector
                sys_selector = None
                if compartment_context is not None:
                    sys_selector = self.convert_bio_sel_to_sys(sys, arg, compartment_override=compartment_context)
                else:
                    sys_selector = self.convert_bio_sel_to_sys(sys, arg)

                func = func.subs(arg, sys_selector)
            elif isinstance(arg, Constant):
                # term is a constant, don't change it
                sys_constant = arg
                func = func.subs(arg, sys_constant)
            else:
                # term is a composite, recurse
                result = None
                if compartment_context is not None:
                    result = self.convert_bio_func_to_sys(sys, arg, compartment_context=compartment_context)
                else:
                    result = self.convert_bio_func_to_sys(sys, arg)
                func = func.subs(arg, result)
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
                    sys.add_relationship('death', a, a, -Constant('k_death', 5)*a)

            # proliferation - for each compartment, join each state to daughter_state with birth coefficient
            for compartment in self.compartments:
                for state in self.gen_states_for_tx_cell(tx_cell):
                    a = sys.get_element_state('tx_cell', tx_cell['name'], compartment, state)
                    b = sys.get_element_state('tx_cell', tx_cell['name'], compartment, tx_cell['daughter_state'])
                    func = a*Constant('k_proliferation', 10)
                    sys.add_relationship('proliferation', a, b, func)

            # state changes - for each compartment, join each state to relevant other states based on state_linkages
            for compartment in self.compartments:
                for state in self.gen_states_for_tx_cell(tx_cell):
                    a = sys.get_element_state('tx_cell', tx_cell['name'], compartment, state)
                    for linkage in tx_cell['state_linkages']:
                        # if the linkage stems from a, add the linkage
                        link_a_sys = self.convert_bio_sel_to_sys(sys, linkage['a'], compartment_override=compartment)
                        if a.selector == link_a_sys.selector:
                            # The func will be in terms of elements in the local compartment.
                            # It is possible that a term will be undefined in the local compartment (e.g. a+b+ cells)
                            b = self.convert_bio_sel_to_sys(sys, linkage['b'], compartment_override=compartment)
                            func = self.convert_bio_func_to_sys(sys, linkage['func'], compartment_context=compartment)
                            sys.add_relationship('state_link', a, b, func)

            # killing - for each compartment, join the killer state to the killed state for the tx cell
            for compartment in self.compartments:
                for kill_link in tx_cell['killing_linkages']:
                    # Get the compartmental target
                    b = self.convert_bio_sel_to_sys(sys, kill_link['target'], compartment_override=compartment)
                    for killer_state in kill_link['killer_states']:
                        a = sys.get_element_state('tx_cell', tx_cell['name'], compartment, killer_state)
                        func = -a*b*Constant('k_kill', 10)
                        sys.add_relationship('killing', a, b, func)

            # cytokine modulation - for each compartment, join the action state to the cytokine
            for compartment in self.compartments:
                for cytokine_link in tx_cell['cytokine_linkages']:
                    # Get the compartmental cytokine
                    # TODO: THIS IS WRONG THE LINKAGE NEEDS TO BE A -> A NOT TX -> CYTOKINE
                    b = self.convert_bio_sel_to_sys(sys, cytokine_link['target_cytokine'], compartment_override=compartment)
                    for action_state in cytokine_link['states']:
                            a = sys.get_element_state('tx_cell', tx_cell['name'], compartment=compartment, state=action_state)
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
            pro_func = None
            if cell['growth_type'] == 'logistic':
                print("celltx BioLayer: Adding a logistic growth term the cells.")
                pro_func = a*Constant('k_cell_proliferation', 10)
            else:
                print("celltx BioLayer: Adding an exponential growth term to the cells.")
                k_cell_prolif = Constant('k_cell_prolif', 10)
                k_cell_carrycap = Constant('k_cell_carrycap', 10)
                pro_func = k_cell_prolif * a * (1-(1/k_cell_carrycap)*a)

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
