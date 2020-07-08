from ..functions import Selector
from .. import graphlayer


class SysLayer:

    def __init__(self):
        self.compartments = []  # unique compartment names
        self.compartment_linkages = []  # tuples of compartment names

        self.elements = []
        # elements represented as dictionaries with:
        # 'type' : kind (tx_cell, cytokine, cell, element)
        # 'name' : name unique within kind
        # 'comparment' : (optional) comparment of membership
        # 'states' : (optional) array of binary states

        self.relationships = []
        # relationships between element singlets as dictionaries with:
        # 'kind' : kind (proliferation, death, migration, circuitry, relationship
        # 'a' : source singlet selector
        # 'b' : destination singlet selector
        # 'func' : sympy expression in terms of singlet selectors describing

    def add_compartment(self, name):
        self.compartments.append(name)

    def add_compartment_linkage(self, a, b):
        self.compartment_linkages.append((a, b))

    def add_element(self, kind, name, compartment=None, states=[]):
        e = {}
        e['type'] = kind
        e['name'] = name
        e['compartment'] = compartment
        e['states'] = states

        self.elements.append(e)

    def add_relationship(self, name, a, b, function):
        r = {}
        r['name'] = name
        r['type'] = kind
        r['a'] = a
        r['b'] = b
        r['func'] = function
        self.relationships.append(r)

    def validate_selector(self, selector):
        return True

    def get_element(self, kind, name, compartment):
        selector_type = 'element'

        cs = {}
        cs['type'] = selector_type
        cs['target_type'] = kind
        cs['target_name'] = name
        cs['target_compartment'] = compartment

        selector_name = '[%s].[%s].[%s].[%s]' % (selector_type, kind, name, compartment)

        selector = Selector(selector_name, selector_type, cs)

        if self.validate_selector(selector):
            return selector
        else:
            # TODO: Raise an exception here
            pass

    def get_element_state(self, kind, name, compartment, state):
        selector_type = 'element_state'

        cs = {}
        cs['type'] = selector_type
        cs['target_type'] = kind
        cs['target_name'] = name
        cs['target_compartment'] = compartment
        cs['target_state'] = state

        selector_name = '[%s].[%s].[%s].[%s].[%s]' % (selector_type, kind, name, compartment, state)

        selector = Selector(selector_name, selector_type, cs)

        if self.validate_selector(selector):
            return selector
        else:
            # TODO: Raise an exception here
            pass

    def compose(self):
        graph = graphlayer.GraphLayer()

        return graph
