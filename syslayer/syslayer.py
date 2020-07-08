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
        # 'kind' : kind (proliferation, death, migration, circuitry,
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

    def get_element(self, kind, name, compartment, state):
        pass

    def get_element_state(self, kind, name, compartment, state):
        pass

    def compose(self):
        graph = graphlayer.GraphLayer()

        return graph
