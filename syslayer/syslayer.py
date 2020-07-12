# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell Lim
# University of California, San Francisco

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
        # 'kind' : kind (proliferation, death, migration, circuitry, relationship)
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

    def add_relationship(self, kind, a, b, function):
        r = {}
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
    def convert_sys_sel_to_graph(self, graph, selector, state_override=None):
        cs = selector.selector
        sys_sel_type = cs['type'] # 'element' or 'element_state'
        if state_override is not None:
            cs['target_state'] = state_override # store the overriden state

        # syslayer selector types: 'element' and 'element_state'
        # graphlayer selector type: 'node'

        if sys_sel_type == 'element': # selecting an element without a state
           sel = graph.get_node(type=cs['target_type'], name=cs['target_name'], compartment=cs['target_compartment'])
           return sel
        elif sys_sel_type == 'element_state': # selecting a state of an element
            sel = graph.get_node(type=cs['target_type'], name=cs['target_name'], compartment=cs['target_compartment'], state=cs['target_state'])
            return sel
        else:
            warn('SysLayer was unable to convert selector tyep %s to graph' % sys_sel_type)



    def compose(self):
        # Generate a graph layer from this systems layer
        graph = graphlayer.GraphLayer()

        # CREATE ALL THE NODES
        for element in self.elements:
            # If the element has states, iterate through the states
            if 'states' in element:
                if type(element['states']) == type([]):
                    for state in element['states']:
                        graph.add_node(type=element['type'], name=element['name'], compartment=element['compartment'], state=state)
            else:
                graph.add_node(type=element['type'], name=element['name'], compartment=element['compartment'])

        # CREATE ALL THE EDGES
        for relationship in self.relationships:
            # The relationship selectors 'a' and 'b' as well as the non-constants in 'func' need to be converted to node selectors
            a = self.convert_sys_sel_to_graph(self, graph, selector)
            # IF THERE ARE MULTIPLE STATES....


        return graph
