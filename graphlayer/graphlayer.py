# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell Lim
# University of California, San Francisco

from warnings import warn

import networkx as nx

from ..functions import Selector, Constant


class GraphLayer():

	def __init__(self):
		self.nodes = []
		# nodes represent entities as dictionaries with:
		# 'type' : kind (tx_cell, cytokine, cell, entity)
		# 'name' : name unique within kind
		# 'compartment' : (optional) compartment of which entity is a member
		# 'state' : (optional) array of binary states describing state

		self.edges = []
		# edges represent transfer functions between entities as dictionaries with:
		# 'type' : kind (proliferation, death, migration, circuitry, influence, state_change)
		# 'a' : source entity selector
		# 'b' : destination entity selector
		# 'func' : sympy expression in terms of singlet selectors

	def add_node(self, type, name, compartment, state=None):
		n = {}
		n['type'] = type
		n['name'] = name
		n['compartment'] = compartment
		n['state'] = state
		self.nodes.append(n)

	def add_edge(self, type, a, b, func):
		e = {}
		e['type'] = type
		e['a'] = a
		e['b'] = b
		e['func'] = func
		self.edges.append(e)
	def validate_selector(self, selector):
		# TODO: Implement validation
		return True
	def get_node(self, type, name, compartment, state=None):
		selector_type = 'node'

		cs = {}
		cs['type'] = selector_type
		cs['target_type'] = type
		cs['target_name'] = name
		cs['target_compartment'] = compartment
		if state is not None:
			cs['state'] = state

		selector_name = '[%s].[%s].[%s].[%s].[%s]' % (selector_type, type, name, compartment, state)

		selector = Selector(selector_name, selector_type, cs)

		if self.validate_selector(selector):
			return selector
		else:
			warn('Failed to get node selector')
			pass

	def get_edge(self):
		pass

	def generate_graph(self):
		G = nx.MultiDiGraph()

		for node in self.nodes:
			node_sel = self.get_node(node['type'], node['name'], node['compartment'], node['state'])
			G.add_node(node_sel.name) # unique node identifier

		for edge in self.edges:
			a = edge['a'].name
			b = edge['b'].name
			G.add_edge(a,b)

		return G

