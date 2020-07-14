# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell Lim
# University of California, San Francisco

from warnings import warn

import networkx as nx
import sympy as sy

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
		selector_latex = '\\text{[%s].[%s].[%s].[%s].[%s]}' % (selector_type, type, name, compartment, state)
		selector = Selector(selector_name, selector_type, cs, latex=selector_latex)

		if self.validate_selector(selector):
			return selector
		else:
			warn('Failed to get node selector')
			pass

	def get_edge(self):
		pass

	def generate_graph(self, inc_labels=True):
		G = nx.MultiDiGraph()

		for node in self.nodes:
			node_sel = self.get_node(node['type'], node['name'], node['compartment'], node['state'])
			G.add_node(node_sel.name, data=node_sel)  # unique node identifier

		for edge in self.edges:
			a = edge['a'].name
			b = edge['b'].name
			if inc_labels:
				G.add_edge(a,b, label=str(edge['func']), func=edge['func'])
			else:
				G.add_edge(a, b)
		return G

	def selfloops_for_node(self, graph, node):
		sl = nx.selfloop_edges(graph, data=True)
		out = []
		for l in sl:
			if l[0] is node:
				out.append(l)
		return out

	def generate_ode_model(self):
		G = self.generate_graph()

		funcs = nx.get_edge_attributes(G, 'func')
		node_sels = nx.get_node_attributes(G, 'data')

		ODEs = []

		# For each node, generate an equation
		for node in G:
			try:
				node_selector = node_sels[node]
			except:
				warn('failed to get node selector for node %s' % node)
				continue
			node_equation = sy.sympify(0)

			# For self loop edges:
			for selfloop in self.selfloops_for_node(G, node):
				node_equation = node_equation + selfloop[2]['func']

			# For edge pointing to the node, add the edge functions to the equation
			for edge in G.in_edges(node):
				# There may be multiple edges connecting the same two nodes (esp for self loops)
				num_edges = G.number_of_edges(edge[0], edge[1])
				for i in range(num_edges):
					node_equation = node_equation + funcs[(edge[0], edge[1], i)]

			# For each edge originating from the node, subtract edge functions if destination node is same type and name
			for edge in G.out_edges(node):
				try:
					origin_node = node_sels[edge[0]]
					destination_node = node_sels[edge[1]]

					if origin_node.selector['target_type'] is destination_node.selector['target_type'] and \
							origin_node.selector['target_name'] == destination_node.selector['target_name']:
						# There may be multiple edges connecting the same two nodes (esp for self loops)
						num_edges = G.number_of_edges(edge[0], edge[1])
						for i in range(num_edges):
							node_equation = node_equation - funcs[(edge[0], edge[1], i)]
				except:
					warn('tried and failed to process an out_edge for node %s' % node)

			equation = sy.Eq(sy.Derivative(node_selector, sy.sympify('t')), node_equation)

			# ODEs[node] = sy.Eq(node_equation)
			ODEs.append(equation)

		return ODEs
