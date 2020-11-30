# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell A. Lim
# University of California, San Francisco

from warnings import warn

import networkx as nx
import sympy as sy

from ..functions import Selector, Constant
from .. import odelayer


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

		selector_name = '[%s].[%s].[%s].[%s]' % (type, name, compartment, state)
		selector_latex = '(\\text{[%s].[%s].[%s].[%s]})' % (type, name, compartment, state)

		if state is None:
			selector_name = '[%s].[%s].[%s]' % (type, name, compartment)
			selector_latex = '(\\text{[%s].[%s].[%s]})' % (type, name, compartment)

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
				G.add_edge(a, b, label=str(edge['func']), func=edge['func'])
			else:
				G.add_edge(a, b)
		return G

	def write_graph_dotfile(self, path, inc_labels=True, display=False):
		graph = self.generate_graph(inc_labels)
		nx.drawing.nx_agraph.write_dot(graph, path)
		if display:
			import subprocess
			cmd = 'dot -Tpng -Gdpi=400 graph.dot > graph.png; open graph.png'
			subprocess.run(cmd, shell=True)

	def selfloops_for_node(self, graph, node):
		sl = nx.selfloop_edges(graph, data=True)
		out = []
		for l in sl:
			if l[0] is node:
				out.append(l)
		return out

	def in_edges_for_node(self, graph, node):
		edges = graph.edges(data=True)
		out = []

		for edge in edges:
			if edge[1] is node:
				out.append(edge)
		return out

	def out_edges_for_node(self, graph, node):
		edges = graph.edges(data=True)
		out = []

		for edge in edges:
			if edge[0] is node:
				out.append(edge)
		return out

	def compose_ode_system(self):
		G = self.generate_graph()

		funcs = nx.get_edge_attributes(G, 'func')
		node_sels = nx.get_node_attributes(G, 'data')

		equations = []

		# For each node, generate an equation
		for node in G:
			try:
				node_selector = node_sels[node]
			except:
				warn('failed to get node selector for node %s' % node)
				continue
			node_equation = sy.sympify(0)
			# print("------------")
			# print("Addressing Node: %s" % node)
			# print("This node has %i in-edges" % len(G.in_edges(node)))
			# print("This node has %i out-edges" % len(G.out_edges(node)))

			# For edge pointing to the node, add the edge functions to the equation
			for edge in self.in_edges_for_node(G, node):
				node_equation = node_equation + edge[2]['func']
				# print("Added Term (IN): %s" % edge[2]['func'])

			# For each edge originating from the node, subtract edge functions if destination node is same type and name
			for edge in self.out_edges_for_node(G, node):
				try:
					if edge[2]['func'].args[0].name == 'k_proliferation' or edge[2]['func'].args[0].name == 'tx_activ_prolif' \
							or edge[2]['func'].args[0].name == 'k_activ_prolif':
						warn("WARNING: Celltx GraphLayer used extremely hacked up protection clause to not subtract proliferation")
						continue
				except:
					pass
				try:
					origin_node = node_sels[edge[0]]
					destination_node = node_sels[edge[1]]

					# If it's a state change (same type, name, compartment), subtract the term.
					if origin_node.selector['target_type'] == destination_node.selector['target_type'] and \
						origin_node.selector['target_name'] == destination_node.selector['target_name'] and \
						origin_node.selector['target_compartment'] == destination_node.selector['target_compartment']:

						if destination_node.selector != origin_node.selector:

							node_equation = node_equation - edge[2]['func']
							# print("Added Term (OUT, statechange rule): %s" % edge[2]['func'])
				except:
					pass
				try:
					origin_node = node_sels[edge[0]]
					destination_node = node_sels[edge[1]]


					# If it's migration (same type, name, state, different compartments), subtract the term.
					if origin_node.selector['target_type'] == destination_node.selector['target_type'] and \
							origin_node.selector['target_name'] == destination_node.selector['target_name'] and \
							origin_node.selector['state'] == destination_node.selector['state'] and \
							origin_node.selector['target_compartment'] != destination_node.selector['target_compartment']:
						if destination_node.selector != origin_node.selector:
							node_equation = node_equation - edge[2]['func']
							# print("Added Term (OUT, migration rule): %s" % edge[2]['func'])
				except:
					pass

			equation = sy.Eq(sy.Derivative(node_selector, sy.sympify('t')), node_equation)

			equations.append(equation)

		ode_layer = odelayer.ODELayer(equations)
		return ode_layer
