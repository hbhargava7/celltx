# Copyright 2020 Hersh K. Bhargava
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
		# 'type' : kind (proliferation, death, migration, circuitry, edge
		# 'a' : source entity selector
		# 'b' : destination entity selector
		# 'func' : sympy expression in terms of singlet selectors

	def add_node(self):
		pass

	def add_edge(self):
		pass

	def get_node(self):
		pass

	def get_edge(self):
		pass

	def generate_graph(self):
		pass

