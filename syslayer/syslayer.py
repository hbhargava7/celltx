from .. import graphlayer

class SysLayer():

	def __init__(self):
		self.compartments = []
		self.compartment_linkages = []
		self.elements = []
		self.relationships = []

	def add_element(self, kind, name, compartment=None, states=None):
		e = {}
		e['type'] = kind
		e['name'] = 'name'
		if (compartment):
			e['compartment'] = compartment
		if (state):
			e['states'] = state

		self.elements.append(e)


	def add_compartment(self, name):
		self.compartments.append(name)

	def add_compartment_linkage(self, a, b):
		self.compartment_linkages.append((a,b))

	def add_relationship(self, name, a, b, function):
		r = {}
		r['name'] = name
		r['a'] = a
		r['b'] = b
		r['func'] = function
		self.relationships.append(r)

	def compose(self):
		graph = graphlayer.GraphLayer()

		return graph