from warnings import warn
import networkx as nx
from . import biolayer
from . import syslayer
from . import graphlayer

class Model():

	### LIFECYCLE
	def __init__(self, name):
		self.name = name
		
		self.bio = biolayer.BioLayer()
		self.sys = syslayer.SysLayer()
		self.graph = graphlayer.GraphLayer()

		print("Model \"%s\" Instantiated." % self.name)

	def compose_biolayer(self):
		self.sys = self.bio.compose()

	def compartment_graph(self):
		G = nx.Graph()
		for compartment in self.compartments:
			G.add_node(compartment)
		for linkage in self.compartment_linkages:
			G.add_edge(linkage[0], linkage[1])
		return G

	def visualize_compartments(self):
		G = self.compartment_graph()
		nx.draw(G, with_labels=True, node_size=2000)

	### COMPARTMENT STATES
	def add_compartment_state(self, compartment, name):
		# WARNING - need to implement redundancy checking here
		self.compartments[compartment].append(name)

	### STATE GETTER
	def cell_state(self, celltype, state):
		# addresses cells on the basis of their circuitry state (e.g. primed unactivated tx cells)
		warn('Have not implemented any state syntax validation.')
		cs = {}
		cs['type'] = 'cell_state'
		cs['celltype'] = celltype
		cs['state'] = state
		return cs

	def compartment_state(self, compartment, state):
		# addresses a compartment-specific variable (e.g. cytokine conc)
		warn('Have not implemented any state syntax validation.')
		cs = {}
		cs['type'] = 'compartment_state'
		cs['compartment'] = compartment
		cs['state'] = state
		return cs

	def ipsi_compartment_state(self, state):
		# tells the interpreter to address the linkage to the version of 
		# the addressed state that is in the same compartment as the other
		# state being linked.
		warn('Have not implemented any state syntax validation.')
		cs = {}
		cs['type'] = 'ipsi_compartment_state'
		cs['state'] = state
		return cs

	### STATE COMPUTATIONS
	def sum_states(self, states):
		cs = {}
		cs['type'] = 'summation_of_states'
		cs['states'] = states

	### STATE LINKAGES
	def link_states(self, a, b, function):
		# I think that a can select multiple entities but b has to select single entities.
		self.state_linkages.append((a,b, function))

	def add_migration_linkages(self):
		for celltype in self.celltypes:
			for state in self.celltypes[celltype]:
				# state looks like ['state1', 'state2', ...]
				
				
				stateName = state[0]
				# stateVal = 

	### EDGE FUNCTIONS
	def edge_func_peg(self, x, name):
		return

	def edge_func_linear(self, x, name):
		return

	def edge_func_hill(self, x, name):
		return

	### CELLS
	def create_celltype(self, name, states):
		if name in self.celltypes:
			raise ValueError("celltype %s already exists" % name)
		self.celltypes[name] = states

	### GRAPH ASSEMBLY
	def gen_circuitry_states(self, state_names):
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
			state = list(map(int,list(s)))
			stateDescription = [(state_names[j], state[j]) for j in range(n)]
			out.append(stateDescription)
		return out

	def gen_model_graph(self):
		# 1. Generate the entity nodes
		self.entities = [] # array of dictionaries
		self.edges = [] # array of dictionaries (keys: 'source' 'destination' 'function')

		# 1.1. Generate the cell state entity nodes
		for celltype in self.celltypes:
			for compartment in self.compartments:
				circuitry_states = self.gen_circuitry_states(self.celltypes[celltype])
				for cstate in circuitry_states:
					entity = {}
					entity['type'] = 'cell_state'
					entity['compartment'] = compartment
					entity['celltype'] = celltype
					entity['circuitry_state'] = cstate
					self.entities.append(entity)

		# 1.2. Generate the compartment state entity nodes
		for compartment in self.compartments:
			for compartment_state in self.compartments[compartment]:
				entity = {}
				entity['type'] = 'compartment_state'
				entity['compartment'] = compartment
				entity['state_name'] = compartment_state
				self.entities.append(entity)

		# 2. Generate the edges
		for linkage in self.state_linkages:
			# need a class reference for linkages
			a = linkage[0]
			b = linkage[1]
			function = linkage[3]

			# a and b are entity selectors.
			# the goal now is to algorithmically find all the entities specified by a selector
			# the current selector types are: 
			# cell_state, compartment_state, and ipsi_compartment_state
			source_entities = self.get_entities_from_selector(a)
			dest_entities = self.get_entities_from_selector(a)

			if len(source_entities) > 1:
				warn('Got multiple source entities for a selector. Taking the first one.')

			# for example, we might be linking  



	def compare_circuitry_states(self, a,b):
		# a and b are arrays of tuples
		# each tuple is ('statename', statevalue)
		warn('Not sure whether compare_circuitry_states implementation is valid.')
		for stateTuple in a:
			if stateTuple not in b: 
				return False
		return True


	def get_entities_from_selector(self, selector, context=''):
		results = []

		if selector['type'] == 'cell_state':
			# cell_state selectors have: type, celltype, state
			# state is an array of tuples like: [('primed',1), ('activated',0)]
			for entity in self.entities:
				if entity['type'] == 'cell_state':
					if entity['celltype'] == selector['celltype']:
						if self.compare_circuitry_states(entity['cstate'], selector['state']):
							results.append(entity)
			return results

		elif selector['type'] == 'compartment_state':
			# compartment_state selectors have: type, compartment, state
			for entity in self.entities:
				if entity['type'] == 'compartment_state':
					if entity['compartment'] == selector['compartment']:
						if entity['state'] == compartment['state']:
							results.append(entity)
			return results

		elif selector['type'] == 'ipsi_compartment_state':
			target_compartment = context['compartment']
			for entity in self.entities:
				if entity['compartment'] == selector['compartment']:
					if entity['state'] == selector['state']:
						results.append(entity)

			return results













