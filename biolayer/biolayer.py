from .. import syslayer

class BioLayer():

	def __init__(self):
		self.compartments = [] # unique compartment names
		self.compartment_linkages = [] # tuples of compartment names

		self.tx_cells = []
		# tx celltypes represented as dictionaries with:
		# 'name' : name
		# 'states' : ['state1', 'state2', ...]

		self.cells = []
		# celltypes represented as dictionaries with:
		# 'name' : name
		# 'compartment' : compartment_name

		self.cytokines = []
		# cytokines represented as dictionaries with:
		# 'name' : name

		self.tx_state_linkages = []
		# tx cell state linkages represented as dictionaries with:
		# 'a' : state selector
		# 'b' : state selector
		# 'func' : function placeholder

		self.interactions = []
		# interactions between 

		# Overview of linkages:
		# Migration linkages - auto added for any tx_celltype
		# Birth and Death linkages - auto added for any tx_cells, cells
		# Degradation linkage - auto added for cytokines

		# State conversion linkages - manually specified for tx cells at biolayer level
		# depend on compartmental levels of other species
		# s1 and s2 are selected by celltype.name and a cellstate
		# function is an expression in terms of tx_cells, cells, or cytokines
		# function is in terms of selectors with no compartment; interpreter
		# will fill in the compartment when generating the interactions
	
	### Adding Elements
	def add_compartment(self, name):
		self.compartments.append(name)

	def link_compartments(self, a, b):
		self.compartment_linkages.append((a,b))

	def add_cells(self, name, compartment):
		c = {}
		c['name'] = name
		c['compartment'] = compartment
		self.cells.append(c)

	def add_tx_cells(self, name, states):
		c = {}
		c['name'] = name
		c['states'] = states
		self.tx_cells.append(c)

	def add_cytokine(self, name):
		c = {}
		c['name'] = name
		self.cytokines.append(c)

	# Getters
	def get_tx_celltype(self, name):
		cs = {}
		cs['type'] = 'tx_celltype'
		cs['name'] = name

	def get_tx_cellstate(self, name, state):
		cs = {}
		cs['type'] = 'tx_cellstate'
		cs['name'] = name
		cs['state'] = state
		return cs
	
	def get_cells(self, name, compartment=None):
		cs = {}
		cs['type'] = 'cells'
		cs['name'] = name
		if compartment:
			cs['compartment'] = compartment
		return cs

	def get_cytokine(self, name, compartment=None):
		cs = {}
		cs['type'] = 'cytokine'
		cs['name'] = name
		if compartment:
			cs['compartment'] = compartment
		return cs

	# Linkage Specifiers
	def add_tx_state_linkage(self, a, b, func):
		csl = {}
		csl['a'] = a
		csl['b'] = b
		csl['func'] = func
		self.tx_state_linkages.append(csl)

	def add_tx_cell_killtarget(tx, target):
		return

	def add_tx_cellprolif(daughter_state):
		return

	def add_tx_cytokine_interaction(tx, target, states, action):
		return

	# Linkage Helpers
	def arithmetic(self, operation, selectors):
		cs = {}
		cs['type'] = 'arithmetic'
		cs['operation'] = operation
		cs['selectors'] = selectors
		return cs

	def func(self, terms):
		# Function of n terms where f = n[0]*n[1]
		f = {}
		f['type'] = 'function'
		f['terms'] = terms 

	def func_kxz(self, k_name, x, z):
		# function f = k_name*x*z
		f = {}
		f['type'] = 'kxz'
		f['x'] = x
		f['z'] = z
		f['k_name'] = k_name
		return f

	def func_kx(self, k_name, x):
		# function f = k_name*x
		f = {}
		f['type'] = 'kx'
		f['x'] = x
		f['k_name'] = k_name
		return f

	# Core
	def compose(self):
		# Generate a system layer from this biology layer
		sys = syslayer.SysLayer()

		### Passthrough the compartments
		sys.compartments = self.compartments
		sys.compartment_linkages = self.compartment_linkages

		### Create elements corresponding to the tx cells
		for tx in self.tx_cells:
			sys.add_element(kind='tx_cell', name=tx['name'], states=tx['states'])

		for cell in self.cells:
			sys.add_element(kind='cell', name=cell['name'], compartment=cell['compartment'], )

		for cytokine in self.cytokines:
			sys.add_element(kind='cytokine', name=cytokine['name'])
		


		### Add linkages
		## tx cell linkages:
		# proliferation
		# death
		# migration
		# state changes
		# killing
		# cytokine synthesis

		## cell linkages:
		# proliferation
		# death

		## cytokines
		# degradation

		return sys

		# convert the specifications into a syslayer

		# first create all of the elements
		# for 


		# then, create all of the relationships














