# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell Lim
# University of California, San Francisco

from warnings import warn
import networkx as nx
from . import biolayer
from . import syslayer
from . import graphlayer


class Model:

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
