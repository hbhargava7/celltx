# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell A. Lim
# University of California, San Francisco

"""Design, construct, and simulate multiscale models of cell therapies."""

from warnings import warn
import networkx as nx
from . import biolayer
from . import syslayer
from . import graphlayer
from . import odelayer


class Model:
    """
    Class which holds abstraction layers that allow model definition, ODE system extraction, and simulation.

    Attributes
    ----------
    name : str
        Name of the model
    bio : BioLayer
        BioLayer instance used for high-level model specification
    sys : SysLayer
        SysLayer instance used for mid-level model specification. Generated by a BioLayer, or manually specified.
    graph : GraphLayer
        GraphLayer instance used for low-level model specification. Generated by a SysLayer, or manually specified.
    ode : ODELayer
        ODELayer instance containing a system of differential equations and logic for their modification and evaluation.
    """

    def __init__(self, name):
        """
        Parameters
        ----------
        name : str
            Name of the model for reference purposes.
        """
        self.name = name

        self.bio = biolayer.BioLayer()
        self.sys = syslayer.SysLayer()
        self.graph = graphlayer.GraphLayer()
        self.ode = None

        print("Model \"%s\" Instantiated." % self.name)

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
