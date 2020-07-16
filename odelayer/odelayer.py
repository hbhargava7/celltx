import numpy as np
from scipy.integrate import odeint

from ..functions import Selector, Constant


class ODELayer():

    def __init__(self, equations):
        self.equations = equations

    def ravel_expression(self, expr):
        args = []
        for arg in expr.args:
            if isinstance(arg, Selector) or isinstance(arg, Constant):
                args.append(arg)
            else:
                args = args + self.ravel_expression(arg)
        return args

    def list_parameters(self):
        constants = []
        for eq in self.equations:
            for term in self.ravel_expression(eq.rhs):
                if isinstance(term, Constant):
                    constants.append(term)
        return list(set(constants))
    def model(self, x, t, p):
        # x is an array of n variables
        # t is the time
        # p is a dict with the param values
        pass

    def split_parameter(self, parameter):
        pass

    def link_parameters(self, a, b):
        pass

    def set_initial_value(self, param, value):
        pass

    def integrate(self):
        pass