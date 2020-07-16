import numpy as np
import sympy as sy
from scipy.integrate import odeint

from ..functions import Selector, Constant


class ODELayer():

    def __init__(self, equations):
        self.equations = equations
        self.f_model = None
        self.args = None
        self.params = None

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

    def gen_ode_model(self):
        # Prune any term from an equation that is not either a constant or the LHS of another equation.
        pruned_eqs = []
        lhs = [eq.lhs.args[0] for eq in self.equations]
        for eq in self.equations:
            for arg in self.ravel_expression(eq.rhs):
                if arg not in lhs and not isinstance(arg, Constant):
                    eq = eq.subs(arg, 0)
            pruned_eqs.append(eq)
        self.equations = pruned_eqs

        # Process
        rhss = [i.rhs for i in self.equations]

        all_args = []
        for eq in rhss:
            args = self.ravel_expression(eq)
            for a in args:
                all_args.append(a)

        # remove duplicates
        unique_sels = []
        unique_consts = []
        seen_sel = []
        seen_const = []

        for val in all_args:
            if isinstance(val, Selector):
                if val.selector not in seen_sel:
                    seen_sel.append(val.selector)
                    unique_sels.append(val)
            elif isinstance(val, Constant):
                if val.name not in seen_const:
                    seen_const.append(val.name)
                    unique_consts.append(val)
            else:
                print('unable to handle %s' % val)
        self.args = unique_sels
        self.params = unique_consts

        unique_args = unique_sels + unique_consts

        self.f_model = sy.lambdify(unique_args, rhss)

    def model(self, X, t, args):
        # print("Evaluating Model")
        # print("X = %s" % X)
        # print("args = %s" % args)
        # print("t = %s" % t)
        in_vals = np.concatenate((X, args))
        # print("results = %s" % self.f_model(*in_vals))
        out = self.f_model(*in_vals)
        out = [0 if foo < 0 else foo for foo in out]
        return out

    def split_parameter(self, parameter):
        pass

    def link_parameters(self, a, b):
        pass

    def set_initial_value(self, param, value):
        pass

    def integrate(self):
        pass