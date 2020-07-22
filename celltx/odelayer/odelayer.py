import numpy as np
import sympy as sy
from sympy.utilities.lambdify import lambdastr
from scipy.integrate import odeint

from ..functions import Selector, Constant


class ODELayer():

    def __init__(self, equations):
        self.equations = equations
        self.f_model = None
        self.args = None
        self.params = None
        self.lambda_string = None

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

        # We need the rhss in the same order as the same order as the unique_sels to pass to lambdify.
        ordered_rhss = []

        for arg in unique_sels:
            for eq in self.equations:
                if eq.lhs.args[0].selector is arg.selector:
                    ordered_rhss.append(eq.rhs)

        unique_args = unique_sels + unique_consts

        self.f_model = sy.lambdify(unique_args, ordered_rhss)
        self.lambda_string = lambdastr(unique_args, ordered_rhss)

    def model(self, X, t, args):
        # X = [0 if foo < 0 else foo for foo in X]
        in_vals = np.concatenate((X, args))

        out = self.f_model(*in_vals)

        # If the current value of a var is 0, don't let the differential be less than zero
        new_out = []
        for i, val in enumerate(out):
            if X[i] == 0 and val < 0:
                new_out.append(0)
            else:
                new_out.append(val)

        return new_out

    def display_args(self):
        print("MODEL ARGUMENTS (index | name)")
        for i, arg in enumerate(self.args):
            print("%i | %s" % (i, arg))

        print("\nMODEL PARAMETERS (index | name | value)")
        for i, param in enumerate(self.params):
            print("%i | %s | %.2f" % (i, param, param.expr))

    def split_parameter(self, parameter):
        pass

    def link_parameters(self, a, b):
        pass

    def set_initial_value(self, param, value):
        pass

    def integrate(self):
        pass
