import numpy as np
from scipy.integrate import odeint

from ..functions import Selector, Constant


class ODELayer():

    def __init__(self, equations):
        self.equations = equations
        self.f_model = None

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
        rhss = [i.rhs for i in model.ode.equations]
        lhss = [i.lhs for i in model.ode.equations]

        all_args = []
        for eq in rhss:
            args = model.ode.ravel_expression(eq)
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

        unique_args = unique_sels + unique_consts

        self.f_model = sy.lambdify(unique_args, rhss)

    def model(self, X, t, args):
        in_vals = np.concatenate((X, args))
        return self.f_model(*in_vals)

    def split_parameter(self, parameter):
        pass

    def link_parameters(self, a, b):
        pass

    def set_initial_value(self, param, value):
        pass

    def integrate(self):
        pass