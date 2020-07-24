# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell A. Lim
# University of California, San Francisco

import numpy as np
from numba import njit
import sympy as sy
from sympy.utilities.lambdify import lambdastr
from scipy.integrate import odeint

from ..functions import Selector, Constant


class ODELayer():

    def __init__(self, equations):
        self.equations = equations
        self.f_model = None
        self.species = None
        self.params = None
        self.lambda_string = None
        self.x0 = None
        self.search_ranges = {}

    def ravel_expression(self, expr):
        """
        Given a Sympy expression, return an array of all arguments (selectors and constant) in the expression.

        Parameters
        ----------
        expr : Sympy.core.expr.Expr
            Expression to be unravelled.

        Returns
        -------
        list[Selector or Constant]
        """
        args = []
        for arg in expr.args:
            if isinstance(arg, Selector) or isinstance(arg, Constant):
                args.append(arg)
            else:
                args = args + self.ravel_expression(arg)
        return args

    def gen_ode_model(self):
        """
        Generate a lambda (self.f_model) from the Sympy expressions in self.equations that takes the values for all
        species and parameters and returns dX/dt.

        Returns
        -------
        None
        """

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

        self.species = unique_sels
        self.params = unique_consts

        # We need the rhss in the same order as the same order as the unique_sels to pass to lambdify.
        ordered_rhss = []

        # In the process, we will rearrange self.equations to be in the same order.
        ordered_equations = []

        for arg in unique_sels:
            for eq in self.equations:
                if eq.lhs.args[0].selector is arg.selector:
                    ordered_rhss.append(eq.rhs)
                    ordered_equations.append(eq)

        self.equations = ordered_equations

        unique_args = unique_sels + unique_consts

        self.f_model = sy.lambdify(unique_args, ordered_rhss)
        self.f_model = njit(self.f_model)
        self.lambda_string = lambdastr(unique_args, ordered_rhss)

        # Also generate starting conditions self.x0 (zeros for each species)
        self.x0 = np.zeros(len(self.species))

    def model(self, X, t, args):
        """
        Return the derivative of the system based on current state, desired timepoint, and param values.
        Effectively a wrapper for the lambda self.f_model that, for instance, prevents species from having negative values.

        Parameters
        ----------
        X : list[float]
            List of values for all species in the model, in the same order as self.species.
        t : float
            Time at which the differential is being evaluated.
        args : list[float]
            List of values for all parameters in the model, in the same order as self.params.

        Returns
        -------
        list[float]
            Derivative of the value of each species in the model, in the same order as self.species.
        """

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

    def integrate(self, t):
        """
        Integrate the model at timepoints in t using literal parameter values.

        """
        # Assemble the parameter values into a list.
        params = [float(param.expr) for param in self.params]

        x = odeint(self.model, self.x0, t, args=(params,))
        return x

    def set_initial_value(self, idx, val):
        self.x0[idx] = val

    def set_param_value(self, name, val):
        for i, param in enumerate(self.params):
            if param.name is name:
                new = param
                new.expr = val
                self.params[i] = new

    def set_search_range(self, param, rnge):
        self.search_ranges[param.name] = rnge

    def execute_paramspace_search(self, t, n_samples, method='LHS'):
        """
        Sample parameter values from parameter-specific ranges specified in self.search_ranges (dict) and simulate.
        Parameters that don't have an entry in self.search_ranges are not to be sampled.
        """



    def display_args(self):
        print("MODEL ARGUMENTS (index | name | initial value)")
        for i, arg in enumerate(self.species):
            print("%i | %s | %s" % (i, arg, str(self.x0[i])))

        print("\nMODEL PARAMETERS (index | name | value)")
        for i, param in enumerate(self.params):
            print("%i | %s | %.2f" % (i, param, param.expr))

    def split_parameter(self, parameter):
        pass

    def link_parameters(self, a, b):
        pass


