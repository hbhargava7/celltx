# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell A. Lim
# University of California, San Francisco

import numpy as np
from numba import njit
import sympy as sy
from sympy.utilities.lambdify import lambdastr
from scipy.integrate import odeint
from smt.sampling_methods import LHS
from tqdm import tqdm
import math, time
import sys
import multiprocessing as mp
from warnings import warn

from ..functions import Selector, Constant
from ..util import format_timedelta


class ODELayer():

    def __init__(self, equations):
        self.equations = equations
        self.f_model = None
        self.species = None
        self.params = None
        self.lambda_string = None
        self.unique_args = None
        self.ordered_rhss = None
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
        self.unique_args = unique_args
        self.ordered_rhss = ordered_rhss
        self.f_model = sy.lambdify(unique_args, ordered_rhss)
        self.f_model = njit(self.f_model)
        print("celltx ODELayer: Successfully compiled self.f_model to C via njit.")
        self.lambda_string = lambdastr(unique_args, ordered_rhss, dummify=True)

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

        # Hacky solution to the quantization problem
        # If the sign of the derivative for a species is negative and the value <1, set it and its derivative to 0.
        # for idx, dx in enumerate(out):
        #     x_val = X[idx]
        #     if x_val < 1 and x_val != 0 and 0 > dx > -1e-100:
        #         out[idx] = -1e200
        #         print('celltx odelayer - quantization rule activated for species at index %i (value is %f). dx = %f. \
        #          Time = %f'%(idx, x_val, dx, t))

        # If the current value of a var is 0, don't let the differential be less than zero
        new_out = []
        for i, val in enumerate(out):
            if X[i] <= 0 and val < 0:
                new_out.append(0)
            else:
                new_out.append(val)
        return new_out

    def index_of_species(self, species_name):
        """
        Get the index within `self.species` of the species named `species_name`.

        Parameters
        ----------
        species_name

        Returns
        -------
        int
        """

        for i, species in enumerate(self.species):
            if species.name == species_name:
                return i
        warn('Celltx ODELayer index_of_species was unable to find species named %s in the model.' % species_name)

    def integrate_quash_species(self, t, species_idx, threshold=1, target=0):
        """
        Goal is to account for the fact that a species in the model should not be able to resurge from a value less than 1.

        Algorithm:
            1. Integrate `self.model` over timespace `t` as normal.
            2. Find (if it happens) the timepoint `t_crit` at which species at `species_idx` falls below `threshold`.
            (after first rising above `threshold` at least once).
            3. Conduct a new simulation from `t_crit` to `t[-1]`, with species at `species_idx` starting at `target` and all
            other species in the model starting from their value in the initial simulation at `t_crit`.
        """
        from copy import deepcopy

        # Get the initial integral
        Xi = self.integrate(t)
        # Iterate through the timecourse of the species of interest.
        has_risen = False # keep track of whether the species has been > threshold yet (don't trigger if x0 was 0)
        t_crit = -1
        for i, val in enumerate(Xi[:, species_idx]):
            if not has_risen:
                if val >= threshold:
                    # Once the species rises above threshold, we are now watching for it to fall back down.
                    has_risen = True
            else:
                # Check if the species has fallen back down
                if val < threshold:
                    t_crit = i
                    # print('found critical value %.2f at timepoint %i' % (val, i))
                    break
        # if t_crit is still -1, the species never fell back down.
        if t_crit == -1:
            return Xi

        # otherwise, setup a new simulation
        # store current x0 values
        x0_archive = deepcopy(self.x0)

        # now assign self.x0 values from t_crit
        for i in range(Xi.shape[1]):
            timecourse = Xi[:,i]
            target_val = timecourse[t_crit]
            self.x0[i] = target_val

            if i == species_idx:
                self.x0[i] = target

        X = self.integrate(t[t_crit:])

        Xi[t_crit:] = X

        self.x0 = x0_archive

        return Xi

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

    def set_initial_value_byname(self, name, val):
        for i, species in enumerate(self.species):
            if species.name == name:
                self.x0[i] = val
                return

        print('celltx.ODELayer: Error finding species %s.' % name)

    def set_param_value(self, name, val):
        for i, param in enumerate(self.params):
            if param.name == name:
                new = param
                new.expr = val
                self.params[i] = new

    def get_param_value(self, name):
        for param in self.params:
            if param.name == name:
                return param.expr
        warn('celltx ODELayer could not find param with name %s' % name)

    def set_search_range(self, param_name, rnge):
        self.search_ranges[param_name] = rnge

    def execute_paramspace_search(self, t, n_samples, parallel, method='LHS'):
        """
        Sample parameter values from parameter-specific ranges specified in self.search_ranges (dict) and simulate.
        Parameters that don't have an entry in self.search_ranges are not to be sampled.
        """
        print('celltx ODELayer: Generating %i samples from the %i dimensional parameter space.' % (n_samples, len(self.params)))
        parameter_sets = self.gen_paramspace_samples(int(n_samples))

        print('celltx ODELayer: Running parallel simulations on %i processors.' % parallel)
        tic = time.time()

        chunked_paramsets = self.chunks(parameter_sets, parallel)

        manager = mp.Manager()
        reservoir = manager.list()

        jobs = []

        for i, chk in enumerate(chunked_paramsets):
            proc = mp.Process(target=self.process_paramset_chunk, args=(chk, reservoir, t, i))
            jobs.append(proc)
            proc.start()

        for process in jobs:
            process.join()
            process.terminate()

        print("celltx ODELayer: Finished all %i simulations in: %s." % (int(n_samples), format_timedelta(time.time() - tic)))
        a = list(reservoir)

        return a

    def process_paramset_chunk(self, chk, out, t, i):
        pbar = tqdm(chk, position=i,file=sys.stdout)
        pbar.set_description('Processor %i Progress' % (i+1))
        for parameter_set in pbar:
            output = []
            try:
                result = odeint(self.model, self.x0, t, args=(parameter_set,))
                output = [parameter_set, result]
            except Exception as e:
                print('ODELayer encountered exception while integrating: %s' % e)
                output = [parameter_set, e]
            out.append(output)


    def chunks(self, lst, nChunks):
        """Divide a list into n lists where all chunks have even size, except for the last one, which is smaller."""
        out = []
        chunkSize = math.floor(len(lst) / nChunks)
        for i in range(0, nChunks):
            start = i * chunkSize
            if i == nChunks - 1:
                out.append(lst[start:])
            else:
                out.append(lst[start:start + chunkSize])
        return out

    def gen_paramspace_samples(self, n_samples):
        """Use Latin Hypercube Sampling to generate samples of the parameter space (looking at self.search_ranges)"""
        ranges = []
        for param in self.params:
            if param.name in self.search_ranges:
                ranges.append(self.search_ranges[param.name])
            else:
                ranges.append([param.expr, param.expr])
        ranges = np.array(ranges)
        s = LHS(xlimits=ranges)
        param_sets = s(n_samples)
        return param_sets

    def display_equations(self, display_func, substitute=False):
        """
        Print out the model equations. Display_func is IPython.display.display.
        If substituting, for each differential equation, substitute the X in dX/dt in the rhs with 'X' for concision.
        """
        if not substitute:
            for eq in self.equations:
                display_func(eq)
        else:
            substitute_eqs = []
            for eq in self.equations:
                X = eq.lhs.args[0]
                X_symbol = sy.sympify('X_self')
                rhs = eq.rhs.subs(X, X_symbol)
                sub_eq = sy.Eq(sy.Derivative(X, sy.sympify('t')), rhs)
                substitute_eqs.append(sub_eq)
            for eq in substitute_eqs:
                display_func(eq)

    def display_args(self):
        print("MODEL SPECIES (index | name | initial value)")
        for i, arg in enumerate(self.species):
            print("%i | %s | %s" % (i, arg, "{:.2e}".format(self.x0[i])))

        print("\nMODEL PARAMETERS (index | name | value)")
        for i, param in enumerate(self.params):
            print("%i | %s | %s" % (i, param, "{:.2e}".format(float(param.expr))))

    def profile_parameter(self, param_name, values, t):
        """
        Profile the model behavior for timeframe t across values of a parameter in `values` array.

        Parameters
        ----------
        param_name : str
            Name of the parameter to investigate

        values : numpy.ndarray[numpy.float64]
            Values of the parameter to investigate

        t : np.ndarray[numpy.float64]
            Timepoint(s) at which to report the state of the model

        Returns
        -------
        np.ndarray[np.ndarray]
            Array of arrays, one for each parameter value in `values`, reporting the model state at times in `t`.

        """
        # Setup the default parameter array
        idx_of_target_param = None
        params = []
        for idx, param in enumerate(self.params):
            params.append(param.expr)
            if param.name == param_name:
                idx_of_target_param = idx

        # For each parameter value
        final = []
        print('celltx ODELayer executing parameter profile simulations for parameter %s.' % param_name)
        for value in tqdm(values):
            try:
                params[idx_of_target_param] = value
                params = [float(val) for val in params]
                result = odeint(self.model, self.x0, t, args=(params,))
                output = [value, result]
                final.append(output)
            except Exception as e:
                # output = [value, e]
                # final.append(output)
                print("celltx ODELayer encountered exception while evaluating model: %s" % e)
                pass

        return final

    def inspect_terms_for_simulation(self, idx, X, t):
        """
        Given the output X of a simulation (e.g. from self.integrate), and the index of an equation in `self.equations`,
        return a dictionary keyed by the `args` of the equation (linearly combined terms) and with values that are
        arrays with the numerical values of the arg for every timepoint represented in X.

        Parameters
        ----------
        idx : int
            the index of the equation whose args to examine
        X : list
            list (1 x n) of lists (m x 1) where n is the number of timepoints and m is the number of equations
            This has the same format as is returned by self.integrate.
        t : numpy.ndarray(float64)
            1 x n array describing the time at each of the timepoints in X. Only really needed if the model takes time
            as an argument.

        Returns
        -------
        list[tuple]
            List of 2-tuples, in which the first value is the arg and the second is a list of values (1 x n)
        """

        equation = self.equations[idx]
        expression = equation.rhs

        output = []
        for arg in expression.args:
            output.append((arg, []))

        for tp in tqdm(range(len(X))):  # for each timepoint
            for j,arg in enumerate(expression.args):
                # substitute each term in the arg with the appropriate value.
                for term in self.ravel_expression(arg):
                    if isinstance(term, Constant):
                        arg = arg.subs(term, term.expr)
                    elif isinstance(term, Selector):
                        # get the value of the selector at the current timepoint
                        index = -1
                        for i, species in enumerate(self.species):
                            if species.selector == term.selector:
                                index = i
                        if index == -1:
                            warn("An internal error occurred; index == -1")
                        value = X[tp, index]
                        arg = arg.subs(term, value)
                    else:
                        warn("Found term in equation that is neither selector nor arg.")
                output[j][1].append(arg)

        return output

