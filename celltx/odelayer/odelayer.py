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
import copy

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
        self.param_search_ranges = {}
        self.x0_search_ranges = {}
        self.pinned_params = []

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

    def index_of_parameter(self, parameter_name):
        """
        Get the index within `self.parameters` of the parameter named `parameter_name`

        Parameters
        ----------
        parameter_name

        Returns
        -------
        int
        """

        for i, parameter in enumerate(self.params):
            if parameter.name == parameter_name:
                return i
        warn('Celltx ODELayer index_of_parameter was unable to find parameter named %s in the model.' % parameter_name)

    def integrate_quash_species(self, t, species_idxs, start_from=None, override_params=None, override_x0=None, r=1):
        """
        Goal is to account for the fact that species in a model should not be able to resurge from a value less than 1.
        In addition, since species at the `odelayer` abstraction level often actually represent individual states
        of true biological species (e.g. activated and unactivated T cells), it may be important to consider multiple
        species when attempting to quash a biological species.

        Recursive algorithm:
            1. Integrate `self.model` over timespace `t`.
            2. Find (if it exists) the first timeopint `t_crit` at which all species in a group S (from species_idxs)
                have fallen below 1.
            3. Recurse, but start the simulation at `t_crit` and feed a modified x0 where species in S are set to 0.
            4. Return the result if there is no t_crit.

        Parameters
        ----------
        t : np.linspace
            The timeframe to integrate.
        species_idxs : list[tuple]
            List of tuples, where each tuple defines a group of species (addressed by indices in self.species) that
            should be quashed together.
        start_from : np.ndarray or None
            Numpy array with a row for each species and a column for each timepoint. If it is provided, the simulation
            will start at timepoint t[start_from.shape[0]] and the first start_from.shape[0] cols will be pasted in.

        Returns
        -------
        np.ndarray with a row for each species in the model and a column for each timepoint (in `t`).

        """
        # print('multiquash called on species_idxs: %s'%species_idxs)
        # Make it work with the old input format for single species quash
        if not isinstance(species_idxs, list):
            species_idxs = [species_idxs]

        # Warn about not fully tested functionality
        # max_len = max([len(s) for s in species_idxs])
        # if max_len > 1:
        #     warn('celltx.odelayer.multiquash: Warning - the quash_species functionality has only been tested for '
        #          'individual species; not multiple in the same group.')

        if start_from is not None:
            # print('celltx.multiquash starting from timepoint %i with r index %i.' % (start_from.shape[0], r))
            calculated_block = self.integrate(t[start_from.shape[0]:], override_x0=start_from[-1, :],
                                              override_params=override_params)
            X_initial = np.concatenate((start_from, calculated_block), axis=0)
        else:
            # print('celltx.multiquash starting with a fresh integral; start_from was %s. rec = %i'%(start_from, r))
            X_initial = self.integrate(t, override_params=override_params,
                                       override_x0=override_x0)  # initial integration

        d = [False for _ in
             species_idxs]  # list to track the t_crit value for each species group (where it fell below 1)

        # If there are no species left to quash, return
        if len(species_idxs) == 0:
            # print('celltx multiquash returning trivial integrand since species_idxs was empty')
            return X_initial

        # For each species_group, find the first timepoint t_crit where all subspecies were < 1 and store in d.
        # If no such timepoint exists, store False in d.
        # print('searching with species_idxs: %s'%species_idxs)
        for i, species_group in enumerate(species_idxs):
            for timepoint, _X in enumerate(X_initial):  # for each slice of the starting integration
                t_crit = timepoint
                for species in species_group:  # if any of the subspecies are >=1, this ain't it.
                    if _X[species] >= 1:
                        # print('rejecting timepoint %i for species %i bc value %.2f'%(timepoint, species, _X[species]))
                        t_crit = False

                if not isinstance(t_crit, type(False)):  # If t_crit is a number, we found our instance.
                    # print('found t_crit %.2f for species group idx %i. '%(t_crit, i))
                    break

            d[i] = t_crit  # Store the value in d
            # print('found d matrix %s' % d)

        # If all values in d are False, return.
        if sum(d) == 0:
            all_booleans = True
            for val in d:
                if not isinstance(val, type(False)):
                    all_booleans = False
                if all_booleans:
                    # print('celltx.odelayer.multiquash found nothing to quash!')
                    # print('multiquash returning since all values in d are False %s'%d)
                    return X_initial

        # If all values are not False, find the first t_crit and corresponding species.
        dd = [val if not isinstance(val, type(False)) else max(d) + 1 for val in
              d]  # convert False values to large values
        idx_to_quash = dd.index(min(dd))  # index of the species group in species_idxs
        timepoint_to_quash = dd[idx_to_quash]  # timepoint where quashing starts
        species_idxs_to_quash = species_idxs[idx_to_quash]  # tuple of subspecies to actually quash.

        # Now, generate the start_from matrix for the recursion.
        _done = X_initial[0:timepoint_to_quash, :]

        # Set the value of each species in the group we're quashing to 0
        for s in species_idxs_to_quash:
            _done[timepoint_to_quash - 1, s] = 0

        # Generate a species_idxs list that is missing the one we just quashed so that we don't get an infinite loop.
        species_idxs.pop(idx_to_quash)

        # Recurse
        # print('celltx.odelayer.multiquash recursing')
        # print('about to recurse (r=%i)'%r)
        return self.integrate_quash_species(t, species_idxs, start_from=_done, override_params=override_params, r=r+1)

    def integrate(self, t, override_x0=None, override_params=None):
        """
        Integrate the model at timepoints in t using literal parameter values.
        """
        # Assemble the parameter values into a list.
        params = [float(param.expr) for param in self.params]

        if override_params is not None:
            params = override_params

        _x0 = self.x0

        if override_x0 is not None:
            _x0 = override_x0

        x = odeint(self.model, _x0, t, args=(params,))
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

    def set_param_search_range(self, param_name, rnge):
        self.param_search_ranges[param_name] = rnge

    def set_x0_search_range(self, x0_idx, rnge):
        self.x0_search_ranges[x0_idx] = rnge

    def execute_argspace_search(self, t, n_samples, parallel, method='LHS', quash_species=None):
        """
        Explore the space of parameter values and initial conditions specified via `self.set_param_search_range` and
        `self.set_x0_search_range`. First, generate `n_samples` sets of input values using `method`, and then simulate
        (possibly in n-`parallel`), the model timecourse for each sample.

        Parameters
        ----------
        t : np.ndarray
            The timeframe over which to integrate for each set
        n_samples : int
            How many samples to generate
        parallel : int
            Number of processing cores to use
        method : str
            Sampling method
        quash_species : None or int
            Index of a species to quash during the simulation (same functionality as `self.integrate_quash_species`).

        Returns
        -------
        Returns a list of length n_simulations where each element is a 2-lists where 0 is a dictionary of
        x0 and parameter values and 1 is a numpy array of result timecourses X.
        """

        print('celltx ODELayer: Generating %i samples from the %i dimensional argument space.' % (
        n_samples, len(self.params) + len(self.x0)))
        argspace_samples = self.gen_argspace_samples(int(n_samples))

        print('celltx ODELayer: Running parallel simulations on %i processors.' % parallel)
        tic = time.time()

        chunked_argspace_samples = self.chunks(argspace_samples, parallel)

        manager = mp.Manager()
        reservoir = manager.list()

        jobs = []

        for i, chk in enumerate(chunked_argspace_samples):
            proc = mp.Process(target=self.process_argspace_chunk, args=(chk, reservoir, t, i, quash_species))
            jobs.append(proc)
            proc.start()

        for process in jobs:
            process.join()
            process.terminate()

        print('celltx ODELayer: Finished all %i simulations in %s.' % (
        int(n_samples), format_timedelta(time.time() - tic)))
        a = list(reservoir)

        for l in a:

            # Assort the arg vals from the simulation into a dictionary keyed by their names.
            # As per gen_argspace_samples, order is all species followed by all params.

            arg_vals = l[0]
            d = {}
            keys = []
            for species in self.species:
                keys.append(species.name)
            for param in self.params:
                keys.append(param.name)

            for i, pv in enumerate(arg_vals):
                d[keys[i]] = pv
            l[0] = d

        return a

    def process_argspace_chunk(self, chk, out, t, i, quash_species):
        """
        Process a chunk of argspace samples. This function is used for parallelizations.

        Parameters
        ----------
        chk : list
            list of lists, where each list has a value for each species and each parameter (in order of self.x)
        out : multiprocess.reservoir
            reservoir to send the outputs
        t : np.ndarray
            timeframe to integrate
        i : int
            Index of the processor (for reporting realtime progress)
        quash_species : list of lists
            groups of species to quash
        """
        # print('processing chunk with quashed species: %s' % quash_species)

        pbar = tqdm(chk, position=i, file=sys.stdout)
        pbar.set_description('Processor %i Progress' % (i + 1))

        for arg_set in pbar:
            output = []
            try:
                # Divide the arg_set into x0 and params. Order is from self.species and self.params.
                x0 = arg_set[0:len(self.species)]
                params = arg_set[len(self.species):]
                # print('calling iqs with quash_species=%s'%quash_species)
                X = self.integrate_quash_species(t=t, species_idxs=copy.deepcopy(quash_species), override_params=params, override_x0=x0)

                output = [arg_set, X]

            except Exception as e:
                print('ODELayer encountered exception while integrating: %s' % e)
                output = [arg_set, e]
            out.append(output)

    def pin_params_for_search(self, param_names):
        """
        Create pairs of parameters which will have the same value when samples are created, e.g. by `self.gen_argspace_samples`.

        Parameters
        ----------
        params : tuple of parameters to pin (should correspond to elements of `self.params`.

        Returns
        -------
        None
        """

        self.pinned_params.append(param_names)

    def gen_argspace_samples(self, n_samples):
        """
        Use Latin Hypercube Sampling to generate samples of the parameter space (looking at self.param_search_ranges
        and self.x0_search_ranges)
        """
        # Order of ranges will be all species followed by all params.
        ranges = []
        for idx, species in enumerate(self.species):
            ranges.append([self.x0[idx], self.x0[idx]])

        for species_idx in self.x0_search_ranges:
            ranges[species_idx] = self.x0_search_ranges[species_idx]

        for param in self.params:
            if param.name in self.param_search_ranges:
                ranges.append(self.param_search_ranges[param.name])
            else:
                ranges.append([param.expr, param.expr])

        ranges = np.array(ranges)
        s = LHS(xlimits=ranges)
        arg_sets = s(n_samples)

        # now handle the `self.pinned_params`: For each pin tuple, get the indices of the params and make the parameter
        # sample values the same (we'll take the generated samples from the first one.

        # the structure of arg_sets is n lists of a args each.

        for pin_set in self.pinned_params:
            pin_indices = [self.index_of_parameter(p) + len(self.species) for p in pin_set]

            for i, sample in enumerate(arg_sets):
                for j, pin_index in enumerate(pin_indices):
                    # skip over the first instance
                    if j == 0:
                        continue
                    # for subsequent instances, copy the 0th instance value to the jth value.
                    arg_sets[i][pin_index] = sample[pin_indices[0]]

        return arg_sets

    def execute_paramspace_search(self, t, n_samples, parallel, method='LHS', quash_species=None):
        """
        Sample parameter values from parameter-specific ranges specified in self.search_ranges (dict) and simulate.
        Parameters that don't have an entry in self.search_ranges are not to be sampled.

        Parameters
        ----------
        t : np.ndarray
            space of timepoints to profile

        Returns
        -------
        list of 2-lists, where list[0] is parameter values (in same order as self.params), and list[1] is a timecourse
        from self.integrate()
        """
        print('celltx ODELayer: Generating %i samples from the %i dimensional parameter space.' % (
        n_samples, len(self.params)))
        parameter_sets = self.gen_paramspace_samples(int(n_samples))

        print('celltx ODELayer: Running parallel simulations on %i processors.' % parallel)
        tic = time.time()

        chunked_paramsets = self.chunks(parameter_sets, parallel)

        manager = mp.Manager()
        reservoir = manager.list()

        jobs = []

        for i, chk in enumerate(chunked_paramsets):
            proc = mp.Process(target=self.process_paramset_chunk, args=(chk, reservoir, t, i, quash_species))
            jobs.append(proc)
            proc.start()

        for process in jobs:
            process.join()
            process.terminate()

        print("celltx ODELayer: Finished all %i simulations in: %s." % (
        int(n_samples), format_timedelta(time.time() - tic)))
        a = list(reservoir)

        # Comprehend l[0] for l in a into a dictionary keyed by self.params names.

        for l in a:
            param_vals = l[0]
            d = {}
            for i, pv in enumerate(param_vals):
                d[self.params[i].name] = pv
            l[0] = d

        return a

    def process_paramset_chunk(self, chk, out, t, i, quash_species):
        pbar = tqdm(chk, position=i, file=sys.stdout)
        pbar.set_description('Processor %i Progress' % (i + 1))
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
            for j, arg in enumerate(expression.args):
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
