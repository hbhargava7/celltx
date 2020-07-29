# Copyright 2020 Hersh K. Bhargava (https://hershbhargava.com)
# Laboratories of Hana El-Samad and Wendell A. Lim
# University of California, San Francisco

import sympy as sy


class Selector(sy.Symbol):

    def __new__(cls, name, selector_type, selector, latex=None, domain=None, **assumptions):
        var = super().__new__(cls, name, **assumptions)
        var.selector = selector
        var.selector_type = selector_type
        var.latex = latex
        var.domain = domain
        return var

    def __getnewargs__(self):
        return self.name, self.selector_type, self.selector, self.latex, self.domain

    def _latex(self, printer):
        if self.latex is not None:
            return self.latex
        return printer._print_Symbol(self)

    def __getstate__(self):
        state = super().__getstate__()
        state.update(latex=self.latex, domain=self.domain)
        return state


class Constant(sy.Symbol):

    def __new__(cls, name, expr, latex=None, domain=None, **assumptions):
        var = super().__new__(cls, name, **assumptions)
        expr = sy.sympify(expr)
        var.expr = expr
        if latex is None:
            if len(name) > 2:
                if name[:2] == 'k_':
                    var.latex = 'k_{\\text{%s}}' % name[2:]
                else:
                    var.latex = latex
            else:
                var.latex = latex
        else:
            var.latex = latex
        var.domain = domain
        var.selector = name
        return var

    def __getnewargs__(self):
        return self.name, self.expr, self.latex, self.domain

    def _latex(self, printer):
        if self.latex is not None:
            return self.latex
        return printer._print_Symbol(self)

    def __getstate__(self):
        state = super().__getstate__()
        state.update(latex=self.latex, domain=self.domain)
        return state


def ravel_expression(expr):
    args = []
    for arg in expr.args:
        if isinstance(arg, Selector) or isinstance(arg, Constant):
            args.append(arg)
        else:
            args = args + ravel_expression(arg)
    return args


def hill(x, kmin, kmax, x50, n):
    return kmin + (kmax - kmin) / (1 + ((x50) / (x)) ** n)
