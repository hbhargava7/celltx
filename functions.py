import sympy as sy


class Selector(sy.Symbol):

    def __new__(cls, name, selector_type, selector, latex=None, domain=None, **assumptions):
        var = super().__new__(cls, name, **assumptions)
        var.selector = selector
        var.selector_type = selector_type
        return var


class Constant(sy.Symbol):

    def __new__(cls, name, expr, latex=None, domain=None, **assumptions):
        var = super().__new__(cls, name, **assumptions)
        expr = sy.sympify(expr)
        var.expr = expr
        return var


def hill(x, kmin, kmax, x50, n):
    return kmin + (kmax - kmin) / (1 + ((x50) / (x)) ** n)
