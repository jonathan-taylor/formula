import numpy as np
import sympy
from formula import Formula, I
from terms import Factor, is_term, is_factor
from utils import factor_codings



class ANCOVA(object):

    """
    Instantiated with a dictionary: dict([(expr, ((factor1, factor2), (factor1,)))]).

    >>> e = Factor('E', ['B', 'M', 'P']) # "B": bachelors, "M":masters, "P":phd
    >>> p = Factor('P', ['M', 'L']) # "M":management, "L":labor
    >>> x = Term('X')
    >>> f = ANCOVA({x:[(e,),(p,)]})
    >>> f.formula
    Formula([1, P_M*X, E_B*X, E_M*X, E_P*X])

    The resulting formula depends on the order of the factors
    in the specification (as it does in R).

    >>> f2 = ANCOVA({x:[(e,),(p,)]})
    >>> f2.formula
    Formula([1, P_M*X, E_B*X, E_M*X, E_P*X])
    >>> 

    It also depends on the sorted order of the levels
    of the factor (as it does in R).

    >>> e2 = Factor('E', ['P', 'M', 'B'])
    >>> f = ANCOVA({x:[(p,),(e2,)]})
    >>> f.formula
    Formula([1, P_M*X, P_L*X, E_M*X, E_P*X])


    """

    add_intercept=True
    add_main_effects=False
    aliases_for_intercept=[1, 'constant', '(Intercept)']

    def __init__(self, *expr_factor_tuples):
        self.graded_dict = {}

        # create a copy of graded_dict
        # removing duplicate tuples of factors in the values
        for expr_factors in expr_factor_tuples:
            # each element of the sequence should have the form
            # (sympy, [factors]) or sympy
            try:
                expr, factors = tuple(expr_factors)
            except TypeError: # not a sequence
                expr, factors = expr_factors, ()

            if is_factor(factors):
                factors = [factors]
            factors = tuple(factors)
            self.graded_dict.setdefault(expr, {}).setdefault(len(factors), []).append(factors)

            ## for sequence_of_factors in graded_dict[expr]:
            ##     if list(sequence_of_factors) not in self.graded_dict[expr]:
            ##         self.graded_dict[expr].append(list(sequence_of_factors))

        # aliases for the intercept

        for s in ANCOVA.aliases_for_intercept:
            if s in self.graded_dict:
                for k in self.graded_dict[s].keys():
                    self.graded_dict.setdefault(sympy.Number(1), {}).setdefault(k, []).extend(self.graded_dict[s][k])
                del(self.graded_dict[s])

        # This is here because (()) == ()
        # so to specify a numeric variable 'x' with no categorical variables
        # we need an entry like {'x':(())} or {'x':[()]} or {'x':[None]}
        ## for expr in self.graded_dict:
        ##     if not self.graded_dict[expr]:
        ##         self.graded_dict[expr] = [()]

        if ANCOVA.add_intercept:
            self.graded_dict.setdefault(sympy.Number(1), {})[0] = [[]]

        if ANCOVA.add_main_effects:
            for expr in self.graded_dict:
                self.graded_dict[expr][0] = [[]] 

    def __repr__(self):
        terms = []
        for expr in self.graded_dict:
            for order in self.graded_dict[expr]:
                for factors in self.graded_dict[expr][order]:
                    terms.append(`(expr, tuple([f.name for f in factors]))`)

        return "ANCOVA(%s)" % ','.join(terms)

    @property
    def sorted_factors(self):
        """
        Take every factor appearing in the formula
        and sort its levels. The sort matters because
        this is the column that is dropped when
        constructing a design in which this factor
        is to be coded as "contrast".
        """
        if not hasattr(self, "_sorted_factors"):
            self._sorted_factors = {}
            for expr in self.graded_dict:
                for order in self.graded_dict[expr]:
                    for factors in self.graded_dict[expr][order]:
                        for factor in factors:
                            if is_factor(factor) and factor.name not in self._sorted_factors:
                                self._sorted_factors[factor.name] = Factor(factor.name, sorted(factor.levels))
        return self._sorted_factors

    @property
    def codings(self):
        """
        Return R's interpretation of how each (expr, [factors])
        instance in the Formula should be coded.
        """
        if not hasattr(self, "_codings"):
            self._codings = {}
            for expr in sorted(self.graded_dict.keys()):
                self._codings[expr] = get_factor_codings(self.graded_dict[expr])
        return self._codings

    @property
    def contrasts(self):
        """
        Return the canonical contrasts of the ANCOVA.
        The order is determined by the sorted order of 
        numerical terms in the ANCOVA.
        Each numerical term yields several 
        """
        if not hasattr(self, "_contrasts"):
            self._contrasts = {}
            self._contrast_reps = []
            self._formulae = []
            for expr in sorted(self.graded_dict.keys()):
                formulae = get_contributions(self.codings[expr],
                                             self.sorted_factors)
                for formula, srep in formulae:
                    v = formula * Formula([expr])
                    sexpr = "I(%s)" % str(expr)
                    if srep:
                        crep = ':'.join([sexpr,srep])
                    else:
                        crep = sexpr
                    self._contrasts[crep] = v
                    self._contrast_reps.append(crep)
                    self._formulae.append(v)
        return self._contrasts

    @property
    def slices(self):
        if not hasattr(self, '_formulae'):
            self.contrasts
        idx = 0
        result = {}
        for crep in self._contrast_reps:
            l = len(self.contrasts[crep].terms)
            result[crep] = slice(idx, idx + l)
            idx += l
        return result

    @property
    def contrast_matrices(self):
        p = len(self.formula.terms)
        matrices = {}
        for crep in self._contrast_reps:
            s = self.slices[crep]
            l = s.stop - s.start
            array = np.zeros((l,p), np.float)
            for i, j in enumerate(range(l)):
                array[i,j+s.start] = 1.
            matrices[crep] = array
        return matrices

    @property
    def formulae(self):
        """
        Return the canonical formulae
        """
        if not hasattr(self, '_formulae'):
            self.contrasts
        return self._formulae

    @property
    def Rstr(self):
        """
        Right-hand side for a formula in R. Any "*" in the
        sympy expressions are replaced with ":".
        """

        results = []
        for expr in self.graded_dict:
            _expr = str(expr).replace("*", ":")
            if _expr != '1':
                for order in self.graded_dict[expr]:
                    for factors in self.graded_dict[expr][order]:
                        results.append(':'.join([_expr] + [f.name for f in factors]))
            else:
                for order in self.graded_dict[expr]:
                    for factors in self.graded_dict[expr][order]:
                        if factors:
                            results.append(':'.join([f.name for f in factors]))
                        else:
                            results.append('1')
        return "+".join(results)

    @property
    def formula(self):
        """
        Create a Formula using R's rules for
        coding factors. 

        If add_intercepts is True,
        then, a "1" is added to the entire expression in R,
        and each sympy expression,s, also appears as 1*s.

        If add_intercepts is False,
        build a Formula strictly from the terms given.
     
        """

        f = self.formulae[0]
        for ff in self.formulae[1:]:
            f = f + ff
        return f.delete_terms(Formula([0]))

    def __mul__(self, other):
        if is_factor(other):
            return self.multiply_by_factor(other)
        else:
            return self.multiply_by_expression(other)

    def multiply_by_factor(self, factor):
        """
        Create a new ANCOVA with each
        existing factor multiplied by factor.
        """
        graded_dict = self.graded_dict.copy()
        for expr in graded_dict:
            result = []
            for factors in graded_dict[expr]:
                result.append(list(factors) + [factor])
            graded_dict[expr] = result
        return ANCOVA(graded_dict)

    def multiply_by_expression(self, expr):
        """
        Create a new ANCOVA with each
        existing expression multiplied by
        expr.
        """
        graded_dict = {}
        for expr in self.graded_dict:
            graded_dict[expr * expr] = self.graded_dict[expr]
        return ANCOVA(graded_dict)


def get_contributions(codings, sorted_factors):
    """
    Determine which columns a subset of factors
    """
    if codings:
        formulae = []
        for prod_of_factors in codings:
            cur_formula = I
            for n, c in codings[prod_of_factors]:
                if c == 'indicator':
                    cur_formula = cur_formula * sorted_factors[n].formula
                else:
                    cur_formula = cur_formula * sorted_factors[n].drop_first
            formulae.append((cur_formula, ':'.join(prod_of_factors)))
    else:
        formulae = [(Formula([1]),'1')]
    return formulae

def get_factor_codings(graded_subsets_of_factors):
    """
    Given a sequence of subsets of factors, determine
    which will be coded with all their degrees of freedom ("indicator")
    and which will be coded as contrasts ("contrast").
    """
    formula = Formula([])
    all_factors = set([])
    graded_subsets_of_names = []

    for order in graded_subsets_of_factors:
        graded_subsets_of_names.extend([sorted([f.name for f in factors]) for
                                          factors in graded_subsets_of_factors[order]])
    if graded_subsets_of_names != [[]]:
        codings = factor_codings(*[sorted(f) for f in graded_subsets_of_names])
    else:
        codings = {}
    return codings