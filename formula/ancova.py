import numpy as np
import sympy
from formula import Formula, I
from terms import Factor, is_term, is_factor, Term # Term for docstrings
from utils import factor_codings


class ANCOVA(object):

    """
    Instantiated with a sequence of (expr, [factor]) tuples.
    If there are no factors, entries can be of the form "expr".
    Similarly, if there is just one factor, entries of the
    sequence can be of the form (expr, factor).

    >>> e = Factor('E', ['B', 'M', 'P']) # "B": bachelors, "M":masters, "P":phd
    >>> p = Factor('P', ['M', 'L']) # "M":management, "L":labor
    >>> x = Term('X')
    >>> f = ANCOVA((x,e),(x,p))
    >>> f.formula
    Formula([1, P_M*X, E_B*X, E_M*X, E_P*X])

    The resulting formula depends on the order of the factors
    in the specification (as it does in R).

    >>> f2 = ANCOVA((x,p),(x,e))
    >>> f2.formula
    Formula([1, P_M*X, E_B*X, E_M*X, E_P*X])
    >>> 

    It also depends on the sorted order of the levels
    of the factor (as it does in R).

    >>> e2 = Factor('E', ['P', 'M', 'B'])
    >>> f = ANCOVA((x,p),(x,e2))
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
            # (sympy, [factors]) or sympy or (sympy, factor)
            try:
                expr, factors = tuple(expr_factors)
            except TypeError: # not a sequence
                expr, factors = expr_factors, ()

            if is_factor(factors):
                factors = [factors]
            factors = tuple(factors)
            self.graded_dict.setdefault(expr, {}).setdefault(len(factors), []).append(factors)

        # aliases for the intercept

        for s in ANCOVA.aliases_for_intercept:
            if s in self.graded_dict:
                for k in self.graded_dict[s].keys():
                    self.graded_dict.setdefault(sympy.Number(1), {}).setdefault(k, []).extend(self.graded_dict[s][k])
                del(self.graded_dict[s])

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

        >>> x = Term('x'); f = Factor('f', [2,1,3]) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,f), (x,(f,h)))
        >>> a.sorted_factors
        {'h': Factor('h', [0, 1]), 'f': Factor('f', [1, 2, 3])}
        >>>

        The ordering of the levels of the factors
        changes which columns are produced when
        a factor is coded as a contrast.

        >> x = Term('x'); f = Factor('f', [2,1,3]) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,h), (x,(f,h)))

        In this example, in the "x:f:h" term, "f" is coded
        as a contrast and its "first" level is dropped. This
        is the "first" of the sorted levels of "f".
        

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

        >>> x = Term('x'); f = Factor('f', range(3)) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,f), (x,(f,h)))
        >>> a.codings
        {1: {}, x: {('f',): [('f', 'indicator')], ('f', 'h'): [('f', 'indicator'), ('h', 'contrast')]}}

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

        >>> x = Term('x'); f = Factor('f', [2,1,3]) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,f), (x,(f,h)))
        >>> a.contrasts
        {'I(x):f': Formula([f_1*x, f_3*x, f_2*x]), 'I(1):1': Formula([1]), 'I(x):f:h': Formula([f_2*h_1*x, f_1*h_1*x, f_3*h_1*x])}

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
        """
        The column slices for corresponding contrast matrices.
        See the docstring of `ANCOVA.contrast_matrices`.

        >>> x = Term('x'); f = Factor('f', [2,1,3]) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,f), (x,(f,h)))
        >>> a.slices['I(x):f']
        slice(1, 4, None)

        """
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
        """
        Return the canonical contrast matrices of the ANCOVA.
        The order is determined by the sorted order of 
        numerical terms in the ANCOVA.
        Each numerical term yields several contrasts.

        >>> x = Term('x'); f = Factor('f', [2,1,3]) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,f), (x,(f,h)))

        >>> a.contrast_matrices['I(x):f']

        array([[ 0.,  1.,  0.,  0.,  0.,  0.,  0.],
               [ 0.,  0.,  1.,  0.,  0.,  0.,  0.],
               [ 0.,  0.,  0.,  1.,  0.,  0.,  0.]])


        Note
        ====

        Not all contrasts are estimable depending on the design
        matrix. Hence, when these contrasts are used to compute F-statistics
        the actual "degrees of freedom" of the F-statistic depends on
        the projection of the rows of the contrast matrices onto the
        row space of the design matrix. Consistent contrast matrices
        can be found using `formula.utils.contrast_from_cols_or_rows`
                                                                    
        """
        
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
        Return the canonical formulae, one for each item
        in self.sequence

        >>> x = Term('x'); f = Factor('f', [2,1,3]) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,f), (x,(f,h)))
        >>> a.formulae
        [Formula([1]), Formula([f_1*x, f_3*x, f_2*x]), Formula([f_2*h_1*x, f_1*h_1*x, f_3*h_1*x])]

        """
        if not hasattr(self, '_formulae'):
            self.contrasts
        return self._formulae

    @property
    def Rstr(self):
        """
        Right-hand side for a formula in R. Any "*" in the
        sympy expressions are replaced with ":".

        >>> x = Term('x'); f = Factor('f', [2,1,3]) ; h =Factor('h', range(2))
        >>> a = ANCOVA((x,f), (x,(f,h)))
        >>> a.Rstr
        '1+x:f+x:f:h'
        >>>

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
    def sequence(self):
        """
        A sequence that can be used to construct an equivalent model.

        >>> x = Term('x'); f = Factor('f', range(3)) ; h =Factor('h', range(2))
        >>> i = Factor('i', range(3))
        >>> a = ANCOVA((x,f), (x,(f,h)),(x,(i,h)))
        >>> a.sequence
        [(1, []), (x, (Factor('f', [0, 1, 2]),)), (x, (Factor('f', [0, 1, 2]), Factor('h', [0, 1]))), (x, (Factor('i', [0, 1, 2]), Factor('h', [0, 1])))]

        """
        result = []
        for expr in self.graded_dict:
            for order in self.graded_dict[expr]:
                for factors in self.graded_dict[expr][order]:
                    result.append((expr, factors))
        return result
    
    @property
    def formula(self):
        """
        Create a Formula using R's rules for
        coding factors. 


        >>> x = Term('x'); f = Factor('f', range(3)) ; h =Factor('h', range(2))
        >>> a=ANCOVA((x,f), (x,(f,h)))
        >>> a.formula
        Formula([f_0*x, f_1*h_1*x, f_2*x, 1, f_0*h_1*x, f_1*x, f_2*h_1*x])
        >>>

        """
        terms = []
        for ff in self.formulae:
            terms += list(ff.terms)
        return Formula(terms)

    # Methods

    def delete_terms(self, *terms):
        """
        >>> x = Term('x'); f = Factor('f', range(3)) ; h =Factor('h', range(2))
        >>> a=ANCOVA((x,f), (x,(f,h)))
        >>> a
        ANCOVA((1, ()),(x, ('f',)),(x, ('f', 'h')))
        >>> a.delete_terms((x,f))
        ANCOVA((1, ()),(x, ('f', 'h')))

        """
        result = self.sequence
        for term in ANCOVA(*terms).sequence:
            result.remove(term)
        return ANCOVA(*result)

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
