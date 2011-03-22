import numpy as np
import sympy
from formula import Formula, I
from terms import Factor, is_term, is_factor, Term # Term for docstrings
from utils import factor_codings, simplicial_complex

from scipy.stats import f
from scikits.statsmodels.regression import OLS

class ANCOVA(object):

    # This flag is defined to avoid using isinstance
    _ancova_flag = True

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

    def __init__(self, *expr_factor_tuples, **keywords):
        # set intercept / main_effect behaviour
        
        add_main_effects = False
        if "add_main_effects" in keywords:
            add_main_effects = keywords['add_main_effects']

        add_intercept=True
        if "add_intercept" in keywords:
            add_intercept = keywords['add_intercept']
        
        self.default_contrast='drop_reference'
        if "default_contrast" in keywords:
            self.default_contrast = keywords['default_contrast']
        
        self.graded_dict = {}

        # create a copy of graded_dict
        # removing duplicate tuples of factors in the values
        for expr_factors in expr_factor_tuples:
            # each element of the sequence should have the form
            # (sympy, [factors]) or sympy or (sympy, factor)
            if is_factor(expr_factors):
                expr_factors = (1, expr_factors)

            try:
                expr, factors = tuple(expr_factors)
            except TypeError: # not a sequence
                expr, factors = expr_factors, ()

            if is_factor(factors):
                factors = [factors]
            factors = tuple(factors)
            l = self.graded_dict.setdefault(sympy.sympify(expr), {}).setdefault(len(factors), [])
            # ensure uniqueness
            if factors not in l:
                l.append(factors)

        # aliases for the intercept

        aliases_for_intercept=[1, 'constant', '(Intercept)']
        for s in aliases_for_intercept:
            if s in self.graded_dict:
                for k in self.graded_dict[s].keys():
                    self.graded_dict.setdefault(sympy.Number(1), {}).setdefault(k, []).extend(self.graded_dict[s][k])
                del(self.graded_dict[s])

        if add_intercept:
            self.graded_dict.setdefault(sympy.Number(1), {})[0] = [()]

        if add_main_effects:
            for expr in self.graded_dict:
                self.graded_dict[expr][0] = [()] 

    def __repr__(self):
        return "ANCOVA(%s)" % str(self.sequence())

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
    def contrast_names(self):
        """
        Return the canonical order of the contrasts of the ANCOVA
        formula that will be used, for instance, in type I sum
        of squares.

        The numeric variables are sorted, then within
        each numeric, the graded order of the
        factors is preserved.
        """
        if not hasattr(self, '_contrast_names'):
            self.contrasts
        return self._contrast_names

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
            self._contrast_names = []
            self._formulae = []
            for expr in sorted(self.graded_dict.keys()):
                formulae = get_contributions(self.codings[expr],
                                             self.sorted_factors,
                                             contrast=self.default_contrast)
                for formula, srep in formulae:
                    v = formula * Formula([expr])
                    if str(expr) != '1':
                        sexpr = "I(%s)" % str(expr)
                        if srep and srep != '1':
                            crep = ':'.join([sexpr,srep])
                        else:
                            crep = sexpr
                    elif srep:
                        crep = ':'.join([srep])
                    else:
                        crep = '1'
                    self._contrasts[crep] = v
                    self._contrast_names.append(crep)
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
        for crep in self._contrast_names:
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
        for crep in self._contrast_names:
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
        in self.sequence()

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
        '1+I(x):f+I(x):f:h'
        >>>

        """

        if not hasattr(self, '_contrast_names'):
            self.contrasts
        return '+'.join(self.contrast_names)

    def sequence(self, expr=None):
        """
        A sequence that can be used to construct an equivalent model
        if expr is None. If expr is not None, then this
        model is equivalent to the given model but only
        has numeric terms matching expr.

        >>> x = Term('x'); f = Factor('f', range(3)) ; h =Factor('h', range(2))
        >>> i = Factor('i', range(3))
        >>> a = ANCOVA((x,f), (x,(f,h)),(x,(i,h)))
        >>> a.sequence()
        [(1, []), (x, (Factor('f', [0, 1, 2]),)), (x, (Factor('f', [0, 1, 2]), Factor('h', [0, 1]))), (x, (Factor('i', [0, 1, 2]), Factor('h', [0, 1])))]

        """
        result = []
        if expr is None:
            exprs = sorted(self.graded_dict.keys())
        else:
            exprs = [sympy.sympify(expr)]
        for expr in exprs:
            for order in self.graded_dict[expr]:
                for factors in self.graded_dict[expr][order]:
                    result.append((expr, factors))
        return result
    
    def all_but_above(self, expr, factors):
        """
        Create an ANCOVA for formula
        that has all variables except everything
        in self that is "greater" than (expr, factors).
        That is, all pairings of (expr, factors) where
        factors is not a superset of factors
        as well as everything else in the model.

        This is used in creating type II sums of squares tables.
        """
        sequence = []
        for k in self.graded_dict:
            if k != expr:
                sequence += self.sequence(k)
            else:
                for order in self.graded_dict[k]:
                    for fs in self.graded_dict[k][order]:
                        if not set(fs).issuperset(factors):
                            sequence.append((expr, fs))
        return ANCOVA(*sequence)

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
        if self.formulae:
            terms = []
            for ff in self.formulae:
                terms += list(ff.terms)
            return Formula(terms)
        else:
            return Formula([0])

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
        result = self.sequence()
        for term in ANCOVA(*terms).sequence():
            result.remove(term)
        return ANCOVA(*result)

    def __mul__(self, other):
        if is_factor(other):
            return self.multiply_by_factor(other)
        elif is_ancova(other):
            result = []
            oseq = other.sequence
            sseq = self.sequence

            for se, sf in self.sequence():
                for oe, of in other.sequence():
                    result.append((se*oe, sf+of))
            return ANCOVA(*result)
        else:
            return self.multiply_by_expression(other)

    def __add__(self, other):
        """
        Note: not commutatitive.

        >>> x = Term('x'); f = Factor('f', range(3)); h = Factor('h', range(4))
        >>> a1 = ANCOVA((x,f))
        >>> a2 = ANCOVA((x,h))
        >>> a1 + a2
        >>> a1 + a2
        ANCOVA([(1, ()), (x, (Factor('f', [0, 1, 2]),)), (x, (Factor('h', [0, 1, 2, 3]),))])

        """
        if not is_ancova(other):
            raise ValueError('other should be an ANCOVA formula')
        return concat(self, other)

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

def get_contributions(codings, sorted_factors, contrast='main_effect'):
    """
    Determine which columns a subset of factors
    """
    if codings:
        formulae = []
        lens = []
        for prod_of_factors in codings:
            cur_formula = I
            for n, c in codings[prod_of_factors]:
                if c == 'indicator':
                    cur_formula = cur_formula * sorted_factors[n].indicator
                else:
                    cur_formula = cur_formula * getattr(sorted_factors[n], contrast)
            formulae.append((cur_formula, ':'.join(prod_of_factors)))
            lens.append(len(prod_of_factors))
        formulae = [s[1] for s in sorted(zip(lens, formulae))]
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

def maximal(ancova):
    """
    Return an ANCOVA formula with only the maximal elements
    for each expression. 
    """
    result = []
    for expr in ancova.graded_dict:
        maximal = simplicial_complex(*[s for _, s in ancova.sequence(expr)])[1]
        for m in maximal:
            result.append((expr, m))
    return ANCOVA(*result)

def concat(*ancovas):
    """
    Create a new ANCOVA formula by concatenating a sequence
    of ANCOVA formulae.

    Note: this is not commutatitive because the
    order in which the (expr, [factors]) appear in the
    initiating sequence changes the resulting formula.

    >>> x = Term('x'); f = Factor('f', range(3)); h = Factor('h', range(4))
    >>> a1 = ANCOVA((x,f))
    >>> a2 = ANCOVA((x,h))
    >>> concat(a1,a2).formula
    Formula([h_3*x, f_0*x, h_1*x, f_2*x, 1, f_1*x, h_2*x])
    >>> concat(a2,a1).formula
    Formula([h_3*x, h_1*x, f_2*x, 1, f_1*x, h_2*x, h_0*x])
    >>> 

    """
    result = []
    for ancova in ancovas:
        result += ancova.sequence()
    return ANCOVA(*result)

def is_ancova(obj):
    """ Is obj an ANCOVA?
    """
    return hasattr(obj, "_ancova_flag")


def typeI(response, ancova, recarray):
    """
    Produce an ANCOVA table
    from a given ANCOVA formula
    with type I sums of squares
    where the order is based on the order of terms
    in the contrast_names of ancova.

    Inputs
    ------

    response: str
              field name of response in recarray

    ancova: ANCOVA
            specifies the model to be fit

    recarray: np.ndarray
              should contain all field names in the terms of ancova
              as well as response
    """

    Y = recarray[response]
    X = ancova.formula.design(recarray, return_float=True)
    model = OLS(Y, X)
    results = model.fit()
    SSE_F = np.sum(results.resid**2)
    df_F = results.df_resid

    model = OLS(Y, ancova.formulae[0].design(recarray, return_float=True))
    results = model.fit()
    SSE_old = np.sum(results.resid**2)
    df_old = results.df_resid

    names = []
    sss = []
    fs = []
    dfs = []
    pvals = []

    names.append(ancova.contrast_names[0])
    fs.append(((np.sum(Y**2) - SSE_old) / (Y.shape[0] - df_old)) / (SSE_F / df_F))
    sss.append((np.sum(Y**2) - SSE_old))
    dfs.append(Y.shape[0] - df_old)
    pvals.append(f_dbn.sf(fs[-1], Y.shape[0]-df_old, df_F))

    for d in range(1,len(ancova.formulae)):
        terms = []
        for f in ancova.formulae[:(d+1)]:
            terms += list(f_dbn.terms)

        # JT: this is not numerically efficient
        # could be done by updating some factorization of the full X

        X = Formula(terms).design(recarray, return_float=True)
        model = OLS(Y, X)
        results = model.fit()
        SSE_new = np.sum(results.resid**2)
        df_new = results.df_resid

        sss.append(SSE_old - SSE_new)
        dfs.append(df_old - df_new)
        fs.append(((SSE_old-SSE_new) / (df_old - df_new)) / (SSE_F / df_F))
        pvals.append(f_dbn.sf(fs[-1], df_old-df_new, df_new))
        names.append(ancova.contrast_names[d])
        SSE_old = SSE_new
        df_old = df_new

    result = np.array(names, np.dtype([('contrast','S%d' % max([len(n) for n in names]))]))
    result = ML.rec_append_fields(result, ['SS', 'F', 'df', 'p_value'], [sss, fs, dfs, pvals])
    return result

def typeII(response, ancova, recarray):
    """
    Produce an ANCOVA table
    from a given ANCOVA formula
    with type II sums of squares.

    Inputs
    ------

    response: str
              field name of response in recarray

    ancova: ANCOVA
            specifies the model to be fit

    recarray: np.ndarray
              should contain all field names in the terms of ancova
              as well as response
    """

    Y = recarray[response]
    X = ancova.formula.design(recarray, return_float=True)
    model = OLS(Y, X)
    results = model.fit()
    SSE_F = np.sum(results.resid**2)
    df_F = results.df_resid

    names = []
    sss = []
    Fs = []
    dfs = []
    pvals = []

    for name, expr_factors in zip(ancova.contrast_names,
                                  ancova.sequence()):
        expr, factors = expr_factors
        R = ancova.all_but_above(expr, factors)
        F = ANCOVA(*(R.sequence() + [(expr, factors)]))
        if R.sequence():
            XR = R.formula.design(recarray, return_float=True)
        else:
            XR = np.zeros((Y.shape[0], 1))
        XF = F.formula.design(recarray, return_float=True)
        modelF = OLS(Y, XF)
        modelR = OLS(Y, XR)
        resultsR = modelR.fit()
        resultsF = modelF.fit()

        SSEF = np.sum(resultsF.resid**2)
        dfF = resultsF.df_resid

        SSER = np.sum(resultsR.resid**2)
        dfR = resultsR.df_resid

        sss.append(SSER-SSEF)
        dfs.append(dfR - dfF)
        Fs.append(((SSER-SSEF) / (dfR-dfF)) / (SSE_F / df_F))
        pvals.append(f_dbn.sf(Fs[-1], dfR-dfF, df_F))
        names.append(name)

    result = np.array(names, np.dtype([('contrast','S%d' % max([len(n) for n in names]))]))
    result = ML.rec_append_fields(result, ['SS', 'F', 'df', 'p_value'], [sss, Fs, dfs, pvals])
    return result

        
if __name__ == "__main__":

    import urllib
    import matplotlib.mlab as ML
    from StringIO import StringIO
    from terms import fromrec
    
    s= StringIO()
    s.write(urllib.urlopen('http://stats191.stanford.edu/data/salary.csv').read())
    s.seek(0)
    d = ML.csv2rec(s)
    terms = fromrec(d)

    f = ANCOVA(1, terms['e'],terms['p'],(1,(terms['e'],terms['p'])))
    ANCOVA.add_intercept = False
    f2 = ANCOVA(terms['e'],(1,(terms['e'],terms['p'])))
    ANCOVA.add_intercept = True
    f3 = ANCOVA((1,(terms['e'],terms['p'])))
    
def typeIII(response, ancova, recarray):
    """
    Produce an ANCOVA table
    with type III sum of squares
    from a given ANCOVA formula.

    Inputs
    ------

    response: str
              field name of response in recarray

    ancova: ANCOVA
            specifies the model to be fit

    recarray: np.ndarray
              should contain all field names in the terms of ancova
              as well as response
    """

    X = ancova.formula.design(recarray, return_float=True)
    Y = recarray[response]
    model = OLS(Y, X)

    results = model.fit()
    names = []
    fs = []
    dfs = []
    sss = []
    pvals = []
    for contrast in ancova.contrast_names:
        r = results.f_test(ancova.contrast_matrices[contrast])
        names.append(contrast)
        fs.append(r.fvalue)
        dfs.append(r.df_num)
        pvals.append(r.pvalue)
        sss.append(r.fvalue * results.scale * r.df_num)
    result = np.array(names, np.dtype([('contrast','S%d' % max([len(n) for n in names]))]))
    result = ML.rec_append_fields(result, ['SS', 'F', 'df', 'p_value'], [sss, fs, dfs, pvals])
    return result

## from pkg_resources import resource_stream
## data = scipy.io.loadmat(resource_stream("formula","data/salary.csv"))
