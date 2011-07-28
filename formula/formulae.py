# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

'''
Formula objects
===============

A formula is basically a sympy expression for the mean of something of
the form::

   mean = sum([Beta(e)*e for e in expr])

Or, a linear combination of sympy expressions, with each one multiplied
by its own "Beta". The elements of expr can be instances of Term (for a
linear regression formula, they would all be instances of Term). But, in
general, there might be some other parameters (i.e. sympy.Symbol
instances) that are not Terms.

The design matrix is made up of columns that are the derivatives of mean
with respect to everything that is not a Term, evaluted at a recarray
that has field names given by [str(t) for t in self.terms].

For those familiar with R's formula syntax, if we wanted a design matrix
like the following::

    > s.table = read.table("http://www-stat.stanford.edu/~jtaylo/courses/stats191/data/supervisor.table", header=T)
    > d = model.matrix(lm(Y ~ X1*X3, s.table)
    )
    > d
       (Intercept) X1 X3 X1:X3
    1            1 51 39  1989
    2            1 64 54  3456
    3            1 70 69  4830
    4            1 63 47  2961
    5            1 78 66  5148
    6            1 55 44  2420
    7            1 67 56  3752
    8            1 75 55  4125
    9            1 82 67  5494
    10           1 61 47  2867
    11           1 53 58  3074
    12           1 60 39  2340
    13           1 62 42  2604
    14           1 83 45  3735
    15           1 77 72  5544
    16           1 90 72  6480
    17           1 85 69  5865
    18           1 60 75  4500
    19           1 70 57  3990
    20           1 58 54  3132
    21           1 40 34  1360
    22           1 61 62  3782
    23           1 66 50  3300
    24           1 37 58  2146
    25           1 54 48  2592
    26           1 77 63  4851
    27           1 75 74  5550
    28           1 57 45  2565
    29           1 85 71  6035
    30           1 82 59  4838
    attr(,"assign")
    [1] 0 1 2 3
    >

With the Formula, it looks like this:

First read the same data as above:

>>> from os.path import dirname, join as pjoin
>>> import numpy as np
>>> import formula
>>> fname = pjoin(dirname(formula.__file__), 'data', 'supervisor.table')
>>> r = np.recfromtxt(fname, names=True)

Define the formula

>>> from formula import terms, Formula
>>> X1, X3 = terms(('X1', 'X3'))
>>> f = Formula([X1, X3, X1*X3, 1])
>>> f.mean
_b0*X1 + _b1*X3 + _b2*X1*X3 + _b3

The 1 is the "intercept" term, I have explicity not used R's default of adding
it to everything.

>>> f.design(r)
array([(51.0, 39.0, 1989.0, 1.0), (64.0, 54.0, 3456.0, 1.0),
       (70.0, 69.0, 4830.0, 1.0), (63.0, 47.0, 2961.0, 1.0),
       (78.0, 66.0, 5148.0, 1.0), (55.0, 44.0, 2420.0, 1.0),
       (67.0, 56.0, 3752.0, 1.0), (75.0, 55.0, 4125.0, 1.0),
       (82.0, 67.0, 5494.0, 1.0), (61.0, 47.0, 2867.0, 1.0),
       (53.0, 58.0, 3074.0, 1.0), (60.0, 39.0, 2340.0, 1.0),
       (62.0, 42.0, 2604.0, 1.0), (83.0, 45.0, 3735.0, 1.0),
       (77.0, 72.0, 5544.0, 1.0), (90.0, 72.0, 6480.0, 1.0),
       (85.0, 69.0, 5865.0, 1.0), (60.0, 75.0, 4500.0, 1.0),
       (70.0, 57.0, 3990.0, 1.0), (58.0, 54.0, 3132.0, 1.0),
       (40.0, 34.0, 1360.0, 1.0), (61.0, 62.0, 3782.0, 1.0),
       (66.0, 50.0, 3300.0, 1.0), (37.0, 58.0, 2146.0, 1.0),
       (54.0, 48.0, 2592.0, 1.0), (77.0, 63.0, 4851.0, 1.0),
       (75.0, 74.0, 5550.0, 1.0), (57.0, 45.0, 2565.0, 1.0),
       (85.0, 71.0, 6035.0, 1.0), (82.0, 59.0, 4838.0, 1.0)], 
      dtype=[('X1', '<f8'), ('X3', '<f8'), ('X1*X3', '<f8'), ('1', '<f8')])
'''

import numpy as np

import sympy

from .utils import contrast_from_cols_or_rows, make_dummy


class Beta(sympy.symbol.Dummy):
    ''' A dummy symbol tied to a Term `term` '''
    _beta_flag = True

    def __new__(cls, name, term):
        new = sympy.symbol.Dummy.__new__(cls, name)
        new._term = term
        return new


class Formula(object):
    """ A Formula is a model for a mean in a regression model.

    It is often given by a sequence of sympy expressions, with the mean
    model being the sum of each term multiplied by a linear regression
    coefficient.

    The expressions may depend on additional Symbol instances,
    giving a non-linear regression model.
    """
    # This flag is defined to avoid using isinstance
    _formula_flag = True

    def __init__(self, seq, char = 'b'):
        """
        Inputs:
        -------
        seq : [``sympy.Basic``]
        char : character for regression coefficient

        """
        self._terms = np.asarray(list(seq))
        if list(self._terms) == []:
            self._terms = np.array([sympy.Number(0)])
        self._counter = 0
        self.char = char

    # Properties
    def _getcoefs(self):
        if not hasattr(self, '_coefs'):
            self._coefs = []
            for term in self.terms:
                self._coefs.append(Beta("%s%d" % (self.char, self._counter), term))
                self._counter += 1
        return np.array(self._coefs)

    def _setcoefs(self, coefs):
        for i, term in enumerate(self.terms):
            if not is_beta(coefs[i]):
                raise ValueError('should be a linear regression coefficient')
            coefs[i].term = term
        self._coefs = np.array(coefs)
    coefs = property(_getcoefs, _setcoefs, doc='Coefficients in the linear regression formula.')

    def _getterms(self):
        return self._terms
    terms = property(_getterms, doc='Terms in the linear regression formula.')

    def _getmean(self):
        """ Expression for mean

        Expression for the mean, expressed as a linear combination of
        terms, each with dummy variables in front.
        """
        return np.sum(self.coefs*self.terms)

    mean = property(_getmean, doc="Expression for the mean, expressed "
                    "as a linear combination of terms, each with dummy "
                    "variables in front.")

    @property
    def unique(self):
        """
        Return a Formula(np.unique(self.terms))
        """
        return self.__class__(np.unique(self.terms))

    @property
    def params(self):
        """
        The parameters in the Formula.
        """
        from .parts import getparams
        return getparams(self.mean)

    def _getdiff(self):
        from .parts import getparams
        params = list(set(getparams(self.mean)))
        params.sort()
        return [sympy.diff(self.mean, p).doit() for p in params]
    design_expr = property(_getdiff)

    def _getdtype(self):
        vnames = [str(s) for s in self.design_expr]
        return np.dtype([(n, np.float) for n in vnames])
    dtype = property(_getdtype, doc='The dtype of the design matrix of the Formula.')

    def __getitem__(self, key):
        """ Return the term such that str(term) == key.

        Parameters
        ----------
        key : str
            name of term to retrieve

        Returns
        -------
        term : sympy.Expression
        """
        names = [str(t) for t in self.terms]
        try:
            idx = names.index(key)
        except ValueError:
            raise ValueError('term %s not found' % key)
        return self.terms[idx]

    def subs(self, old, new):
        """ Perform a sympy substitution on all terms in the Formula

        Returns a new instance of the same class

        Parameters
        ----------
        old : sympy.Basic
           The expression to be changed
        new : sympy.Basic
           The value to change it to.

        Returns
        -------
        newf : Formula

        Examples
        --------
        >>> from formula import terms
        >>> s, t = terms('s, t')
        >>> f, g = [sympy.Function(l) for l in 'fg']
        >>> form = Formula([f(t),g(s)])
        >>> newform = form.subs(g, sympy.Function('h'))
        >>> newform.terms
        array([f(t), h(s)], dtype=object)
        >>> form.terms
        array([f(t), g(s)], dtype=object)
        """
        return self.__class__([term.subs(old, new) for term in self.terms])

    def delete_terms(self, other):
        """
        """
        l1 = list(self.terms)
        l2 = list(other.terms)
        for t in l2:
            if t in l1:
                l1.remove(t)
        return self.__class__(l1)

    def __repr__(self):
        return """Formula(%s)""" % `list(self.terms)`

    def __add__(self, other):
        """
        Create a new Formula by combining terms
        of other with those of self.

        Examples
        --------
        >>> from formula import terms
        >>> x, y, z = terms('x, y, z')
        >>> f1 = Formula([x,y,z])
        >>> f2 = Formula([y])+I
        >>> f3=f1+f2
        >>> sorted(f1.terms)
        [x, y, z]
        >>> sorted(f2.terms)
        [1, y]
        >>> sorted(f3.terms)
        [1, x, y, z]
        >>>
        """
        if hasattr(other, 'formula'):
            other = other.formula
        if self.terms.shape == ():
            return other
        elif other.terms.shape == ():
            return self
        return self.__class__(np.unique(list(self.terms) + list(other.terms)))

    def __mul__(self, other):
        """
        Create a new Formula by combining terms
        of other with those of self.

        Examples
        --------
        >>> from formula import terms
        >>> x, y, z = terms('x, y, z')
        >>> f1 = Formula([x,y,z])
        >>> f2 = Formula([y])+I
        >>> f3=f1+f2
        >>> sorted(f1.terms)
        [x, y, z]
        >>> sorted(f2.terms)
        [1, y]
        >>> sorted(f3.terms)
        [1, x, y, z]
        """
        if hasattr(other, 'formula'):
            other = other.formula
        result = np.unique(np.outer(self.terms, other.terms).reshape(-1))
        return self.__class__([T for T in result if T != 0])

    def __array__(self):
        return self.terms

    def __eq__(self, other):
        if hasattr(other, 'formula'):
            other = other.formula
        return set(np.unique(self.terms)) == set(np.unique(other.terms))

    def _setup_design(self):
        """
        Create a callable object to evaluate the design matrix
        at a given set of parameter values to be specified by
        a recarray and observed Term values, also specified
        by a recarray.
        """
        d = self.design_expr

        # Before evaluating, we recreate the formula
        # with numbered terms, and numbered parameters.

        # This renaming has no impact on the
        # final design matrix as the
        # callable, self._f below, is a lambda
        # that does not care about the names of the terms.

        # First, find all terms in the mean expression,
        # and rename them in the form "__t%d__" with a
        # random offset.
        # This may cause a possible problem
        # when there are parameters named something like "__t%d__".
        # Using the random offset will minimize the possibility
        # of this happening.

        # This renaming is here principally because of the
        # intercept.

        random_offset = np.random.random_integers(low=0, high=2**30)

        from .parts import getterms, getparams, is_factor_term
        terms = getterms(self.mean)

        newterms = []
        for i, t in enumerate(terms):
            newt = sympy.Symbol("__t%d__" % (i + random_offset))
            for j, _ in enumerate(d):
                d[j] = d[j].subs(t, newt)
            newterms.append(newt)

        # Next, find all parameters that remain in the design expression.
        # In a standard regression model, there will be no parameters
        # because they will all be differentiated away in computing
        # self.design_expr. In nonlinear models, parameters will remain.

        params = getparams(self.design_expr)
        newparams = []
        for i, p in enumerate(params):
            newp = make_dummy("__p%d__" % (i + random_offset))
            for j, _ in enumerate(d):
                d[j] = d[j].subs(p, newp)
            newparams.append(newp)

        # These "implemented" functions are used for things like
        # the natural splines, etc. You can represent natural splines
        # with sympy but the expression is pretty awful.

        # Note that  ``d`` here is list giving the differentiation of the
        # expression for the mean.  self._f(...) therefore also returns a list
        self._f = sympy.lambdify(newparams + newterms, d, ("numpy"))

        # The input to self.design will be a recarray of that must
        # have field names that the Formula will expect to see.
        # However, if any of self.terms are FactorTerms, then the field
        # in the recarray will not actually be in the Term.
        #
        # For example, if there is a Factor 'f' with levels ['a','b'],
        # there will be terms 'f_a' and 'f_b', though the input to
        # design will have a field named 'f'. In this sense,
        # the recarray used in the call to self.design
        # is not really made up of terms, but "preterms".

        # In this case, the callable

        preterm = []
        for t in terms:
            if not is_factor_term(t):
                preterm.append(str(t))
            else:
                preterm.append(t.factor_name)
        preterm = list(set(preterm))

        # There is also an argument for parameters that are not
        # Terms.

        self._dtypes = {'param':np.dtype([(str(p), np.float) for p in params]),
                        'term':np.dtype([(str(t), np.float) for t in terms]),
                        'preterm':np.dtype([(n, np.float) for n in preterm])}

        self.__terms = terms

    def design(self,
               input,
               param=None,
               return_float=False,
               contrasts=None):
        """ Construct the design matrix, and optional contrast matrices.

        Parameters
        ----------
        input : np.recarray
           Recarray including fields needed to compute the Terms in
           getparams(self.design_expr).
        param : None or np.recarray
           Recarray including fields that are not Terms in
           getparams(self.design_expr)
        return_float : bool, optional
           If True, return a np.float array rather than a np.recarray
        contrasts : None or dict, optional
           Contrasts. The items in this dictionary should be (str,
           Formula) pairs where a contrast matrix is constructed for
           each Formula by evaluating its design at the same parameters
           as self.design. If not None, then the return_float is set to True.
        """
        from .parts import is_factor_term
        self._setup_design()

        preterm_recarray = input
        param_recarray = param

        # The input to design should have field names for all fields in self._dtypes['preterm']
        if not set(preterm_recarray.dtype.names).issuperset(self._dtypes['preterm'].names):
            raise ValueError("for term, expecting a recarray with "
                             "dtype having the following names: %s"
                             % `self._dtypes['preterm'].names`)
        # The parameters should have field names for all fields in self._dtypes['param']
        if param_recarray is not None:
            if not set(param_recarray.dtype.names).issuperset(self._dtypes['param'].names):
                raise ValueError("for param, expecting a recarray with "
                                 "dtype having the following names: %s"
                                 % `self._dtypes['param'].names`)
        # If the only term is an intercept,
        # the return value is a matrix of 1's.
        if list(self.terms) == [sympy.Number(1)]:
            a = np.ones(preterm_recarray.shape[0], np.float)
            if not return_float:
                a = a.view(np.dtype([('intercept', np.float)]))
            return a
        elif not self._dtypes['term']:
            raise ValueError("none of the expresssions are self.terms "
                             "are Term instances; shape of resulting "
                             "undefined")
        # The term_recarray is essentially the same as preterm_recarray,
        # except that all factors in self are expanded
        # into their respective binary columns.
        term_recarray = np.zeros(preterm_recarray.shape[0],
                                 dtype=self._dtypes['term'])
        for t in self.__terms:
            if not is_factor_term(t):
                term_recarray[t.name] = preterm_recarray[t.name]
            else:
                term_recarray['%s_%s' % (t.factor_name, t.level)] = \
                    np.array(map(lambda x: x == t.level, preterm_recarray[t.factor_name]))
        # The lambda created in self._setup_design needs to take a tuple of
        # columns as argument, not an ndarray, so each column
        # is extracted and put into float_tuple.
        float_array = term_recarray.view(np.float)
        float_array.shape = (term_recarray.shape[0], -1)
        float_array = float_array.T
        float_tuple = tuple(float_array)
        # If there are any parameters, they also must be extracted
        # and put into a tuple with the order specified
        # by self._dtypes['param']
        if param_recarray is not None:
            param = tuple(float(param_recarray[n]) for n in self._dtypes['param'].names)
        else:
            param = ()
        # Evaluate the design at the parameters and tuple of arrays
        D = self._f(*(param+float_tuple))
        # TODO: check if this next stepis necessary
        # I think it is because the lambda evaluates sympy.Number(1) to 1
        # and not an array.
        D_tuple = [np.asarray(w) for w in D]

        need_to_modify_shape = []
        OK_row_shapes = []
        for i, row in enumerate(D_tuple):
            if row.shape in [(),(1,)]:
                need_to_modify_shape.append(i)
            else:
                OK_row_shapes.append(row.shape[0])
        # Make sure that each array has the correct shape.
        # The columns in need_to_modify should just be
        # the intercept column, which evaluates to have shape == ().
        # This makes sure that it has the correct number of rows.
        for i in need_to_modify_shape:
            D_tuple[i].shape = ()
            D_tuple[i] = np.multiply.outer(D_tuple[i], np.ones(preterm_recarray.shape[0]))
        # At this point, all the columns have the correct shape and the
        # design matrix is almost ready to output.
        D = np.array(D_tuple).T
        # If we will return a float matrix or any contrasts,
        # we may have some reshaping to do.
        if contrasts is None:
            contrasts = {}
        if return_float or contrasts:
            # If the design matrix is just a column of 1s
            # return a 1-dimensional array.
            D = np.squeeze(D.astype(np.float))
            # If there are contrasts, the pseudo-inverse of D
            # must be computed.
            if contrasts:
                if D.ndim == 1:
                    _D = D.reshape((D.shape[0], 1))
                else:
                    _D = D
                pinvD = np.linalg.pinv(_D)
        else:
            # Correct the dtype.
            # XXX There seems to be a lot of messing around with the dtype.
            # This would be a convenient place to just add
            # labels like a DataArray.
            D = np.array([tuple(r) for r in D], self.dtype)
        # Compute the contrast matrices, if any.
        if contrasts:
            cmatrices = {}
            for key, cf in contrasts.items():
                if not is_formula(cf):
                    cf = self.__class__([cf])
                L = cf.design(input, param=param_recarray,
                              return_float=True)
                cmatrices[key] = contrast_from_cols_or_rows(L, _D, pseudo=pinvD)
            return D, cmatrices
        else:
            return D

# The intercept formula

I = Formula([sympy.Number(1)])

def is_beta(obj):
    """ Is obj a Beta?
    """
    return hasattr(obj, "_beta_flag")


def is_formula(obj):
    """ Is obj a Formula?
    """
    return hasattr(obj, "_formula_flag")




