import sympy, numpy as np
from formula import Formula

class Term(sympy.Symbol):
    """A sympy.Symbol type to represent a term an a regression model
    
    Terms can be added to other sympy expressions with the single
    convention that a term plus itself returns itself.

    It is meant to emulate something on the right hand side of a formula
    in R. In particular, its name can be the name of a field in a
    recarray used to create a design matrix.

    >>> t = Term('x')
    >>> xval = np.array([(3,),(4,),(5,)], np.dtype([('x', np.float)]))
    >>> f = t.formula
    >>> d = f.design(xval)
    >>> print d.dtype.descr
    [('x', '<f8')]
    >>> f.design(xval, return_float=True)
    array([ 3.,  4.,  5.])
    """
    # This flag is defined to avoid using isinstance in getterms
    # and getparams.
    _term_flag = True

    @property
    def formula(self):
        """
        Return a Formula with only terms=[self].
        """
        return Formula([self])

class FactorTerm(Term):
    """ Boolean Term derived from a Factor.

    Its properties are the same as a Term except that its product with
    itself is itself.
    """
    # This flag is defined to avoid using isinstance in getterms
    # and _poly_from_formula
    _factor_term_flag = True

    def __new__(cls, name, level):
        new = Term.__new__(cls, "%s_%s" % (name, level))
        new.level = level
        new.factor_name = name
        return new

    def __mul__(self, other):

        if is_factor_term(other) and self.factor_name == other.factor_name:
            if self.level == other.level:
                return self
            else:
                return 0
        else:
            return sympy.Symbol.__mul__(self, other)

class Beta(sympy.symbol.Dummy):
    ''' A symbol tied to a Term `term` '''
    _beta_flag = True

    def __new__(cls, name, term):
        new = sympy.symbol.Dummy.__new__(cls, name)
        new._term = term
        return new

class Factor(object):
    """ A qualitative variable in a regression model
    
    A Factor is similar to R's factor. The levels of the Factor can be
    either strings or ints.
    """

    # This flag is defined to avoid using isinstance in getterms
    # and getparams.
    _factor_flag = True

    def __init__(self, name, levels, char='b'):
        """
        Parameters
        ----------
        name : str
        levels : [str or int]
            A sequence of strings or ints.
        char : str

        Returns
        -------
        """
        # Check whether they can all be cast to strings or ints without
        # loss.
        levelsarr = np.asarray(levels)
        if levelsarr.ndim == 0 and levelsarr.dtype.kind == 'S':
            levelsarr = np.asarray(list(levels))
        
        if levelsarr.dtype.kind != 'S': # the levels are not strings
            if not np.alltrue(np.equal(levelsarr, np.round(levelsarr))):
                raise ValueError('levels must be strings or ints')
            levelsarr = levelsarr.astype(np.int)
            
        self.levels = list(levelsarr)
        self.name = name
        self._char = char

    def __getitem__(self, level):
        """
        self.get_term(level)
        """
        return self.get_term(level)

    def __repr__(self):
        return "Factor('%s', %s)" % (self.name, `self.levels`)

    def get_term(self, level):
        """
        Retrieve a term of the Factor...
        """
        if level not in self.levels:
            raise ValueError('level not found')
        return self.formula["%s_%s" % (self.name, str(level))]


    # TODO: allow different specifications of the contrasts
    # here.... main_effect is like R's contr.sum

    def _getmaineffect(self, ref=-1):
        if not hasattr(self, "_main_effect"):
            terms = [FactorTerm(self.name, l) for l in 
                     self.levels]
            ref_term = terms[ref]
            terms.pop(ref)
            self._main_effect = Formula([term - ref_term for term in terms])
        return self._main_effect
    main_effect = property(_getmaineffect)

    @property
    def drop_first(self):
        ref = 0
        if not hasattr(self, "_drop1"):
            terms = [FactorTerm(self.name, l) for l in 
                     self.levels]
            ref_term = terms[ref]
            terms.pop(ref)
            self._drop1 = Formula(terms)
        return self._drop1


    def stratify(self, variable):
        """ Create a new variable, stratified by the levels of a Factor.

        Parameters
        ----------
        variable : str or a simple sympy expression whose string representation
            are all lower or upper case letters, i.e. it can be interpreted
            as a name

        Returns
        -------
        formula : Formula
            Formula whose mean has one parameter named variable%d, for each
            level in self.levels
        
        Examples
        --------
        >>> f = Factor('a', ['x','y'])
        >>> sf = f.stratify('theta')
        >>> sf.mean
        _theta0*a_x + _theta1*a_y
        """
        if not set(str(variable)).issubset(lowercase +
                                           uppercase + '0123456789'):
            raise ValueError('variable should be interpretable as a '
                             'name and not have anything but digits '
                             'and numbers')
        variable = sympy.sympify(variable)
        f = Formula(self.formula.terms, char=variable)
        f.name = self.name
        return f

    @property
    def formula(self):
        """
        Create a formula.
        """
        if not hasattr(self, "_formula"):
            self._formula = Formula([FactorTerm(self.name, l) for l in 
                                     self.levels], char=self._char)
        return self._formula

    @staticmethod
    def fromcol(col, name):
        """ Create a Factor from a column array.

        Parameters
        ----------
        col : ndarray
            an array with ndim==1

        name : str
            name of the Factor

        Returns
        -------

        factor : Factor

        Examples
        --------
        >>> data = np.array([(3,'a'),(4,'a'),(5,'b'),(3,'b')], np.dtype([('x', np.float), ('y', 'S1')]))
        >>> f1 = Factor.fromcol(data['y'], 'y')
        >>> f2 = Factor.fromcol(data['x'], 'x')
        >>> d = f1.design(data)
        >>> print d.dtype.descr
        [('y_a', '<f8'), ('y_b', '<f8')]
        >>> d = f2.design(data)
        >>> print d.dtype.descr
        [('x_3', '<f8'), ('x_4', '<f8'), ('x_5', '<f8')]
        """
        col = np.asarray(col)
        if col.ndim != 1 or (col.dtype.names and len(col.dtype.names) > 1):
            raise ValueError('expecting an array that can be thought '
                             'of as a column or field of a recarray')
        levels = np.unique(col)
        if not col.dtype.names and not name:
            name = 'factor'
        elif col.dtype.names:
            name = col.dtype.names[0]
        return Factor(name, levels)


def is_term(obj):
    """ Is obj a Term?
    """
    return hasattr(obj, "_term_flag")


def is_factor_term(obj):
    """ Is obj a FactorTerm?
    """
    return hasattr(obj, "_factor_term_flag")

def is_factor(obj):
    """ Is obj a Formula?
    """
    return hasattr(obj, "_factor_flag")
