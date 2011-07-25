import numpy as np

import sympy
from sympy.utilities.lambdify import implemented_function, lambdify

from .parts import Factor, Term
from .formulae import Formula
from .utils import SYMPY_0p6


def make_recarray(rows, names, dtypes=None):
    """ Create recarray from `rows` with field names `names`

    Create a recarray with named columns from a list of rows and names
    for the columns. If dtype is None, the dtype is based on rows if it
    is an np.ndarray, else the data is cast as np.float. If dtypes are
    supplied, it uses the dtypes to create a np.dtype unless rows is an
    np.ndarray, in which case dtypes are ignored

    Parameters
    ----------
    rows: array-like
        Rows that will be turned into an recarray.
    names: sequence
        Sequence of strings - names for the columns.
    dtypes: [str or np.dtype]
        Used to create a np.dtype, can be np.dtypes or string.

    Returns
    -------
    v : np.ndarray

    Examples
    --------
    We've used the ``#doctest:+ELLIPSIS option here to allow us not to care
    either about the default integer size (i4 or i8), or the machine byte order
    (<,>), in the record array dtype.

    >>> arr = np.array([[3,4],[4,6],[6,8]])
    >>> make_recarray(arr, ['x','y']) # doctest:+ELLIPSIS
    array([[(3, 4)],
           [(4, 6)],
           [(6, 8)]], 
          dtype=[('x', ...), ('y', ...)])
    >>> r = make_recarray(arr, ['w', 'u'])
    >>> make_recarray(r, ['x','y']) # doctest:+ELLIPSIS
    array([[(3, 4)],
           [(4, 6)],
           [(6, 8)]], 
          dtype=[('x', ...), ('y', ...)])
    >>> make_recarray([[3,4],[4,6],[7,9]], 'wv', [np.float, np.int]) # doctest:+ELLIPSIS
    array([(3.0, 4), (4.0, 6), (7.0, 9)], 
          dtype=[('w', ...), ('v', ...)])
    """
    # XXX This function is sort of one of convenience
    # Would be nice to use DataArray or something like that
    # to add axis names.
    if isinstance(rows, np.ndarray):
        if rows.dtype.isbuiltin:
            dtype = np.dtype([(n, rows.dtype) for n in names])
        else:
            dtype = np.dtype([(n, d[1]) for n, d in zip(names, rows.dtype.descr)])
        if dtypes is not None:
            raise ValueError('dtypes not used if rows is an ndarray')
        return rows.view(dtype)
    if dtypes is None:
        dtype = np.dtype([(n, np.float) for n in names])
    else:
        dtype = np.dtype([(n, d) for n, d in zip(names, dtypes)])
    nrows = []
    vector = -1
    for r in rows:
        if vector < 0:
            a = np.array(r)
            if a.shape == ():
                vector = True
            else:
                vector = False
        if not vector:
            nrows.append(tuple(r))
        else:
            nrows.append(r)
    if vector:
        if len(names) != 1: # a 'row vector'
            nrows = tuple(nrows)
            return np.array(nrows, dtype)
        else:
            nrows = np.array([(r,) for r in nrows], dtype)
    return np.array(nrows, dtype)


def terms_from_rec(rec, keep=[], drop=[]):
    """ Construct terms from recarray

    For fields with a string-dtype, it is assumed that these are
    qualtiatitve regressors, i.e. Factors.

    Parameters
    ----------
    rec: recarray
        Recarray whose field names will be used to create a formula.
    keep: []
        Field names to explicitly keep, dropping all others.
    drop: []
        Field names to drop.
    """
    spec = {}
    for n in rec.dtype.names:
        if rec[n].dtype.kind == 'S':
            spec[n] = Factor.fromcol(rec[n], n)
        else:
            spec[n] = Term(n).formula
    for d in drop:
        del(spec[d])
    return spec


def natural_spline(t, knots=None, order=3, intercept=False):
    """ Return a Formula containing a natural spline

    Spline for a Term with specified `knots` and `order`.

    Parameters
    ----------
    t : ``Term``
    knots : None or sequence, optional
       Sequence of float.  Default None (same as empty list)
    order : int, optional
       Order of the spline. Defaults to a cubic (==3)
    intercept : bool, optional
       If True, include a constant function in the natural
       spline. Default is False

    Returns
    -------
    formula : Formula
         A Formula with (len(knots) + order) Terms
         (if intercept=False, otherwise includes one more Term), 
         made up of the natural spline functions.

    Examples
    --------
    The following results depend on machine byte order
       
    >>> x = Term('x')
    >>> n = natural_spline(x, knots=[1,3,4], order=3)
    >>> xval = np.array([3,5,7.]).view(np.dtype([('x', np.float)]))
    >>> n.design(xval, return_float=True)
    array([[   3.,    9.,   27.,    8.,    0.,   -0.],
           [   5.,   25.,  125.,   64.,    8.,    1.],
           [   7.,   49.,  343.,  216.,   64.,   27.]])
    >>> d = n.design(xval)
    >>> print d.dtype.descr
    [('ns_1(x)', '<f8'), ('ns_2(x)', '<f8'), ('ns_3(x)', '<f8'), ('ns_4(x)', '<f8'), ('ns_5(x)', '<f8'), ('ns_6(x)', '<f8')]
    >>>                    
                    
    """
    if knots is None:
        knots = {}
    fns = []
    for i in range(order+1):
        n = 'ns_%d' % i
        def f(x, i=i):
            return x**i
        s = implemented_function(n, f)
        fns.append(s(t))

    for j, k in enumerate(knots):
        n = 'ns_%d' % (j+i+1,)
        def f(x, k=k, order=order):
            return (x-k)**order * np.greater(x, k)
        s = implemented_function(n, f)
        fns.append(s(t))

    if not intercept:
        fns.pop(0)

    ff = Formula(fns)
    return ff


def define(name, expr):
    """ Create function of t expression from arbitrary expression `expr`

    Take an arbitrarily complicated expression `expr` of 't' and make it
    an expression that is a simple function of t, of form ``'%s(t)' %
    name`` such that when it evaluates (via ``lambdify``) it has the
    right values.

    Parameters
    ----------
    expr : sympy expression
       with only 't' as a Symbol
    name : str

    Returns
    -------
    nexpr: sympy expression

    Examples
    --------
    >>> from sympy import lambdify
    >>> t = Term('t')
    >>> expr = t**2 + 3*t
    >>> print expr
    t**2 + 3*t
    >>> newexpr = define('f', expr)
    >>> print newexpr
    f(t)
    >>> f = lambdify(t, newexpr)
    >>> f(4)
    28
    >>> 3*4+4**2
    28
    """
    # make numerical implementation of expression
    t = Term('t')
    v = lambdify(t, expr)
    # convert numerical implementation to sympy function
    f = implemented_function(name, v)
    # Return expression that is function of time
    return f(t)


def terms(names, **kwargs):
    ''' Return list of terms with names given by `names`

    This is just a convenience in defining a set of terms, and is the
    equivalent of ``sympy.symbols`` for defining symbols in sympy.

    We enforce the sympy 0.7.0 behavior of returning symbol "abc" from input
    "abc", rthan than 3 symbols "a", "b", "c".

    Parameters
    ----------
    names : str or sequence of str
       If a single str, can specify multiple ``Term``s with string
       containing space or ',' as separator. 
    \\**kwargs : keyword arguments
       keyword arguments as for ``sympy.symbols``

    Returns
    -------
    ts : ``Term`` or tuple
       ``Term`` instance or list of ``Term`` instance objects named from `names`

    Examples
    --------
    >>> terms(('a', 'b', 'c'))
    (a, b, c)
    >>> terms('a, b, c')
    (a, b, c)
    >>> terms('abc')
    abc
    '''
    if (SYMPY_0p6
        and isinstance(names, basestring)
        and not set(', ').intersection(names)):
        if not kwargs.get('each_char', False):
            # remove each_char (or no-op if absent)
            kwargs.pop('each_char', None)
            names = (names,)
    syms = sympy.symbols(names, **kwargs)
    try:
        len(syms)
    except TypeError:
        return Term(syms.name)
    return tuple(Term(s.name) for s in syms)

