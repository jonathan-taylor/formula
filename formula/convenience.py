import numpy as np
from itertools import combinations

from .parts import Factor, Term
from .formulae import Formula


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
    The following tests depend on machine byte order to pass
    
    >>> arr = np.array([[3,4],[4,6],[6,8]])
    >>> make_recarray(arr, ['x','y'])
    array([[(3, 4)],
           [(4, 6)],
           [(6, 8)]], 
          dtype=[('x', '<i8'), ('y', '<i8')])
    >>> r = make_recarray(arr, ['w', 'u'])
    >>> make_recarray(r, ['x','y'])
    array([[(3, 4)],
           [(4, 6)],
           [(6, 8)]], 
          dtype=[('x', '<i8'), ('y', '<i8')])
    >>> make_recarray([[3,4],[4,6],[7,9]], 'wv', [np.float, np.int])
    array([(3.0, 4), (4.0, 6), (7.0, 9)], 
          dtype=[('w', '<f8'), ('v', '<i8')])
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
        s = aliased_function(n, f)
        fns.append(s(t))

    for j, k in enumerate(knots):
        n = 'ns_%d' % (j+i+1,)
        def f(x, k=k, order=order):
            return (x-k)**order * np.greater(x, k)
        s = aliased_function(n, f)
        fns.append(s(t))

    if not intercept:
        fns.pop(0)

    ff = Formula(fns)
    return ff


def define(name, expr):
    """
    Take an expression of 't' (possibly complicated)
    and make it a '%s(t)' % name, such that
    when it evaluates it has the right values.

    Parameters
    ----------
    expr : sympy expression, with only 't' as a Symbol
    name : str

    Returns
    -------
    nexpr: sympy expression

    Examples
    --------
    >>> t = Term('t')
    >>> expr = t**2 + 3*t
    >>> print expr
    3*t + t**2
    >>> newexpr = define('f', expr)
    >>> print newexpr
    f(t)
    >>> import aliased
    >>> f = aliased.lambdify(t, newexpr)
    >>> f(4)
    28
    >>> 3*4+4**2
    28
    """
    v = vectorize(expr)
    return aliased_function(name, v)(Term('t'))


def terms(*names):
    ''' Return list of terms with names given by `names`

    This is just a convenience in defining a set of terms, and is the
    equivalent of ``sympy.symbols`` for defining symbols in sympy. 

    Parameters
    ----------
    *names : str or sequence of str
       If a single str, can specify multiple ``Term``s with string
       containing space or ',' as separator. 

    Returns
    -------
    ts : list
       list of Term instance objects named from `names`

    Examples
    --------
    >>> terms('a', 'b', 'c')
    (a, b, c)
    >>> terms('a, b, c')
    (a, b, c)
    '''
    # parse separated single string
    if len(names) == 1:
        name = names[0]
        if isinstance(name, basestring):
            for sep in ', ':
                if sep in name:
                    names = (n.strip() for n in name.split(sep))
                    break
    return tuple(Term(n) for n in names)

