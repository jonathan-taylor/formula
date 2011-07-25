from itertools import combinations

import numpy as np

import sympy

try:
    # matrix_rank in numpy >= 1.5.0
    from numpy.linalg import matrix_rank as rank
except ImportError:
    from numpy.lingalg import svd
    def rank(M, tol=None):
        """
        Return matrix rank of array using SVD method

        Rank of the array is the number of SVD singular values of the
        array that are greater than `tol`.

        Parameters
        ----------
        M : array_like
            array of <=2 dimensions
        tol : {None, float}
            threshold below which SVD values are considered zero. If `tol` is
            None, and ``S`` is an array with singular values for `M`, and
            ``eps`` is the epsilon value for datatype of ``S``, then `tol` is
            set to ``S.max() * eps``.
        """
        M = np.asarray(M)
        if M.ndim > 2:
            raise TypeError('array should have 2 or fewer dimensions')
        if M.ndim < 2:
            return int(not all(M==0))
        S = svd(M, compute_uv=False)
        if tol is None:
            tol = S.max() * np.finfo(S.dtype).eps
        return sum(S > tol)


def fullrank(X, r=None):
    """ Return a matrix whose column span is the same as X
    using an SVD decomposition.

    If the rank of X is known it can be specified by r-- no check is
    made to ensure that this really is the rank of X.
    """

    if r is None:
        r = rank(X)

    V, D, U = np.linalg.svd(X, full_matrices=0)
    order = np.argsort(D)
    order = order[::-1]
    value = []
    for i in range(r):
        value.append(V[:,order[i]])
    return np.asarray(np.transpose(value)).astype(np.float64)


def contrast_from_cols_or_rows(L, D, pseudo=None):
    """ Construct a contrast matrix from a design matrix D
    
    (possibly with its pseudo inverse already computed)
    and a matrix L that either specifies something in
    the column space of D or the row space of D.

    Parameters
    ----------
    L : ndarray
       Matrix used to try and construct a contrast.
    D : ndarray
       Design matrix used to create the contrast.
    pseudo : None or array-like, optional
       If not None, gives pseudo-inverse of `D`.  Allows you to pass
       this if it is already calculated. 
       
    Returns
    -------
    C : ndarray
       Matrix with C.shape[1] == D.shape[1] representing an estimable
       contrast.

    Notes
    -----
    From an n x p design matrix D and a matrix L, tries to determine a p
    x q contrast matrix C which determines a contrast of full rank,
    i.e. the n x q matrix

    dot(transpose(C), pinv(D))

    is full rank.

    L must satisfy either L.shape[0] == n or L.shape[1] == p.

    If L.shape[0] == n, then L is thought of as representing
    columns in the column space of D.

    If L.shape[1] == p, then L is thought of as what is known
    as a contrast matrix. In this case, this function returns an estimable
    contrast corresponding to the dot(D, L.T)

    This always produces a meaningful contrast, not always
    with the intended properties because q is always non-zero unless
    L is identically 0. That is, it produces a contrast that spans
    the column space of L (after projection onto the column space of D).
    """
    L = np.asarray(L)
    D = np.asarray(D)
    
    n, p = D.shape

    if L.shape[0] != n and L.shape[1] != p:
        raise ValueError, 'shape of L and D mismatched'

    if pseudo is None:
        pseudo = pinv(D)

    if L.shape[0] == n:
        C = np.dot(pseudo, L).T
    else:
        C = np.dot(pseudo, np.dot(D, L.T)).T
        
    Lp = np.dot(D, C.T)

    if len(Lp.shape) == 1:
        Lp.shape = (n, 1)
        
    if rank(Lp) != Lp.shape[1]:
        Lp = fullrank(Lp)
        C = np.dot(pseudo, Lp).T

    return np.squeeze(C)


def simplicial_complex(*simplices):
  
    """
    Take a list of simplices and compute
    its the simplicial complex generated
    by these simplices, returning
    the maximal simplices of this complex.
    
    >>> faces, maximal, all = simplicial_complex([('a','b','c'), ('c','d'), ('e','f'), ('e',)])
    >>> faces
    {1: set([('a',), ('c',), ('b',), ('e',), ('d',), ('f',)]), 2: set([('b', 'c'), ('a', 'b'), ('e', 'f'), ('c', 'd'), ('a', 'c')]), 3: set([('a', 'b', 'c')])}
    >>> maximal
    [('a', 'b', 'c'), ('e', 'f'), ('c', 'd')]
    >>> all
    set([('b', 'c'), ('a',), ('c',), ('c', 'd'), ('e', 'f'), ('a', 'c'), ('d',), ('b',), ('f',), ('a', 'b'), ('e',), ('a', 'b', 'c')])
    >>> 

    """

    faces = {0:set([])}

    if not simplices:
        return faces, (), set([set([])])
    else:
        l = [len(list(x)) for x in simplices]
        lmax = max(l)
        for i in range(lmax):
            faces[i+1] = set([])

        for simplex, lsimplex in zip(simplices, l):
            simplex = sorted(simplex)
            for k in range(lsimplex+1):
                for v in combinations(simplex, k):
                    faces[k].add(v)
        # now determine the maximal faces

        maximal = list(faces[lmax])
        allf = list(faces[lmax])
        for i in sorted(faces.keys())[-2::-1]:
            allf += list(faces[i])
            for simplex in faces[i]:
                keep = True
                for msimplex in maximal:
                    if set(simplex).issubset(msimplex):
                        keep = False
                if keep:
                    maximal.append(simplex)
        return faces, maximal[::-1], allf[::-1]

def factor_codings(*factor_monomials):
    """ Find which factors to code with indicator or contrast variables

    Determine which factors to code with indicator variables (using
    len(factor.levels) columns of 0s and 1s) or contrast coding (using
    len(factor.levels)-1).  The elements of the sequence should be tuples of
    strings.  Further, the factors are assumed to be in *graded* order, that is
    [len(f) for f in factor_monomials] is assumed non-decreasing.

    Examples
    --------
    >>> factor_codings(('b',), ('a',), ('b', 'c'), ('a','b','c'))
    {('b', 'c'): [('b', 'indicator'), ('c', 'contrast')], ('a',): [('a', 'contrast')], ('b',): [('b', 'indicator')], ('a', 'b', 'c'): [('a', 'contrast'), ('b', 'indicator'), ('c', 'indicator')]}
    >>> factor_codings(('a',), ('b',), ('b', 'c'), ('a','b','c'))
    {('b', 'c'): [('b', 'indicator'), ('c', 'contrast')], ('a',): [('a', 'indicator')], ('b',): [('b', 'contrast')], ('a', 'b', 'c'): [('a', 'contrast'), ('b', 'indicator'), ('c', 'indicator')]}

    Here is a version with debug strings to see what is happening:

    >>> factor_codings(('a',), ('b', 'c'), ('a','b','c')) #doctest: +SKIP
    Adding a from ('a',) as indicator because we have not seen any factors yet.
    Adding b from ('b', 'c') as indicator because set([('c',), ()]) is not a subset of set([(), ('a',)])
    Adding c from ('b', 'c') as indicator because set([(), ('b',)]) is not a subset of set([(), ('a',)])
    Adding a from ('a', 'b', 'c') as contrast because set([('c',), ('b', 'c'), (), ('b',)]) is a subset of set([('b', 'c'), (), ('c',), ('b',), ('a',)])
    Adding b from ('a', 'b', 'c') as indicator because set([('c',), (), ('a', 'c'), ('a',)]) is not a subset of set([('b', 'c'), (), ('c',), ('b',), ('a',)])
    Adding c from ('a', 'b', 'c') as indicator because set([('a', 'b'), (), ('b',), ('a',)]) is not a subset of set([('b', 'c'), (), ('c',), ('b',), ('a',)])
    {('b', 'c'): [('b', 'indicator'), ('c', 'indicator')], ('a',): [('a', 'indicator')], ('a', 'b', 'c'): [('a', 'contrast'), ('b', 'indicator'), ('c', 'indicator')]}

    Notes
    -----
    Even though the elements of factor_monomials are assumed to be in graded
    order, the final result depends on the ordering of the strings of the
    factors within each of the tuples.
    """
    lmax = 0
    from copy import copy
    already_seen = set([])
    final_result = []
    for factor_monomial in factor_monomials:
        result = []
        factor_monomial = list(factor_monomial)
        if len(factor_monomial) < lmax:
            raise ValueError('factors are assumed to be in graded order')
        lmax = len(factor_monomial)

        for j in range(len(factor_monomial)):
            cur = copy(list(factor_monomial))
            cur.pop(j)
            terms = simplicial_complex(cur)[2]
            if already_seen and set(terms).issubset(already_seen):
                result.append((factor_monomial[j], 'contrast'))
            else:
                result.append((factor_monomial[j], 'indicator'))
        already_seen = already_seen.union(simplicial_complex(factor_monomial)[2])
        final_result.append((tuple(factor_monomial), result))
    return dict(final_result)


def make_dummy(name):
    """ Make dummy variable of given name

    Parameters
    ----------
    name : str
        name of dummy variable

    Returns
    -------
    dum : `Dummy` instance

    Notes
    -----
    The interface to Dummy changed between 0.6.7 and 0.7.0
    """
    from distutils.version import LooseVersion
    if LooseVersion(sympy.__version__) >= LooseVersion('0.7.0'):
        return sympy.Dummy(name)
    return sympy.Symbol(name, dummy=True)


class BombError(Exception):
    pass


class Bomber(object):
    """ Raise exception for any attribute access """
    def __init__(self, msg):
        self.msg = msg
    def __getattribute__(self, name):
        raise BombError(self.msg)
