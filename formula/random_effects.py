import numpy as np

from .formulae import Formula
from .utils import make_dummy


class RandomEffects(Formula):
    """ Covariance matrices for common random effects analyses.

    Examples
    --------
    >>> from formula import Factor, make_recarray

    Make a subject factor with two levels

    >>> subj_factor = Factor('s', [2,3])

    Make corresponding data vector

    >>> subj = make_recarray([2,2,2,3,3], 's')

    Make random effects object with default covariance

    >>> c = RandomEffects(subj_factor.terms)
    >>> c.cov(subj)
    array([[1.0*_s2_0, 1.0*_s2_0, 1.0*_s2_0, 0, 0],
           [1.0*_s2_0, 1.0*_s2_0, 1.0*_s2_0, 0, 0],
           [1.0*_s2_0, 1.0*_s2_0, 1.0*_s2_0, 0, 0],
           [0, 0, 0, 1.0*_s2_1, 1.0*_s2_1],
           [0, 0, 0, 1.0*_s2_1, 1.0*_s2_1]], dtype=object)

    Specify the covariance

    >>> c = RandomEffects(subj_factor.terms,
    ...                   sigma=np.array([[4,1],[1,6]]))
    >>> c.cov(subj)
    array([[ 4.,  4.,  4.,  1.,  1.],
           [ 4.,  4.,  4.,  1.,  1.],
           [ 4.,  4.,  4.,  1.,  1.],
           [ 1.,  1.,  1.,  6.,  6.],
           [ 1.,  1.,  1.,  6.,  6.]])
    """
    def __init__(self, seq, sigma=None, char = 'e'):
        """
        Parameters
        ----------
        seq : [``sympy.Basic``]
        sigma : ndarray
             Covariance of the random effects. Defaults
             to a diagonal with entries for each random
             effect.
        char : character for regression coefficient
        """
        # XXX - char not used anywhere
        self._terms = np.asarray(seq)
        q = self._terms.shape[0]
        self._counter = 0
        if sigma is None:
            self.sigma = np.diag([make_dummy('s2_%d' % i) for i in range(q)])
        else:
            self.sigma = np.asarray(sigma)
        if self.sigma.shape != (q,q):
            raise ValueError('incorrect shape for covariance '
                             'of random effects, '
                             'should have shape %s' % repr(q,q))
        self.char = char

    def cov(self, term, param=None):
        """
        Compute the covariance matrix for
        some given data.

        Parameters:
        -----------
        term : np.recarray
             Recarray including fields corresponding to the Terms in 
             getparams(self.design_expr).
        param : np.recarray
             Recarray including fields that are not Terms in 
             getparams(self.design_expr)

        Outputs:
        --------
        C : ndarray
             Covariance matrix implied by design and self.sigma.
        """
        D = self.design(term, param=param, return_float=True)
        return np.dot(D, np.dot(self.sigma, D.T))

