""" Testing random effects module
"""

import numpy as np
from sympy import Symbol

from ..parts import Factor, Term, FactorTerm
from ..convenience import make_recarray, terms
from ..random_effects import RandomEffects

from numpy.testing import (assert_array_almost_equal,
                           assert_almost_equal,
                           assert_array_equal)

from nose.tools import (assert_true, assert_equal,
                        assert_raises)


def test_random_effects():
    # Basic creation tests
    levels = [FactorTerm('s', L) for L in (2,3)]
    re = RandomEffects(levels)
    subj = make_recarray([2,2,2,3,3], 's')
    # The output is a bit difficult to test because of the sympy dummy variables
    c = re.cov(subj)
    # Getting levels via a factor
    subj_factor = Factor('s', [2,3])
    # Specifying covariance directly
    re = RandomEffects(subj_factor.terms, sigma=np.array([[4,1],[1,6]]))
    C = re.cov(subj)
    assert_almost_equal(C, [[4,4,4,1,1],
                            [4,4,4,1,1],
                            [4,4,4,1,1],
                            [1,1,1,6,6],
                            [1,1,1,6,6]])
    a = Symbol('a')
    b = Symbol('b')
    re = RandomEffects(subj_factor.terms, sigma=np.array([[a,0],[0,b]]))
    C = re.cov(subj)
    # we return the symbols * 1.0.  The 1.0 disappears in sympy 0.6.x, but does
    # not in sympy 0.7.0, so we have to put it in for the comparision
    expected = np.array([[a,a,a,0,0],
                        [a,a,a,0,0],
                        [a,a,a,0,0],
                        [0,0,0,b,b],
                        [0,0,0,b,b]]) * 1.0
    assert_array_equal(C, expected)

