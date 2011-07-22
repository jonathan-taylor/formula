""" Testing random effects module
"""

import numpy as np

from ..parts import Factor, Term, FactorTerm
from ..convenience import make_recarray, terms
from ..random_effects import RandomEffects

from numpy.testing import (assert_array_almost_equal,
                           assert_array_equal)

from nose.tools import (assert_true, assert_equal, 
                        assert_raises)


def test_init():
    # Basic creation tests
    levels = [FactorTerm('s', L) for L in (2,3)]
    re = RandomEffects(levels)
    subj = make_recarray([2,2,2,3,3], 's')
    # The output is a bit difficult to test because of the sympy dummy variables
    c = re.cov(subj)
    # Specifying covariance directly
    re2 = RandomEffects(levels, sigma=[[4,1],[1,6]])
    c2 = re2.cov(subj)
    # Specifying covariance directly
    exp = np.ones((5,5))
    exp[:3,:3] = 4
    exp[3:,3:] = 6
    assert_array_equal(c2, exp)


