""" Testing formula utils
"""

import numpy as np

from .. utils import fullrank

from numpy.testing import (assert_array_almost_equal,
                           assert_almost_equal,
                           assert_array_equal)

from nose.tools import assert_true, assert_equal, assert_raises


def test_fullrank():
    X = np.random.standard_normal((40,5))
    X[:,0] = X[:,1] + X[:,2]
    Y1 = fullrank(X)
    assert_equal(Y1.shape, (40,4))
    Y2 = fullrank(X, r=3)
    assert_equal(Y2.shape, (40,3))
    Y3 = fullrank(X, r=4)
    assert_equal(Y3.shape, (40,4))
    assert_almost_equal(Y1, Y3)

