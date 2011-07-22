""" Testing ancova module
"""

import numpy as np

from ..parts import Term, Factor
from ..ancova import ANCOVA

from numpy.testing import (assert_array_almost_equal,
                           assert_array_equal)

from nose.tools import assert_true, assert_equal, assert_raises

X = Term('X')
F = Factor('F', range(3))
H = Factor('H', range(2))

def test_init():
    a0 = ANCOVA((X,F), (X,(F,H)))
    assert_equal(a0.sequence(),
                 [(1,()), (X,(F,)), (X, (F,H))])


def test_delete_terms():
    a0 = ANCOVA((X,F), (X,(F,H)))
    a1 = a0.delete_terms((X,F))
    assert_equal(a1.sequence(),
                 [(1,()), (X, (F,H))])

