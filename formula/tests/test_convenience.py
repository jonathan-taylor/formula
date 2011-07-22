""" Testing
"""

import numpy as np

from ..convenience import terms
from ..parts import Term

from numpy.testing import (assert_array_almost_equal,
                           assert_array_equal)

from nose.tools import assert_true, assert_equal, assert_raises


def test_terms():
    t, = terms('a')
    a, b, c = Term('a'), Term('b'), Term('c')
    assert_equal(t, a)
    ts = terms('a', 'b', 'c')
    assert_equal(ts, (a, b, c))
    # a string without separator chars returns one symbol.  This is the
    # future sympy default. 
    assert_equal(terms('abc'), (Term('abc'),))
    # separators return multiple symbols
    assert_equal(terms('a b c'), (a, b, c))
    assert_equal(terms('a, b, c'), (a, b, c))
    # nothing returns empty
    assert_equal(terms(), ())
