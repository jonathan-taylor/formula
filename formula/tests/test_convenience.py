""" Testing
"""

import numpy as np

from ..convenience import terms
from ..parts import Term

from numpy.testing import (assert_array_almost_equal,
                           assert_array_equal)

from nose.tools import assert_true, assert_equal, assert_raises


def test_terms():
    t = terms('a')
    assert_true(isinstance(t, Term))
    a, b, c = Term('a'), Term('b'), Term('c')
    assert_equal(t, a)
    ts = terms(('a', 'b', 'c'))
    assert_equal(ts, (a, b, c))
    # a string without separator chars returns one symbol.  This is the
    # future sympy default.
    assert_equal(terms('abc'), Term('abc'))
    # separators return multiple symbols
    assert_equal(terms('a b c'), (a, b, c))
    assert_equal(terms('a, b, c'), (a, b, c))
    # no arg is an error
    assert_raises(TypeError, terms)
    # but empty arg returns empty tuple
    assert_equal(terms(()), ())
    # Test behavior of deprecated each_char kwarg
    try:
        res = terms('abc', each_char=False)
    except TypeError:
        return
    assert_equal(res, Term('abc'))
    assert_equal(terms('abc', each_char=True), (a, b, c))
