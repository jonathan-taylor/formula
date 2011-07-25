""" Testing parts of formulas
"""

import numpy as np

from ..parts import fromrec, is_term, is_factor

from numpy.testing import (assert_array_almost_equal,
                           assert_array_equal)

from nose.tools import assert_true, assert_equal, assert_raises


def test_fromrec():
    # Test sttructured arrays to terms and factors
    D = np.rec.array([
            (43, 51, 30, 39, 1, 92, 'blue'),
            (63, 64, 51, 54, 3, 73, 'blue'),
            (71, 70, 68, 69, 6, 86, 'red'),
            (61, 63, 45, 47, 4, 84, 'red'),
            (81, 78, 56, 66, 1, 83, 'blue'),
            (43, 55, 49, 44, 4, 49, 'blue'),
            (58, 67, 42, 56, 6, 68, 'green'),
            (71, 75, 50, 55, 0, 66, 'green'),
            (72, 82, 72, 67, 1, 83, 'blue'),
            (67, 61, 45, 47, 2, 80, 'red'),
            (64, 53, 53, 58, 8, 67, 'blue'),
            (67, 60, 47, 39, 9, 74, 'green'),
            (69, 62, 57, 42, 5, 63, 'blue'),
            (68, 83, 83, 45, 9, 77, 'red'),
            (77, 77, 54, 72, 9, 77, 'red'),
            (81, 90, 50, 72, 0, 54, 'blue'),
            (74, 85, 64, 69, 9, 79, 'green'),
            (65, 60, 65, 75, 5, 80, 'green'),
            (65, 70, 46, 57, 5, 85, 'red'),
            (50, 58, 68, 54, 4, 78, 'red'),
            (50, 40, 33, 34, 3, 64, 'blue'),
            (64, 61, 52, 62, 6, 80, 'blue'),
            (53, 66, 52, 50, 3, 80, 'red'),
            (40, 37, 42, 58, 0, 57, 'red'),
            (63, 54, 42, 48, 6, 75, 'blue'),
            (66, 77, 66, 63, 8, 76, 'blue'),
            (78, 75, 58, 74, 0, 78, 'red'),
            (48, 57, 44, 45, 1, 83, 'blue'),
            (85, 85, 71, 71, 7, 74, 'red'),
            (82, 82, 39, 59, 4, 78, 'blue')], 
                     dtype=[('y', 'i8'),
                            ('x1', 'i8'),
                            ('x2', 'i8'),
                            ('x3', 'i8'),
                            ('x4', 'O'),
                            ('x5', 'i8'),
                            ('x6', '|S5')])
    fts = fromrec(D)
    for t in ['y', 'x1', 'x2', 'x3', 'x5']:
        assert_true(is_term(fts[t]))
    for t in ['x4', 'x6']:
        assert_true(is_factor(fts[t]))
    assert_equal(fts['x4'].levels, range(10))
    assert_equal(fts['x6'].levels, ['blue', 'green', 'red'])
