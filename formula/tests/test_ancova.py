""" Testing ancova module
"""

from os.path import join as pjoin, dirname

import numpy as np
# recfromcsv available from numpy 1.3.0
from numpy import recfromcsv

from ..parts import Term, Factor, fromrec
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


def test_smoke():
    # smoke test, more or less
    salary_fname = pjoin(dirname(__file__), '..', 'data', 'salary.csv')
    d = recfromcsv(salary_fname)
    terms = fromrec(d)
    f = ANCOVA(1, terms['e'],terms['p'],(1,(terms['e'],terms['p'])))
    ANCOVA.add_intercept = False
    f2 = ANCOVA(terms['e'],(1,(terms['e'],terms['p'])))
    ANCOVA.add_intercept = True
    f3 = ANCOVA((1,(terms['e'],terms['p'])))
