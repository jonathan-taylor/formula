""" Testing ancova module
"""

from os.path import join as pjoin, dirname

import numpy as np
# recfromcsv available from numpy 1.3.0
from numpy import recfromcsv

from ..parts import Term, Factor, fromrec
from ..ancova import ANCOVA, typeI, typeII, typeIII

from numpy.testing import (assert_array_almost_equal,
                           assert_array_equal)

from nose.tools import assert_true, assert_equal, assert_raises

X = Term('X')
F = Factor('F', range(3))
H = Factor('H', range(2))
SALARY = recfromcsv(pjoin(dirname(__file__), '..', 'data', 'salary.csv'))

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
    terms = fromrec(SALARY)
    f = ANCOVA(1, terms['e'],terms['p'],(1,(terms['e'],terms['p'])))
    ANCOVA.add_intercept = False
    f2 = ANCOVA(terms['e'],(1,(terms['e'],terms['p'])))
    ANCOVA.add_intercept = True
    f3 = ANCOVA((1,(terms['e'],terms['p'])))

def test_multiply_by_factor():
    terms = fromrec(SALARY)
    f = ANCOVA(1, terms['e'])
    f2 = f.multiply_by_factor(terms['p'])
    assert_equal(ANCOVA(1, terms['p'], (1,(terms['e'], terms['p']))), f2)
    


def test_types123():
    terms = fromrec(SALARY)
    x = terms['x']; e = terms['e']; p = terms['p']
    ancova = ANCOVA((x,e),(x,p),(x,(p,e)))
    res1 = typeI('s', ancova, SALARY)
    res2 = typeII('s', ancova, SALARY)
    res3 = typeIII('s', ancova, SALARY)
    # Reversing the order changes the ANOVA tables, in particular
    # the degrees of freedom associated to each contrast. This is
    # because the codings change when the order of the factors change.
    ancova2 = ANCOVA((x,p),(x,e), (x,(p,e)))
    rres2 = typeII('s', ancova2, SALARY)
    rres3 = typeIII('s', ancova2, SALARY)
