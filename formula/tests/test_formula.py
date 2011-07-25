# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Test functions in formula.formulae
"""

import numpy as np

import sympy
from sympy.utilities.lambdify import implemented_function

from .. import formulae as F
from ..parts import Term, Factor, stratify, fromrec
from ..convenience import make_recarray

from nose.tools import (assert_true, assert_equal, assert_false,
                        assert_raises)
from numpy.testing import assert_almost_equal 


def test_getparams_terms():
    t = Term('t')
    x, y, z = [sympy.Symbol(l) for l in 'xyz']
    assert_equal(set(F.getparams(x*y*t)), set([x,y]))
    assert_equal(set(F.getterms(x*y*t)), set([t]))
    matrix_expr = np.array([[x,y*t],[y,z]])
    assert_equal(set(F.getparams(matrix_expr)), set([x,y,z]))
    assert_equal(set(F.getterms(matrix_expr)), set([t]))


def test_formula_params():
    t = Term('t')
    x, y = sympy.symbols('x, y')
    f = F.Formula([t*x,y])
    assert_equal(set(f.params), set([x,y] + list(f.coefs)))


def test_contrast1():
    x = Term('x')
    assert_equal(x, x+x)
    y = Term('y')
    z = Term('z')
    f = F.Formula([x,y])
    arr = make_recarray([[3,5,4],[8,21,-1],[4,6,-2]], 'xyz')
    D, C = f.design(arr, contrasts={'x':x.formula,
                                    'diff':F.Formula([x-y]),
                                    'sum':F.Formula([x+y]),
                                    'both':F.Formula([x-y,x+y])})
    assert_almost_equal(C['x'], np.array([1,0]))
    assert_almost_equal(C['diff'], np.array([1,-1]))
    assert_almost_equal(C['sum'], np.array([1,1]))
    assert_almost_equal(C['both'], np.array([[1,-1],[1,1]]))
    f = F.Formula([x,y,z])
    arr = make_recarray([[3,5,4],[8,21,-1],[4,6,-2]], 'xyz')
    D, C = f.design(arr, contrasts={'x':x.formula,
                                    'diff':F.Formula([x-y]),
                                    'sum':F.Formula([x+y]),
                                    'both':F.Formula([x-y,x+y])})
    assert_almost_equal(C['x'], np.array([1,0,0]))
    assert_almost_equal(C['diff'], np.array([1,-1,0]))
    assert_almost_equal(C['sum'], np.array([1,1,0]))
    assert_almost_equal(C['both'], np.array([[1,-1,0],[1,1,0]]))


def test_design_expression():
    t1 = Term("x")
    t2 = Term('y')
    f = t1.formula + t2.formula
    assert_true(str(f.design_expr) in ['[x, y]', '[y, x]'])


def test_formula_property():
    # Check that you can create a Formula with one term
    t1 = Term("x")
    f = t1.formula
    assert_equal(f.design_expr, [t1])


def test_mul():
    f = Factor('t', [2,3])
    f2 = Factor('t', [2,3,4])
    t2, t3 = f.terms
    x = Term('x')
    assert_equal(t2, t2*t2)
    assert_equal(f, f*f)
    assert_false(f == f2)
    assert_equal(set((t2*x).atoms()), set([t2,x]))


def test_str_formula():
    t1 = Term('x')
    t2 = Term('y')
    f = F.Formula([t1, t2])
    assert_equal(str(f), "Formula([x, y])")


def test_design():
    # Check that you get the design matrix we expect
    t1 = Term("x")
    t2 = Term('y')

    n = make_recarray([2,4,5], 'x')
    assert_almost_equal(t1.formula.design(n)['x'], n['x'])

    f = t1.formula + t2.formula
    n = make_recarray([(2,3),(4,5),(5,6)], 'xy')

    assert_almost_equal(f.design(n)['x'], n['x'])
    assert_almost_equal(f.design(n)['y'], n['y'])

    f = t1.formula + t2.formula + F.I + t1.formula * t2.formula
    assert_almost_equal(f.design(n)['x'], n['x'])
    assert_almost_equal(f.design(n)['y'], n['y'])
    assert_almost_equal(f.design(n)['1'], 1)
    assert_almost_equal(f.design(n)['x*y'], n['x']*n['y'])
    # drop x field, check that design raises error
    ny = np.recarray(n.shape, dtype=[('x', n.dtype['x'])])
    ny['x'] = n['x']
    assert_raises(ValueError, f.design, ny)
    n = np.array([(2,3,'a'),(4,5,'b'),(5,6,'a')], np.dtype([('x', np.float),
                                                            ('y', np.float),
                                                            ('f', 'S1')]))
    f = Factor('f', ['a','b'])
    ff = t1.formula * f + F.I
    assert_almost_equal(ff.design(n)['f_a*x'], n['x']*[1,0,1])
    assert_almost_equal(ff.design(n)['f_b*x'], n['x']*[0,1,0])
    assert_almost_equal(ff.design(n)['1'], 1)


def test_implemented():
    x = Term('x')
    f = implemented_function('f', lambda x: 2*x)
    g = implemented_function('g', lambda x: np.sqrt(x))
    ff = F.Formula([f(x), g(x)**2])
    n = make_recarray([2,4,5], 'x')
    assert_almost_equal(ff.design(n)['f(x)'], n['x']*2)
    assert_almost_equal(ff.design(n)['g(x)**2'], n['x'])


def test_factor_getterm():
    fac = Factor('f', 'ab')
    assert_equal(fac['f_a'], fac.get_term('a'))
    fac = Factor('f', [1,2])
    assert_equal(fac['f_1'], fac.get_term(1))
    fac = Factor('f', [1,2])
    assert_raises(ValueError, fac.get_term, '1')
    m = fac.main_effect
    assert_equal(set(m.terms), set([fac['f_1']-fac['f_2']]))


def test_stratify():
    fac = Factor('x', [2,3])
    y = sympy.Symbol('y')
    f = sympy.Function('f')
    assert_raises(ValueError, stratify, fac, f(y))


def test_nonlin1():
    # Fit an exponential curve, with the exponent stratified by a factor
    # with a common intercept and multiplicative factor in front of the
    # exponential
    x = Term('x')
    fac = Factor('f', 'ab')
    f = F.Formula([sympy.exp(stratify(fac, x).mean)]) + F.I
    params = F.getparams(f.mean)
    assert_equal(set([str(p) for p in params]), set(['_x0', '_x1', '_b0', '_b1']))
    test1 = set(['1',
                 'exp(_x0*f_a + _x1*f_b)',
                 '_b0*f_a*exp(_x0*f_a + _x1*f_b)',
                 '_b0*f_b*exp(_x0*f_a + _x1*f_b)'])
    test2 = set(['1',
                 'exp(_x0*f_a + _x1*f_b)',
                 '_b1*f_a*exp(_x0*f_a + _x1*f_b)',
                 '_b1*f_b*exp(_x0*f_a + _x1*f_b)'])
    assert_true(test1 or test2)
    n = make_recarray([(2,3,'a'),(4,5,'b'),(5,6,'a')], 'xyf', ['d','d','S1'])
    p = make_recarray([1,2,3,4], ['_x0', '_x1', '_b0', '_b1'])
    A = f.design(n, p)
    print A, A.dtype


def test_intercept():
    dz = make_recarray([2,3,4],'z')
    v = F.I.design(dz, return_float=False)
    assert_equal(v.dtype.names, ('intercept',))


def test_nonlin2():
    dz = make_recarray([2,3,4],'z')
    z = Term('z')
    t = sympy.Symbol('th')
    p = make_recarray([3], ['tt'])
    f = F.Formula([sympy.exp(t*z)])
    assert_raises(ValueError, f.design, dz, p)


def test_Rintercept():
    x = Term('x')
    y = Term('x')
    xf = x.formula
    yf = y.formula
    newf = (xf+F.I)*(yf+F.I)
    assert_equal(set(newf.terms), set([x,y,x*y,sympy.Number(1)]))


def test_return_float():
    x = Term('x')
    f = F.Formula([x,x**2])
    xx= make_recarray(np.linspace(0,10,11), 'x')
    dtype = f.design(xx).dtype
    assert_equal(set(dtype.names), set(['x', 'x**2']))
    dtype = f.design(xx, return_float=True).dtype
    assert_equal(dtype, np.float)


def test_subtract():
    x, y, z = [Term(l) for l in 'xyz']
    f1 = F.Formula([x,y])
    f2 = F.Formula([x,y,z])
    f3 = f2 - f1
    assert_equal(set(f3.terms), set([z]))
    f4 = F.Formula([y,z])
    f5 = f1 - f4
    assert_equal(set(f5.terms), set([x]))


def test_subs():
    t1 = Term("x")
    t2 = Term('y')
    z = Term('z')
    f = F.Formula([t1, t2])
    g = f.subs(t1, z)
    assert_equal(list(g.terms), [z, t2])


