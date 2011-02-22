import numpy as np
from formula.terms import Factor
import sympy
from string import uppercase
from StringIO import StringIO
import matplotlib.mlab as ML
import nose.tools as nt

from os import remove
from formula.terms import Term, Factor
from formula.categorical import CategoricalFormula
import matplotlib.mlab as ML
import tempfile, rpy2.robjects

def random_letters(size, nlevels=6):
    return [uppercase[i] for i in np.random.random_integers(0,nlevels-1,size=size)]

def random_subset(subset, size):
    s = sorted(subset)
    return [s[i] for i in np.random.random_integers(0,len(s)-1,size=size)]

def random_recarray(size):
    initial = np.empty(size, dtype=np.dtype([('Y', np.float)]))
    initial['Y'] = np.random.standard_normal(size)
    numeric_vars = [np.random.standard_normal(size) for _ in range(10)]
    categorical_vars = [random_letters(size, l) for l in [3,4,7,6,4,5,8]]
    inter = ML.rec_append_fields(initial, ['n%s' % l for l in uppercase[:10]], numeric_vars)
    final = ML.rec_append_fields(inter, ['c%s' % l for l in uppercase[:len(categorical_vars)]], categorical_vars)
    return final, sympy.symbols(['n%s' % l for l in uppercase[:10]]), [Factor('c%s' % s, np.unique(l)) for s, l in zip(uppercase[:len(categorical_vars)],
                                                                                                                       categorical_vars)]

def random_categorical_formula(size=500):
    nterms = np.random.poisson(5)
    X, n, c = random_recarray(size)
    d = {}
    for _ in range(nterms):
        expr = random_subset(n, np.random.binomial(len(n), 0.5))
        f = []
        for _ in range(np.random.poisson(2)):
            factors = random_subset(c, np.random.poisson(1))
            f.append(np.unique(factors))
        d[np.product(np.unique(expr))] = f
    return CategoricalFormula(d)

def random_from_factor(factor, size):
    return random_subset(factor.levels, size)

def random_from_terms_factors(terms, factors, size):
    dtype = np.dtype([(str(t), np.float) for t in terms] + [(f.name,'S30') for f in factors])
    data = np.empty(size, 
                    np.dtype([(str(terms[0]), np.float)]))
    data[str(terms[0])] = np.random.standard_normal(size)
    for t in terms[1:]:
        data = ML.rec_append_fields(data, str(t), 
                                    np.random.standard_normal(size))
    for f in factors:
        data = ML.rec_append_fields(data, f.name, random_from_factor(f, size))
    return data

def random_from_categorical_formula(cf, size):
    exprs = []
    factors = []
    for key, value in cf.expr_factor_dict.items():
        if str(key) != '1':
            exprs += str(key).split("*")
        for fs in value:
            factors += list(fs)
    return random_from_terms_factors(list(set(exprs)), list(set(factors)), size)

def simple():

    x = Term('x'); y = Term('y') ; z = Term('z')
    f = Factor('f', ['a','b','c'])
    g = Factor('g', ['aa','bb','cc'])
    h = Factor('h', ['a','b','c','d','e','f','g','h','i','j'])
    d = CategoricalFormula({x*y:((g,f),),x:([f],[g,f]),
                            1:([g],),
                            z:([h],[h,g,f])})

    return d

def testR(d=simple(), size=500):

    X = random_from_categorical_formula(d, size)

    X = ML.rec_append_fields(X, 'response', np.random.standard_normal(size))
    fname = tempfile.mktemp()
    ML.rec2csv(X, fname)
    Rstr = '''
    data = read.table("%s", sep=',', header=T)
    cur.lm = lm(response ~ %s, data)
    COEF = coef(cur.lm)
    ''' % (fname, d.Rstr)
    rpy2.robjects.r(Rstr)
    remove(fname)
    nR = list(np.array(rpy2.robjects.r("names(COEF)")))

    nt.assert_true('(Intercept)' in nR)
    nR.remove("(Intercept)")
    nF = [str(t).replace("_","").replace("*",":") for t in d.formula.terms]
             
    nR = sorted([sorted(n.split(":")) for n in nR])

    nt.assert_true('1' in nF)
    nF.remove('1')

    nF = sorted([sorted(n.split(":")) for n in nF])
    nt.assert_equal(nR, nF)

    return d, X, nR, nF
