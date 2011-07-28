# playing with designs - not a test
from ..convenience import terms_from_rec
from ..ancova import ANCOVA
from .test_design import random_recarray

X, n, c =random_recarray(400)
s=terms_from_rec(X)
f=ANCOVA((s['nG']['nG'], (s['cA'], s['cB'])), (s['nG']['nG'],(s['cD'],)))
# This used to exist, but does no more
# f.Rformula()
