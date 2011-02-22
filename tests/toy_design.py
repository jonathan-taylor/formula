
X=random_recarray(400)
_ip.magic("run convenience.py")
s=terms_from_rec(X)
_ip.magic("run categorical.py")
f=CategoricalFormula((s['nG']['nG'], (s['cA'], s['cB'])), (s['nG']['nG'],(s['cD'],)))
f.Rformula()
f=CategoricalFormula((s['nG']['nG'], (s['cA'], s['cB'])), (s['nG']['nG'],(s['cD'],)))
