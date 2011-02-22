Categorical Formulae
====================

Example
~~~~~~~

This is an example with one numeric variable
and 2 categorical variables. It can be found 
`here <http://stats191.stanford.edu/interactions.html>`_
though this example has modified the data to be characters
instead of numeric. This just makes *R* automatically
treat some things as factors.


.. rcode::

   url = 'http://stats191.stanford.edu/data/salary.csv'
   data = read.table(url, header=T, sep=',')

   salary.lm = lm(S ~ X:E + X:P, data=data)
   print(coef(salary.lm))

The order of the two terms above matters

.. rcode::

   salary2.lm = lm(S ~ X:P + X:E, data=data)
   print(coef(salary2.lm))

.. testsetup::

   from formula.CategoricalFormula import CategoricalFormula
   from formula.terms import Term, Factor

We can recreate this example using the *CategoricalFormula* object

.. testcode::

   e = Factor('E', ['B', 'M', 'P'])
   p = Factor('P', ['M', 'L'])
   x = Term("X")

   f = CategoricalFormula({x:[(e,),(p,)]})
   print sorted(f.formula.terms)

This would output:

.. testoutput::

   [1, E_B*X, E_M*X, E_P*X, P_M*X]

Interchanging the order above

.. testcode::

   f2 = CategoricalFormula({x:[(p,),(e,)]})
   print sorted(f2.formula.terms)

.. testoutput::

   [1, E_M*X, E_P*X, P_L*X, P_M*X]

The difference in the order has to do with how *R* constructs
a design matrix from a set of factors (and the corresponding
numeric variable).

Codings
~~~~~~~

The difference manifests itself in how each factor
multiplying *X* above is *coded* in each
product. A factor appearing in a product of factors
can be coded either as in *indicator* which
includes all relevant dummy variables for the
factor or as a *contrast* which includes one less
column. The default behavior of *R* is to drop a term,
which *CategoricalFormula* also follows. *R* drops
the first term in a sorted list of the levels of the factor,
*CategoricalFormula* follows the same behavior by default.

For instance:

.. testcode::

   print 'f', f.codings
   print 'f2', f2.codings

would yield

.. testoutput::

   {1: {}, X: {('E',): [('E', 'indicator')], ('P',): [('P', 'contrast')]}}
   {1: {}, X: {('E',): [('E', 'contrast')], ('P',): [('P', 'indicator')]}}

In the first formula, *M* appears as a contrast while *E* appears
as an indicator. In the second formula, *E* appears as a contrast
and *M* as an indicator.

Two-way interactions
~~~~~~~~~~~~~~~~~~~~

We can also specify two-way interactions in *R*

.. rcode::

   print(coef(lm(S ~ X:E:P, data=data)))

As well as in the *CategoricalFormula*

.. testcode::

   f3 = CategoricalFormula({x:[(p,e)]})
   print f3.codings
   print sorted(f3.formula.terms)

which yields

.. testoutput::

   {1: {}, X: {('E', 'P'): [('E', 'indicator'), ('P', 'indicator')]}}
   [1, E_B*P_L*X, E_B*P_M*X, E_M*P_L*X, E_M*P_M*X, E_P*P_L*X, E_P*P_M*X]

If we add in some parts of the formula, it becomes a little harder to
predict:

.. rcode::

   print(coef(lm(S ~ X:E:P + X:P + X:E, data=data)))

.. testcode::

   f4 = CategoricalFormula({x:[(p,e),(p,),(e,)]})
   print f4.codings
   print sorted(f4.formula.terms)

which yields

.. testoutput::

   {1: {}, X: {('E', 'P'): [('E', 'contrast'), ('P', 'contrast')], ('E',): [('E', 'contrast')], ('P',): [('P', 'indicator')]}}
   [1, E_M*X, E_P*X, P_L*X, P_M*X, E_M*P_M*X, E_P*P_M*X]

which agrees with *R*.

Changing the order above again changes the terms in the
formula

.. testcode::

   f5 = CategoricalFormula({x:[(p,e),(e,),(p,)]})
   print f5.codings
   print sorted(f5.formula.terms)

as it does in *R*:

.. rcode::

   print(coef(lm(S ~ X:E:P + X:E + X:P, data=data)))

What is important is the *graded* order. That is, for the numeric
variable *X*, the first order factors are ordered in *f4* as 
*[set([P]),set([E])]* and its
second order factors are *[set([P,E])]* while it has no zeroth order
factors. The only difference between *f4* and *f5* is the order
of its first order factors.

Adding *X* to the *R* formula adds a zeroth order factor.

.. rcode::

   print(coef(lm(S ~ X + X:E:P + X:E + X:P, data=data)))

With the categorical formula, this can be achieved
by

.. testcode::

   f6 = CategoricalFormula({x:[(), (p,e),(e,),(p,)]})
   print f6.codings
   print sorted(f6.formula.terms)

which yields

.. testoutput::

   {1: {}, X: {(): [], ('E', 'P'): [('E', 'contrast'), ('P', 'contrast')], ('P',): [('P', 'contrast')], ('E',): [('E', 'contrast')]}}
   [1, E_M*X, E_P*X, P_M*X, E_M*P_M*X, E_P*P_M*X, X]

Intercept
~~~~~~~~~

The ubiquitous intercept can be suppressed using
the keyword argument "add_intercept" to the
constructor of CategoricalFormula

.. testcode::

   f7 = CategoricalFormula({x:[(), (p,e),(e,),(p,)]}, add_intercept=False)
   print f7.codings
   print sorted(f7.formula.terms)

which yields

.. testoutput::

   {X: {(): [], ('E', 'P'): [('E', 'contrast'), ('P', 'contrast')], ('P',): [('P', 'contrast')], ('E',): [('E', 'contrast')]}}
   [E_M*X, E_P*X, P_M*X, E_M*P_M*X, E_P*P_M*X, X]

In *R* the intercept can be removed (most of the time) by appending *-1*
to the string specifying the formula:

.. rcode::

   print(coef(lm(S ~ X + X:E:P + X:E + X:P - 1, data=data)))

This is not quite the same as obtained by *CategoricalFormula* and 
this is related to *R*'s treatment of factors and numeric variables
as equal. The *CategoricalFormula* makes a distinction
between these two.

Contrasts
~~~~~~~~~

Each *(expr, factor)* pair in the *CategoricalFormula* specification
maps to a specific contrast.

.. testcode::

   f7 = CategoricalFormula({x:[(), (p,e),(e,),(p,)]}, add_intercept=False)
   print f7.contrasts

.. testoutput::

   {'I(X)': Formula([X]),
    'I(X):E': Formula([E_P*X, E_M*X]),
    'I(X):E:P': Formula([E_M*P_M*X, E_P*P_M*X]),
    'I(X):P': Formula([P_M*X])}

As opposed to

.. testcode::

   f3 = CategoricalFormula({x:[(p,e)]})
   print f3.contrasts

which yields

.. testoutput::

   {'I(1):1': Formula([1]),
    'I(X):E:P': Formula([E_B*P_L*X, E_P*P_M*X, E_M*P_L*X, E_P*P_L*X, E_B*P_M*X, E_M*P_M*X])}

Slices
~~~~~~

Each contrast can be associated with some columns of the
final design matrix. These are also elements
of the *formula* attribute

.. testcode::

   f3 = CategoricalFormula({x:[(p,e)]})
   print f3.slices

which yields

.. testoutput::

   {'I(1):1': slice(0, 1, None), 'I(X):E:P': slice(1, 7, None)}

Contrast Matrices
~~~~~~~~~~~~~~~~~

The slices can be interpreted as contrast matrices

.. testcode::

   f3 = CategoricalFormula({x:[(p,e)]})
   print f3.slices

yielding

.. testoutput::

   {'I(1):1': array([[ 1.,  0.,  0.,  0.,  0.,  0.,  0.]]),
    'I(X):E:P': array([[ 0.,  1.,  0.,  0.,  0.,  0.,  0.],
	  [ 0.,  0.,  1.,  0.,  0.,  0.,  0.],
	  [ 0.,  0.,  0.,  1.,  0.,  0.,  0.],
	  [ 0.,  0.,  0.,  0.,  1.,  0.,  0.],
	  [ 0.,  0.,  0.,  0.,  0.,  1.,  0.],
	  [ 0.,  0.,  0.,  0.,  0.,  0.,  1.]])}

Note, however, that not all these contrasts are estimable. Whether
they are estimable or not depends on the actual
design matrix used to fit an OLS model. Users should keep this in mind.
In this example, the contrast *I(X):E:P* would not be estimable
if we never observed a laborer with a PhD, for example.
