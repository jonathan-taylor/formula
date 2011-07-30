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

Let's try this with the ANCOVA model.  First we do some set up:

.. ipython:: python

   import numpy as np
   # Set precision printing to match R output
   np.set_printoptions(precision=4)
   import matplotlib.mlab as ML
   from pkg_resources import resource_stream
   salary = ML.csv2rec(resource_stream("formula", "data/salary.csv"))
   from formula.parts import fromrec
   from formula.ancova import *

Next we set up the model:

.. ipython:: python

   terms = fromrec(salary)
   x = terms['x']; e = terms['e']; p = terms['p']
   ancova = ANCOVA((x,e),(x,p))
   formula = ancova.formula
   print formula

To fit the model, we use *scikits.statsmodels*

.. ipython:: python

   from scikits.statsmodels.api import OLS
   model = OLS(salary['s'], formula.design(salary, return_float=True))
   results = model.fit()
   print results.params

Interchanging the order above

.. ipython:: python

   ancova2 = ANCOVA((x,p),(x,e))
   ancova2.formula.terms
   f2 = ancova2.formula
   model = OLS(salary['s'], f2.design(salary, return_float=True))
   results = model.fit()
   results.params

The difference in the order has to do with how *R* constructs
a design matrix from a set of factors (and the corresponding
numeric variable).

Factor Codings
~~~~~~~~~~~~~~

The difference manifests itself in how each factor
multiplying *X* above is *coded* in each
product. A factor appearing in a product of factors
can be coded either as in *indicator* which
includes all relevant dummy variables for the
factor or as a *contrast* which includes one less
column. The default behavior of *R* is to drop a term,
which *ANCOVA* also follows. *R* drops
the first term in a sorted list of the levels of the factor,
*ANCOVA* follows the same behavior by default.

For instance:

.. ipython:: python

   ancova.codings
   ancova2.codings

In the first formula, *P* appears as a contrast while *E* appears
as an indicator. In the second formula, *E* appears as a contrast
and *P* as an indicator.

Two-way interactions
~~~~~~~~~~~~~~~~~~~~

We can also specify two-way interactions in *R*

.. rcode::

   print(coef(lm(S ~ X:P:E, data=data)))

As well as in the *ANCOVA*

.. ipython:: python

   ancova3 = ANCOVA((x,(p,e)))
   ancova3.codings
   f3 = ancova3.formula
   model = OLS(salary['s'], f3.design(salary, return_float=True))
   results = model.fit()
   results.params


If we add in some parts of the formula, it becomes a little harder to
predict:

.. rcode::

   print(coef(lm(S ~ X:E:P + X:P + X:E, data=data)))

.. ipython:: python

   ancova4 = ANCOVA((x,(p,e)),(x,p),(x,e))
   ancova4.codings
   ancova4.formula.terms
   model = OLS(salary['s'], ancova4.formula.design(salary, return_float=True))
   results = model.fit()
   results.params


Changing the order above again changes the terms in the
formula

.. ipython:: python

   ancova5 = ANCOVA((x,(p,e)),(x,e),(x,p))
   ancova5.codings
   ancova.formula.terms
   model = OLS(salary['s'], ancova5.formula.design(salary, return_float=True))
   results = model.fit()
   results.params

as it does in *R*:

.. rcode::

   print(coef(lm(S ~ X:E:P + X:E + X:P, data=data)))

What is important is the *graded* order. That is, for the numeric
variable *X*, the first order factors are ordered in *f4* as
*[set([P]),set([E])]* and its
second order factors are *[set([P,E])]* while it has no zeroth order
factors. The only difference between *ancova4* and *ancova5* is the order
of its first order factors.

Adding *X* to the *R* formula adds a zeroth order factor.

.. rcode::

   print(coef(lm(S ~ X + X:E:P + X:E + X:P, data=data)))

With the categorical formula, this can be achieved
by

.. ipython:: python

   ancova6 = ANCOVA(x,(x,e),(x,p),(x,(p,e)))
   ancova6.codings
   ancova6.formula.terms
   model = OLS(salary['s'], ancova6.formula.design(salary, return_float=True))
   results = model.fit()
   results.params

One more example

.. rcode::

   print(coef(lm(S ~ X:E:P + X:E, data=data)))

.. ipython:: python

   ancova6a = ANCOVA((x,(e,p)),(x,e))
   ancova6a.codings
   ancova6a.formula.terms
   model = OLS(salary['s'], ancova6a.formula.design(salary, return_float=True))
   results = model.fit()
   results.params


Intercept
~~~~~~~~~

The ubiquitous intercept can be suppressed using
the keyword argument "add_intercept" to the
constructor of ANCOVA

.. ipython:: python

   ancova7 = ANCOVA(x,(x,(p,e)),(x,e),(x,p), add_intercept=False)
   ancova7.formula.terms
   model = OLS(salary['s'], ancova7.formula.design(salary, return_float=True))
   results = model.fit()
   results.params

In *R* the intercept can be removed (most of the time) by appending *-1*
to the string specifying the formula:

.. rcode::

   print(coef(lm(S ~ X + X:P:E + X:E + X:P - 1, data=data)))

This design matrix is not the same as obtained by *ANCOVA*,
hence, the coefficients are also different.
This is related to *R*'s treatment of factors and numeric variables
as equal. The *ANCOVA* module makes a distinction
between these two. The reason *R* has a missing value in the coefficients
is that its rules for generating design matrices told it that *E* should
be coded with indicators in the term *X:E* which leads
to a linear dependence with *X* already in the model.
The *ANCOVA* implementation treats *X* as *(X,1)* and hence when *(X,E)*
is to be added it sees that there will be a linear dependence if
*E* is added with indicator functions. Effectively, all columns with
*X* in them are the product of the columns of a
purely categorical formula. In this case, the columns
are the same as

.. ipython:: python

   ancova7a = ANCOVA((1,(p,e)), (1,e), (1,p))
   ancova7a.formula.terms
   ancova7a.formula.terms * x
   ancova7.formula.terms

This is how the *ANCOVA* is constructed. For each numeric term,
there is a corresponding pure categorical formula. For example

.. ipython:: python

   z = Term('z')
   ancova7b = ANCOVA((1,e), (z,e), (z,(e,p)), (x*z,e), (x,e), (x,p), x*z)
   ancova7b.sequence(x)
   ancova7b.sequence(z*x)
   ancova7b.sequence(1)
   ancova7b.sequence(z)

Any of those sequences above can be used to create new ANCOVA instances
whose formulae is that numeric expression multiplied by the corresponding
purely categorical formula.

.. ipython:: python

   ANCOVA(*ancova7b.sequence(z)).formula.terms
   purely_categorical = ANCOVA(*[(1, factors) for _, factors in ancova7b.sequence(z)])
   purely_categorical.formula.terms
   purely_categorical.formula.terms * z




Contrasts
~~~~~~~~~

Each *(expr, factor)* pair in the *ANCOVA* specification
maps to a specific contrast.

.. ipython:: python

   ancova7.contrasts

As opposed to

.. ipython:: python

   ancova3.contrasts

These contrasts are the default contrasts that
drop the first level of the factor. This can be changed
with the *default_contrast* keyword argument

.. ipython:: python

   ancova8 = ANCOVA(x,(x,(p,e)),(x,e),(x,p), default_contrast='main_effect')
   ancova8.contrasts

Contrast Matrices & Slices
~~~~~~~~~~~~~~~~~~~~~~~~~~

Each contrast can be associated with some columns of the
final design matrix. These are also elements
of the *formula* attribute

.. ipython:: python

   ancova3.slices


The slices can be interpreted as contrast matrices

.. ipython:: python

   ancova3.contrast_matrices

Note, however, that these contrast matrices depend on the *default_coding*
argument. Generally speaking, they are appropriate for use when
the *default_coding* is "main_effect" rather than "drop_reference".
*TODO: construct these properly for different default coding*

Further, not all these contrasts are estimable.
. Whether
they are estimable or not depends on the actual
design matrix used to fit an OLS model. Users should keep this in mind.
In this example, the contrast *I(X):E:P* would not be estimable
if we never observed a laborer with a PhD, for example.

Sums of squares
~~~~~~~~~~~~~~~

.. ipython:: python

   ancova = ANCOVA((x,e),(x,p),(x,(p,e)))
   print ML.rec2txt(typeI('s', ancova, salary))

Compare this to the R output

.. rcode::

   anova(lm(S ~ X:E + X:P + X:P:E, data=data))


For type II:


.. ipython:: python

   print ML.rec2txt(typeII('s', ancova, salary))


.. rcode::

   library(car)
   Anova(lm(S ~ X:E + X:P + X:P:E, data=data), type='II')

And type III:


.. ipython:: python

   print ML.rec2txt(typeIII('s', ancova, salary))


.. rcode::

   library(car)
   Anova(lm(S ~ X:E + X:P + X:P:E, data=data), type='III')

Reversing the order changes the ANOVA tables, in particular
the degrees of freedom associated to each contrast. This is
because the codings change when the order of the factors change.

.. ipython:: python

   ancova2 = ANCOVA((x,p),(x,e), (x,(p,e)))
   print ML.rec2txt(typeII('s', ancova2, salary))

.. rcode::

   library(car)
   Anova(lm(S ~ X:P + X:E + X:P:E, data=data), type='II')

.. ipython:: python

   print ML.rec2txt(typeIII('s', ancova2, salary))

.. rcode::

   library(car)
   Anova(lm(S ~ X:P + X:E + X:P:E, data=data), type='III')
