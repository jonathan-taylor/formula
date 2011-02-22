import sympy
from numpy.random import random_integers
from random import sample
from string import lowercase
from itertools import combinations

class Monoid(object):

    flags = set(['monoid'])
    MI = 1

    @classmethod
    def product(cls, a, b):
        """
        The product of two elements of the monoid.
        """
        raise NotImplementedError

    @staticmethod
    def validate(a):
        """
        Is a an element of the ring?
        Raise an exception if not,
        else returns a canonical representative for a.
        """
        raise NotImplementedError

    @classmethod
    def generic(cls):
        """
        Return a generic element of the monoid.
        """
        raise NotImplentedError

class Ring(Monoid):

    flags = set(['monoid', 'ring'])

    # the additive identity of the ring
    AI = None

    # the multiplicative identity of the ring
    MI = None

    @classmethod
    def sum(cls, a, b):
        """
        The sum of two elements of the ring.
        """
        raise NotImplementedError

    @classmethod
    def difference(cls, a, b):
        """
        The difference of two elements of the ring.
        """
        return cls.sum(a, cls.inverse(b))

    @classmethod
    def product(cls, a, b):
        """
        The product of two elements of the ring.
        """
        raise NotImplementedError

    @staticmethod
    def validate(a):
        """
        Is a an element of the ring?
        Raise an exception if not,
        else returns a canonical representative for a.
        """
        raise NotImplementedError

    @classmethod
    def inverse(cls, a):
        """
        The additive inverse of a.
        """
        raise NotImplentedError
    
    @classmethod
    def generic(cls):
        """
        Return a generic element of the ring.
        """
        raise NotImplentedError

class CommutativeRing(Ring):

    flags = set(['monoid', 'ring', 'commutative_monoid'])

class Z(CommutativeRing):

    AI = 0
    MI = 1

    @classmethod
    def sum(cls, a, b):
        a = cls.validate(a)
        b = cls.validate(b)
        return a + b

    @classmethod
    def product(cls, a, b):
        a = cls.validate(a)
        b = cls.validate(b)
        return a * b

    @staticmethod
    def validate(a):
        return int(a)

    @classmethod
    def inverse(cls, a):
        a = cls.validate(a)
        return - a

    high = 1e7
    low = -1e7
    @classmethod
    def generic(cls):
        return random_integers(cls.low, cls.high)

class Z2(Z):

    @classmethod
    def sum(cls, a, b):
        a = cls.validate(a)
        b = cls.validate(b)
        return a + b

    @classmethod
    def product(cls, a, b):
        a = cls.validate(a)
        b = cls.validate(b)
        return a * b

    @staticmethod
    def validate(a):
        return int(a) % 2

    @classmethod
    def inverse(cls, a):
        a = cls.validate(a)
        return - a

    @classmethod
    def generic(cls):
        return random_integers(0, 1)


class SympyRing(CommutativeRing):
    """
    The sympy ring.
    """

    AI = sympy.Number(0)
    MI = sympy.Number(1)

    @classmethod
    def sum(cls, a, b):
        a = cls.validate(a)
        b = cls.validate(b)
        return sympy.Add(a, b)

    @classmethod
    def product(cls, a, b):
        a = cls.validate(a)
        b = cls.validate(b)
        return sympy.Mul(a, b)

    @staticmethod
    def validate(a):
        return sympy.sympify(a)

    @classmethod
    def inverse(cls, a):
        a = cls.validate(a)
        return sympy.Mul(-1,a)

    @classmethod
    def generic(cls, monomial=False):
        """
        Return a generic polynomial in sympy variables.
        If cls.monomial_generic is True, then return
        a monomial.
        """

        def _monomial():
            s = sympy.Number(1)
            for v in sample(lowercase, random_integers(0,3)):
                s *= sympy.Symbol(v)**random_integers(1,3)
            return s
        if monomial:
            return _monomial()
        else:
            s = 1
            for _ in range(random_integers(0,4)):
                s += _monomial()
            return s
        
class SubsetMonoid(Monoid):

    """
    Sets are represented as sorted tuples
    because they must be hashable to collect terms in
    the monoid ring.
    """

    MI = ()

    @classmethod
    def product(cls, a, b):
        a = cls.validate(a)
        b = cls.validate(b)
        return tuple(sorted((set(a).union(b))))

    @staticmethod
    def validate(a):
        return tuple(sorted(set(a)))

    @classmethod
    def generic(cls):
        def _random_string():
            return ''.join([lowercase[random_integers(0,25)] for _ in range(4)])
        return tuple([_random_string() for _ in range(random_integers(0,3))])

"""
A monoid ring has a ring as well as a monoid.
This is how we define addition, subtraction
and multiplication.
"""

def sum(monoid_ring, a, b):
    """
    Return the sum of two elements in the monoid ring.
    """
    return simplify(monoid_ring.ring, *collect_terms(monoid_ring.ring, *(list(a)+list(b))))


def difference(monoid_ring, a, b):
    """
    The difference of a and b in the monoid ring. This
    is just sum(a, monoid_ring.scalar_multiply(ring.inverse(ring.multiplicative_identity), b)
    """
    return simplify(monoid_ring.ring,
                    *collect_terms(monoid_ring.ring,
                                   *(list(a)+list([(monoid_ring.ring.inverse(r), m)
                                                   for r, m in b]))))

def product(monoid_ring, a, b):
    """
    Return the product of two elements in the monoid ring.
    """
    result = []
    for ring_a, monoid_a in a:
        for ring_b, monoid_b in b:
            result.append([monoid_ring.ring.product(ring_a, ring_b),
                           monoid_ring.monoid.product(monoid_a, monoid_b)])
    return simplify(monoid_ring.ring, *collect_terms(monoid_ring.ring, *result))

def scalar_multiply(monoid_ring, ring_element, *terms):
    """
    Multiply each term's ring coefficient by f on the left.
    """
    return tuple([(monoid_ring.ring.product(ring_element, r), m) for r, m in terms])

def monoid_multiply(monoid_ring, monoid_element, *terms):
    """
    Multiply each term's monoid coefficient by f on the left.
    """
    return tuple([(r, monoid_ring.monoid.product(monoid_element, m))
                  for r, m in terms])

def generic(monoid_ring, monomial=False):
    """

    """
    if not monomial:
        return tuple([(monoid_ring.ring.generic(),
                       monoid_ring.monoid.generic()) for
                      _ in range(random_integers(0,3))])
    else:
        return (monoid_ring.ring.generic(),
                monoid_ring.monoid.generic())

class MonoidRing(Ring):

    """
    Elements of the monoid ring are of the form

    [(r1,m1), (r2,m2), (r3,m3)]

    where each r* is in the ring of the MonoidRing
    and each m* is in the monoid of the MonoidRing.

    """
    def __new__(cls, monoid, ring):

        assert 'ring' in ring.flags
        assert 'monoid' in monoid.flags

        return type("MonoidRing", (Ring,),
                    {'monoid':monoid,
                     'ring':ring,
                     'sum':classmethod(sum),
                     'difference':classmethod(difference),
                     'product':classmethod(product),
                     'scalar_multiply':classmethod(scalar_multiply),
                     'AI':[(ring.AI, monoid.MI)],
                     'MI':[(ring.MI, monoid.MI)],
                     'generic':classmethod(generic)})


class MonoidAlgebra(Ring):

    """
    Elements of the monoid ring are of the form

    [(r1,m1), (r2,m2), (r3,m3)]

    where each r* is in the ring of the MonoidRing
    and each m* is in the monoid of the MonoidRing.

    """
    def __new__(cls, monoid, ring):

        assert 'ring' in ring.flags
        assert 'commutative_monoid' in ring.flags
        assert 'monoid' in monoid.flags

        cls.flags.add('commutative_monoid')

        return type("MonoidAlgebra", (Ring,),
                    {'monoid':monoid,
                     'ring':ring,
                     'sum':classmethod(sum),
                     'difference':classmethod(difference),
                     'product':classmethod(product),
                     'scalar_multiply':classmethod(scalar_multiply),
                     'AI':[(ring.AI, monoid.MI)],
                     'MI':[(ring.MI, monoid.MI)],
                     'generic':classmethod(generic)})

def collect_terms(ring, *terms):
    result = {}
    for ring_element, monoid_element in terms:
        if monoid_element in result:
            result[monoid_element] = ring.sum(result[monoid_element],
                                              ring_element)
        else:
            result[monoid_element] = ring_element
    return tuple([(r, m) for m, r in result.items()])

def simplify(ring, *terms):
    """
    Simplify terms, discarding those
    whose coefficients are I, the
    additive identity in the corresponding ring.

    >>> mr = MonoidRing(SubsetMonoid, Z)
    >>> simplify(mr, (3, 'a'), (mr.ring.inverse(3), 'a'))
    []
    >>>

    """
    return filter(lambda rm: rm[0] != ring.AI,
                  collect_terms(ring, *terms))

def delete(terms, terms_to_delete, match_coefs=True):
    """
    Remove the terms_to_delete from terms,
    via a set difference, if match_coefs=True.

    If match_coefs=False, then delete just based on the monoid elements
    of terms_to_delete.
    """

    if match_coefs:
        return list(set(terms).difference(terms_to_delete))
    else:
        monoids_to_keep = set([m for _, m in terms]).difference([m for _, m in terms_to_delete])
        return filter(lambda x: x[1] in monoids_to_keep, terms)


class EquivalenceRelation(object):

    def __init__(self, monoid, ring):
        self.monoid = monoid
        self.ring = ring
        self.monoid_ring = MonoidRing(monoid, ring)

    def equals(self, monoid_ring_element_a,
               monoid_ring_element_b):
        """
        Does a == b in the equivalence relation?
        The equivalence relation is

        set([m for _, m in monoid_ring_element_a]) ==
        set([m for _, m in monoid_ring_element_b])
        """
        return self.representative(a)  == self.representative(b)

    def generic(self, *monoid_elements):
        """
        Generate a generic representative
        of a sequence of monoid_elements.
        """
        return tuple([(self.ring.generic(),
                       self.monoid.validate(m)) for m in monoid_elements])

    def representative(self, monoid_ring_element):
        """
        Given an element in a monoid ring, i.e. a sequence of the
        form

        a = [(r1,m1), (r2, m2), ...]

        return a representation

        m = set([m1, m2, ...])

        A generic element of the equivalence class of a
        is

        [(self.ring.generic(), me) for me in m]

        """
        return set([m for _, m in monoid_ring_element])

    def sum(self, monoid_elements_a, monoid_elements_b):
        """
        Return the representative of a
        generic member of the equivalence class of
        a and a generic member of the equivalence
        class of b.
        """
        a = self.generic(*monoid_elements_a)
        b = self.generic(*monoid_elements_b)
        return self.representative(self.monoid_ring.sum(a,b))

    def product(self, monoid_elements_a, monoid_elements_b):
        """
        Return the representative of the sum of a 
        generic member of the equivalence class
        a  and a generic member of the equivalence
        class of b.
        """
        a = self.generic(*monoid_elements_a)
        b = self.generic(*monoid_elements_b)
        return self.representative(self.monoid_ring.product(a,b))

                            
def generic_monomials(sympy_monomials=True):
    """
    Return a sequence of generic monomials in 
    rings.MonoidRing(FactorMonoid, rings.SympyRingMon)

    The elements of the sequence have the form

    (sympy_monomial, factors)

    where sympy_monomial is a generic sympy monomial
    and factors is a sequence of generic Factors.

    Returns
    =======

    monomials : [(sympy_monomial, factors)]

    Rstr : a corresponding R formula string

    """
    result = []

    monoid_ring = rings.MonoidRing(FactorMonoid, rings.SympyRing)

    for _ in range(np.random.random_integers(0,5)):
        rings.SympyRing.monomial_generic = True
        symbol = rings.SympyRing.generic(sympy_monomials)
        for _ in range(np.random.random_integers(0,3)):
            factors = FactorMonoid.generic()
            result.append((symbol, factors))
    Rterms = []
    for s, fs in result:
        Rterms.append(':'.join([('I(%s)' % s)] + [f.name for f in fs]))
    
    return result, " + ".join(Rterms)
