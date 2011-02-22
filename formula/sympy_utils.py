import sympy

def getparams(expression):
    """ Return the parameters of an expression that are not Term 
    instances but are instances of sympy.Symbol.

    Examples
    --------
    >>> x, y, z = [Term(l) for l in 'xyz']
    >>> f = Formula([x,y,z])
    >>> getparams(f)
    []
    >>> f.mean
    _b0*x + _b1*y + _b2*z
    >>> getparams(f.mean)
    [_b0, _b1, _b2]
    >>>                 
    >>> th = sympy.Symbol('theta')
    >>> f.mean*sympy.exp(th)
    (_b0*x + _b1*y + _b2*z)*exp(theta)
    >>> getparams(f.mean*sympy.exp(th))
    [theta, _b0, _b1, _b2]
    """
    atoms = set([])
    expression = np.array(expression)
    if expression.shape == ():
        expression = expression.reshape((1,))
    if expression.ndim > 1:
        expression = expression.reshape((np.product(expression.shape),))
    for term in expression:
        atoms = atoms.union(sympy.sympify(term).atoms())
    params = []
    for atom in atoms:
        if isinstance(atom, sympy.Symbol) and not is_term(atom):
            params.append(atom)
    params.sort()
    return params


def getterms(expression):
    """ Return the all instances of Term in an expression.

    Examples
    --------
    >>> x, y, z = [Term(l) for l in 'xyz']
    >>> f = Formula([x,y,z])
    >>> getterms(f)
    [x, y, z]
    >>> getterms(f.mean)
    [x, y, z]
    """
    atoms = set([])
    expression = np.array(expression)
    if expression.shape == ():
        expression = expression.reshape((1,))
    if expression.ndim > 1:
        expression = expression.reshape((np.product(expression.shape),))
    for e in expression:
        atoms = atoms.union(e.atoms())
    terms = []
    for atom in atoms:
        if is_term(atom):
            terms.append(atom)
    terms.sort()
    return terms


