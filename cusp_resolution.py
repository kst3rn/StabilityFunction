r"""
Resolution of cusps on plane models of plane curves.
====================================================

EXAMPLES:

    sage: from cusp_resolution import resolve_cusp
    sage: R.<z,x,y> = PolynomialRing(QQ, ("z","x","y"))
    sage: v_2 = QQ.valuation(2)
    sage: F = 2*z**4 + 2*x*y*z**2 - x**3*z + y**2*z**2 + x**4 + y**4
    sage: resolve_cusp(F, v_2)

"""


from sage.all import QQ, PolynomialRing, matrix, Infinity, randint, Curve
from approximate_factors import approximate_factorization


def resolve_cusp(F, v_K):
    r""" Return a base change matrix resolving the cusp.

    INPUT:

    - ``F`` -- a trivariat form of degree `\geq 3` over a field `K`
    - ``v_K`` -- a nontrivial discrete valuation on `K`

    It is assumed that 
    - `F` defines a smooth quartic curve `X`,
    - `F` is integral and primitive with respect to `v_K`, so that it defines
      an integral model `\mathcal{X}` of `X`, and
    - the special fiber has a cusp in normal form at the point `P=(1:0:0)`.

    OUTPUT:

    a tripel `(v_L, T, \bar{F})`, where `v_L` is an extension of `v_K` to a finite field
    extension `L/K`, `T` is an upper triangular`(3,3)`-matrix over `L`,
    representing the base change to the plane model resolving the cusp `P`, and
    `\bar{F}` is a semistable cubic, the resulting one-tail.

    """
    # check validity of the input
    assert v_K.is_discrete_valuation()
    K = v_K.domain()
    F = F.change_ring(K)
    # Check that F is a polynomial in three variables
    assert F.nvariables() == 3, "F must be a polynomial in three variables."
    # Check if F is homogeneous by inspecting all monomials, and compute degree
    monomial_degrees = [m.total_degree() for m in F.monomials()]
    assert len(set(monomial_degrees)) == 1, "F must be homogeneous."
    d = monomial_degrees[0]
    assert d >= 3, "F must be of degree at least 3."
    # Check that F is integral and satisfies the cusp condition at (1:0:0)
    for m in F.monomials():
        assert v_K(F.monomial_coefficient(m)) >= 0, "F must be integral."
        k, i, j = m.exponents()[0]
        if 2*i + 3*j < 6:
            assert v_K(F.monomial_coefficient(m)) > 0, "P=(1:0:0) is not a cusp"
        elif 2*i + 3*j == 6:
            assert v_K(F.monomial_coefficient(m)) == 0, "P=(1:0:0) is not a cusp"

    # Initialization
    k = v_K.residue_field()
    p = k.characteristic()
    R = PolynomialRing(K, ("c", "b", "a"), order='lex')
    c, b, a = R.gens()
    F0 = F.change_ring(R)
    z0, x0, y0 = F0.variables()
    F0 = F0(z0, x0 + a*z0, y0 + b*z0 + c*x0)

    # the system of equations in a,b,c we want to (approximately) solve
    if p == 2:
        A = F0[4, 0, 0]
        B = F0[3, 1, 0]
        C = F0[2, 2, 0]
    else:
        A = F0[4, 0, 0]
        B = F0[3, 0, 1]
        C = F0[2, 1, 1]
    # we want to find an approximate solution alpha, beta, gamma for
    # A=B=C=0, with valuations at least v_a,v_b,v_c
    J = R.ideal([A, B, C])
    G = J.groebner_basis()
    assert len(G) == 3, "Unexpected Groebner basis length."

    f = G[2].univariate_polynomial()
    f_factors = approximate_factorization(f, v_K)
    # f is a univariate equation for alpha
    # we have to identify the correct factor of f, whose root alpha
    # leads to a solution with alpha, beta, gamma of positive valuation

    # We need a better and more failsafe way to increase the precision. 
    # There are two precision parameters to play with:
    # 1. the precision of the approximate factors of f
    # 2. the minimal required valuation of H(gamma,beta,alpha), which is
    #    a measure how good alpha, beta, gamma solve the equations.
    # We can arbitrarily increase the first parameter, and then the second
    # will increase as well (not clear how fast).
    # The problem is that we do not know beforehand which factor gives a
    # valid solution; it seems that typically it is only one of two or three,
    # no matter how much we increase the precision.
    # A necessary condition for an approximate solution to be "valid" is that
    # the valuations of alpha, beta, gamma are positive; a sufficient condion
    # is that it is *stably positive*, i.e. remains positive when improving
    #  the precision. It would be good if we had a test for this! 
    
    for g in f_factors:
        # print(f"We try the factor g = {g.approximate_factor()} of degree {g.degree()}")
        # print()
        # We set an initial precision for such a solution which should be
        # sufficient for testing this
        # THIS IS EXPERIMENTAL AND HAS TO BE REPLACED BY A FAILSAFE METHOD LATER!
        prec = 5
        N = 25  # N=20 was insufficient in one example
        while True:
            v_L, alpha, beta, gamma = _solve1(G, g, v_K, prec)
            # print(f"solution over {v_L.domain()} with precision {prec}")
            M = min(v_L(H(gamma, beta, alpha)) for H in G)
            # print(f"M = {M}")
            # print()
            if M >= N:
                # alpha, beta, gamma are solutions up to the desired precision
                # we exit the while loop
                break
            else:
                prec += 5
        V = [v_L(a) for a in [alpha, beta, gamma]]
        # print(f"V = {V}")
        if all(v > 0 for v in V):
            # we have found the correct factor g
            break
    else:
        # no factor was found
        raise ValueError("no solution was found!")
    
    # now g is the correct factor of f, but its current precision may not
    # be sufficient
    z, x, y = F.variables()
    while True:
        F1 = F(z, alpha*z + x, beta*z + gamma*x + y)
        # we compute the matrix V with entries v_L(coef of x^i*y^j)
        # and the "minimal slope" t to check if F1 represents,
        # after scaling, a resolution of the cusp
        V, t = _valuation_matrix(F1, d, v_L)
        if p == 2:
            tests = [(0,0), (1,0), (2,0)]
        else:
            tests = [(0,0), (0,1), (1,1)]
        if all([V[i, j]/(6-2*i-3*j) > t for i, j in tests]):
            break
        else:
            prec += 5
            v_L, alpha, beta, gamma = _solve1(G, g, v_K, prec)

    # we now have to extend L such that it contains an element Pi with
    # valuation t; then the reduction of F_2:=F_1(Pi^2*x,Pi^3*y,z)
    
    if not t in v_L.value_group():
        e = (t*v_L.value_group().denominator()).denominator()
        # we have to replace L by an extension with ramification index e,
        # and compute alpha, beta, gamma again
        v_L, alpha, beta, gamma = _solve1(G, g, v_K, prec, E=e)
    L = v_L.domain()
    Pi = v_L.element_with_valuation(t)
    F1 = F(z, alpha*z + x, beta*z + gamma*x + y)
    F2 = F1(z, Pi**2*x, Pi**3*y)/Pi**6
    F2b = F2.map_coefficients(v_L.reduce, v_L.residue_field())
    Fb = [Gb for Gb, _ in F2b.factor() if Gb.degree() == 3][0]
    T = matrix(L, 3, 3, [1, alpha, beta, 0, Pi**2, Pi**2*gamma, 0, 0, Pi**3])
    return v_L, T, Fb
    

def _solve1(G, g, v_K, prec, E=1):
    r""" Return v_L, alpha, beta, gamma

    This is a helper function for `resolve_cusp`.
    
    """
    b, c, _ = G[0].parent().gens()
    K = v_K.domain()
    if E == 1:
        L = K.extension(g.approximate_factor(prec), "alpha")
        v_L = v_K.extension(L)
        alpha = L.gen()
    else:
        S = PolynomialRing(K, "pi")
        pi = S.gen()
        L = K.extension([g.approximate_factor(prec), pi**E-v_K.p()], ["alpha", "pi"])
        alpha, pi = L.gens()
        v_L = v_K.extension(L)
    beta = v_L.simplify(-G[1](c, b, alpha).univariate_polynomial()[0], prec)
    gamma = v_L.simplify(-G[0](c, b, alpha).univariate_polynomial()[0], prec)
    return v_L, alpha, beta, gamma 
        

def _valuation_matrix(F1, d, v_L):
    V = matrix(QQ, d+1, d+1)
    t = Infinity
    for i in range(d+1):
        for j in range(d-i+1):
            if 2*i + 3*j < 6:
                V[i, j] = v_L(F1.coefficient([d-i-j, i, j]))
                s = V[i, j]/(6-2*i-3*j)
                if s < t:
                    t = s
    return V, t

# ------------------------------------------------------------------------------------------------------

#                    Tests

def random_padic_integer(v_K):
    r""" Return a random element of the ring of integers.
    
    INPUT:

    - ``v_K`` -- a p-adic valuation on a number field `K`

    OUTPUT:

    a random element of the ring of integers of `v_K`.

    """
    K = v_K.domain()
    pi = v_K.uniformizer()
    t = v_K(pi)
    v_K = v_K/t
    a = K.random_element()
    if a == 0:
        return a
    m = randint(0, 2)
    return a*pi**(-v_K(a) + m)


def random_curve_with_cusp(v_K, d=4):
    r""" Return a random plane curve with a cusp.
    
    INPUT:

    - ``v_K`` -- a p-adic valuation on a number field `K`
    - ``d`` -- an integer `\geq 3`

    OUTPUT:

    A form `F` of degree `d` over `K` in `z,x,y,` which represents
    an plane integeral model of a smooth plane curve over `K`. The special
    fiber has a cusp at `(1:0:0)` in normal form, i.e. with leading term
    `y^2-x^3`.

    """
    K = v_K.domain()
    R = PolynomialRing(K, ("z", "x", "y"))
    z, x, y = R.gens()
    pi = v_K.uniformizer()
    while True:
        F = R.zero()
        for i in range(d + 1):
            for j in range(d - i + 1):
                if (i,j) == (0,2):
                    F += y**2*z**(d-2)
                elif (i,j) == (3,0):
                    F += -x**3*z**(d-3)
                elif 2*i+3*j < 6:
                    F += pi*random_padic_integer(v_K)*x**i*y**j*z**(d-i-j)
                else:
                    F += random_padic_integer(v_K)*x**i*y**j*z**(d-i-j)
        X = Curve(F)
        if X.is_smooth():
            return F


def test_suite(v_K, N):
    for _ in range(N):
        F = random_curve_with_cusp(v_K)
        print(f"F = {F}")
        v_L, _, Fb = resolve_cusp(F, v_K)
        print(f" L = {v_L.domain()}")
        print(f"Fb = {Fb}")
        print()



# ---------------------------------------------------------------------------------

#               problematic examples

R = PolynomialRing(QQ, ("z","x","y"))
z, x, y = R.gens()
F = -2*z**4 + 2*z**3*x - 4/83*z**2*x**2 - z*x**3 + 2*z**3*y + 4/3*z*x**2*y - 2*x**3*y + z**2*y**2 + 2*x**2*y**2 - y**4

F = -40*z**3*x - z*x**3 - 21/109*x**4 + 8/13*z**3*y + 8*z**2*x*y + 2*z*x**2*y + 2*x**3*y + z**2*y**2 + z*x*y**2 - 4/7*x**2*y**2 - 2/11*z*y**3 + 2*x*y**3
