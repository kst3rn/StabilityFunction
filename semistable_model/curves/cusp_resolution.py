r"""
Resolution of cusps on plane models of plane curves.
====================================================

EXAMPLES:

    sage: R.<z,x,y> = PolynomialRing(QQ, ("z","x","y"))
    sage: v_2 = QQ.valuation(2)
    sage: F = 2*z**4 + 2*x*y*z**2 - x**3*z + y**2*z**2 + x**4 + y**4

This form represents an integral plane model of a smooth quartic over `\QQ`,
whose special fiber with respect to the `2`-adic valuation has a cusp in
normal form in `(1:0:0)`. We can resolve this cusp as follows:

    sage: v_L, T, Fb = resolve_cusp(F, v_2)

Here `v_L` is the unique extension of `v_2` to a finite extension `L` of `K=\QQ`,

    sage: v_L.domain()
    Number Field in alpha with defining polynomial a^8 + 13464*a^7 + 21040*a^6 - 6832/5*a^5 + 93040*a^4 + 39360/43*a^3 + 2820928*a^2 - 10112/443*a - 10816/367

and `T` is a base change matrix with entries in `L` representing a plane model
whose special fiber consists of a semistable plane cubic in Weierstrass normal form,

    sage: Fb
    x^3 + z^2*y + z*y^2

plus the line at infinity. Thus this cubic is the one-tail of the semistable model
corresponding to the cusp.  

"""


from sage.all import QQ, NumberField, PolynomialRing, matrix, Infinity, randint, Curve, SR
from semistable_model.curves.approximate_factors import approximate_factorization


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

    EXAMPLES:

    The following example gives an error:

        sage: R.<x,y,z> = QQ[]
        sage: F = -16*x^4 - 15*x^3*y - 12*x^2*y^2 - 5*x*y^3 + 15*x^2*z^2 + 12*x*y*z^2 - 4*y^2*z^2 + 8*z^4
        sage: v_K = QQ.valuation(2)
        sage: X = PlaneCurveOverValuedField(F, v_K)
        sage: XX = X.semistable_model_with_rational_cusps()
        sage: Xs = XX.special_fiber()
        sage: C = Xs.rational_cusps()[0]
        sage: v_L = XX.base_ring_valuation()
        sage: L = v_L.domain()
        sage: T = C.move_to_e0_x2()
        sage: M = T.map_coefficients(v_L.lift, L)
        sage: cusp_model = XX.apply_matrix(M)
        sage: resolve_cusp(cusp_model.defining_polynomial(), v_L)
        ---------------------------------------------------------------------------
        NotImplementedError                       Traceback (most recent call last)
        ...
        NotImplementedError: Expected G[1] to be linear in beta (beta - r(alpha)).
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

    f, _, _ = assert_groebner_expected_order(G, a, b, c)
    f = f.univariate_polynomial()
    # proceed with root finding on f and set beta=r(alpha), gamma=s(alpha)

    # old code
    # assert len(G) == 3, "Unexpected Groebner basis length."
    # f = G[2].univariate_polynomial()

    f_factors = approximate_factorization(f, v_K)
    # print(f"f_factors = {f_factors}")
    # f is a univariate equation for alpha
    # we have to identify the correct factor of f, whose root alpha
    # leads to a solution with alpha, beta, gamma of positive valuation

    # beta and gamma are polynomials in alpha; we find the minimal valuation
    # m of the coefficients of these polynomials
    m1 = min(v_K(c) for c in G[0].coefficients())
    m2 = min(v_K(c) for c in G[1].coefficients())
    m = min(m1, m2)
    # print(m1, m2, m)
    
    # we have to compute alpha with a precision strictly greater than 
    # max(0, -m); because then we can already guarantee whether the valuation
    # of beta and gamma are positive
    
    for g in f_factors:
        # print(f"We try the factor g = {g.approximate_factor()} of degree {g.degree()}")
        # print()
        prec = max(0, -m) + 2
        # print(f"precision = {prec}")
        v_L, alpha, beta, gamma = _solve1(G, g, v_K, prec)
        # print(f"solution over {v_L.domain()} with precision {prec}")
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
    c, b, _ = G[0].parent().gens()
    K = v_K.domain()
    if E == 1:
        L = K.extension(g.approximate_factor(prec), "alpha")
        v_L = v_K.extension(L)
        alpha = L.gen()
    else:
        S = PolynomialRing(K, "pi")
        pi = S.gen()
        L = NumberField([g.approximate_factor(prec), pi**E-v_K.p()], ["alpha", "pi"])
        alpha, pi, *_ = L.gens()
        v_L = v_K.extension(L)
    beta = v_L.simplify(-G[1](c, b, alpha).univariate_polynomial()[0], prec)
    gamma = v_L.simplify(-G[0](c, b, alpha).univariate_polynomial()[0], prec)
    return v_L, alpha, beta, gamma 
        

def _valuation_matrix(F1, d, v_L):
    r"""
    Another helper function.
    """
    V = matrix(SR, d+1, d+1)
    t = Infinity
    for i in range(d+1):
        for j in range(d-i+1):
            if 2*i + 3*j < 6:
                V[i, j] = v_L(F1.coefficient([d-i-j, i, j]))
                s = V[i, j]/(6-2*i-3*j)
                if s < t:
                    t = QQ(s)
    return V, t


def _is_univariate_in(p, var):
    # True iff p involves no variables other than var
    vars_ = p.variables()
    return vars_ == () or vars_ == (var,)


def assert_groebner_expected_order(G, alpha, beta, gamma):
    """
    Assert that the Groebner basis G is exactly in the expected triangular form
        [gamma - s(alpha), beta - r(alpha), f(alpha)]
    (with respect to lex(gamma, beta, alpha)), and normalized so the leading
    coefficients of gamma and beta are 1.

    Raises NotImplementedError if not, otherwise returns (f, r, s).
    """
    if len(G) != 3:
        raise NotImplementedError(f"Expected Groebner basis of length 3, got {len(G)}.")

    p_gamma, p_beta, p_alpha = G[0], G[1], G[2]

    # --- check G[2] = f(alpha)
    if not _is_univariate_in(p_alpha, alpha):
        raise NotImplementedError("Expected G[2] to be univariate f(alpha).")

    # --- check G[1] = beta - r(alpha)
    if gamma in p_beta.variables():
        raise NotImplementedError("Expected G[1] to involve only alpha,beta (no gamma).")
    if p_beta.degree(beta) != 1:
        raise NotImplementedError("Expected G[1] to be linear in beta (beta - r(alpha)).")
    if p_beta.coefficient({beta: 1}) != 1:
        raise NotImplementedError("Expected G[1] to be monic in beta (leading coeff 1).")
    r = -p_beta.subs({beta: 0})
    if not _is_univariate_in(r, alpha):
        raise NotImplementedError("Expected r(alpha) to be univariate in alpha.")

    # --- check G[0] = gamma - s(alpha)
    if p_gamma.degree(gamma) != 1:
        raise NotImplementedError("Expected G[0] to be linear in gamma (gamma - s(alpha)).")
    if p_gamma.coefficient({gamma: 1}) != 1:
        raise NotImplementedError("Expected G[0] to be monic in gamma (leading coeff 1).")
    s = -p_gamma.subs({gamma: 0})
    if not _is_univariate_in(s, alpha):
        raise NotImplementedError("Expected s(alpha) to be univariate in alpha.")
    # (Strictly enforces "gamma - s(alpha)" i.e. no beta-dependence at all.)
    if beta in p_gamma.variables():
        raise NotImplementedError("Expected G[0] not to involve beta (gamma - s(alpha)).")

    f = p_alpha  # already univariate

    return f, r, s

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
