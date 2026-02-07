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

    sage: v_L, E, e = resolve_cusp(F, v_2)

Here `v_L` is the unique extension of `v_2` to a finite extension `L` of `K=\QQ`,

    sage: v_L.domain()
    Number Field in alpha with defining polynomial a^8 + 13464*a^7 + 21040*a^6 - 6832/5*a^5 + 93040*a^4 + 39360/43*a^3 + 2820928*a^2 - 10112/443*a - 10816/367

and is *tail*, a semistable plane cubic in Weierstrass normal form,

    sage: E
    Projective Plane Curve over Finite Field of size 2 defined by x^3 + y^2*z + y*z^2,

and `e` is a positive integer, the ramification index needed for an extension `L'/L`
over which the resolution can be performed.

TODO: 

- experiment with different variable orders in R to compare performance.


"""


from sage.all import ZZ, QQ, PolynomialRing, matrix, Infinity, randint, Curve, SR
from semistable_model.curves.approximate_solutions import approximate_solutions


def resolve_cusp(F, v_K, compute_matrix=False):
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

    If ``compute_matrix`` is ``False`` (the default):

    a triple `(v_L, E, e)`, where 
    - `e` is a positive integer,
    - `v_L` is an extension of `v_K` to a finite field extension `L/K` 
      such that the cusp can be resolved over any extension of `L`
      with ramification index `e`, 
    - `E` is a semistable plane cubic over the residue field of `L`, 
      the resulting one-tail. 

    If ``compute_matrix`` is ``True``:

    a tripel `(v_L, E, T)`, where v_L` and `E` are as above, with `e=1`,
    and `T` is an upper triangular`(3,3)`-matrix over `L`,
    representing the base change to the plane model resolving the cusp `P`.
   

    EXAMPLES:

    The following example caused problems before, but is now fixed:

        sage: R.<x,y,z> = QQ[]
        sage: F = -16*x^4 - 15*x^3*y - 12*x^2*y^2 - 5*x*y^3 + 15*x^2*z^2 + 12*x*y*z^2 - 4*y^2*z^2 + 8*z^4
        sage: v_K = QQ.valuation(2)
        sage: X = PlaneCurveOverValuedField(F, v_K)
        sage: XX = X.git_semistable_model_with_rational_cusps()
        sage: Xs = XX.special_fiber()
        sage: C = Xs.rational_cusps()[0]
        sage: v_L = XX.base_ring_valuation()
        sage: L = v_L.domain()
        sage: T = C.move_to_e0_x2()
        sage: M = T.map_coefficients(v_L.lift, L)
        sage: cusp_model = XX.apply_matrix(M)
        sage: v_L, E, e = resolve_cusp(cusp_model.defining_polynomial(), v_L) # long time!
        sage: v_L.domain()
        Number Field in alpha with defining polynomial a^8 - ... over its base field

        sage: E
        Projective Plane Curve over Finite Field of size 2 defined by x^3 + y^2*z + y*z^2

        sage: e
        1
        
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

    # we find *one* solution to A=B=C=0 in the maximal ideal of O_K
    s = approximate_solutions(J, v_K, positive_valuation=True, one_solution=True)
    if s is None:
        raise ValueError("No solution found! this shouldn't have happend..")
     
    L = s.extension()
    v_L = s.extension_valuation() 
    gamma, beta, alpha = s.approximation()
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
            gamma, beta, alpha = s.improve_approximation()
    # the ramification that we would need to resolve the cusp
    e = (t*v_L.value_group().denominator()).denominator()
        
    # compute the tail component
    Ab = PolynomialRing(v_L.residue_field(), 3, ["x", "y", "z"])
    Fb = Ab.zero()
    for i in range(d+1):
        for j in range(d-i+1):
            if 2*i + 3*j <= 6:
                c = _normalized_reduction(v_L, L(F1.coefficient([d-i-j, i, j])), t*(6-2*i-3*j))
                Fb += c * Ab.monomial(i, j, 3-i-j)
    E = Curve(Fb)

    if not compute_matrix:
        return v_L, E, e

    # we now have to extend L such that it contains an element Pi with
    # valuation t; then the reduction of F_2:=F_1(Pi^2*x,Pi^3*y,z)
    L = v_L.domain()
    if e != 1:
        # we have to replace L by an extension with ramification index e
        S = PolynomialRing(L, "pi")
        pi = S.gen()
        L = L.extension(pi**e - v_L.uniformizer(), "pi_e")
        v_L = v_L.extension(L)
        alpha = L(alpha)
        beta = L(beta)
        gamma = L(gamma)
    Pi = v_L.element_with_valuation(t)
    T = matrix(L, 3, 3, [1, alpha, beta, 0, Pi**2, Pi**2*gamma, 0, 0, Pi**3])
    return v_L, E, T


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

def _normalized_reduction(v_K, a, s):
    r""" Return the normalized reduction of an element wrt a discrete valuation.
    
    INPUT:

    - ``v_K`` -- a discrete valuation on a field `K`
    - ``a`` -- an element of `K`
    - ``s`` -- a rational number,

    It is assumed that

    .. MATH::

        v_K(a) \geq s.

    OUTPUT:

    If `v_K(a) = s` then the reduction of `\pi^{-s}\cdot a` is returned; here 
    `\pi^{-s}\in K` denotes the standard element of `K` with valuation `-s`, defined
    as 

    .. MATH::

        pi^{-s} := \pi_K^{s/v_K(\pi_K)},

    where `\pi_K` is a fixed uniformizer for `v_K`.

    If `v_K(a) > s`, then the zero element is returned. If `v_K(a) < s`, an error is
    raised.

    Note that the result depends on the the choice of `\pi_K`. We use the standard 
    generator returned by :meth:`uniformizer`. 
    
    """
    t = v_K(a)
    if t < s:
        raise ValueError("We must have `v_K(a) >= s; a = {a}, s = {s}")
    elif t == s:
        pi = v_K.uniformizer()
        m = ZZ(t/v_K(pi))
        return v_K.reduce(pi**(-m)*a)
    else:
        return v_K.residue_field().zero() 


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
