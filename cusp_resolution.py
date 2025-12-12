r"""
Resolution of cusps on plane models of plane curves.
====================================================

EXAMPLES:

    sage: from cusp_resolution import resolve_cusp
    sage: R.<z,x,y> = PolynomialRing(QQ, ("z","x","y"))
    sage: v_2 = QQ.valuation(2)
    sage: F = 2*z^4 + 2*x*y*z^2 - x^3*z + y^2*z^2 + x^4 + y^4
    sage: resolve_cusp(F, v_2)

"""


from sage.all import QQ, PolynomialRing, matrix, Infinity, lcm

import sys
sys.path.insert(0, '/home/stefan/code/MCLF/mclf')   # parent of the package dir, not .../mclf
import mclf
from mclf.padic_extensions.padic_number_fields import pAdicNumberField
from mclf.padic_extensions.approximate_factorizations import approximate_factorization


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

    a pair `(v_L, T)`, where `v_L` is an extension of `v_K` to a finite field
    extension `L/K` and `T` is an upper triangular`(3,3)`-matrix over `L`,
    representing the base change to the plane model resolving the cusp `P`.

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

    if p == 2:
        A = F0[4, 0, 0]
        B = F0[3, 1, 0]
        C = F0[2, 2, 0]
    else:
        A = F0[4, 0, 0]
        B = F0[3, 0, 1]
        C = F0[2, 1, 1]
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C ={C}")
    print()
    J = R.ideal([A, B, C])
    G = J.groebner_basis()
    assert len(G) == 3, "Unexpected Groebner basis length."

    prec = 5
    while True:
        v_L, alpha, beta, gamma = approximate_solution(G, v_K, prec)
        print(f"alpha = {alpha}, beta={beta}, gamma={gamma}")
        print(f"prec = {prec}")
        print(f"v_L(A)={v_L(A(gamma,beta,alpha))}")
        print(f"v_L(B)={v_L(B(gamma,beta,alpha))}")
        print(f"v_L(C)={v_L(C(gamma,beta,alpha))}")
        print()
        L = v_L.domain()
        z, x, y = F.variables()
        F1 = F(z, alpha*z + x, beta*z + gamma*x + y)
        print(f"F1 = {F1}")
        V = matrix(QQ, d+1, d+1)
        t = Infinity
        for i in range(d+1):
            for j in range(d-i+1):
                if 2*i + 3*j < 6:
                    V[i, j] = v_L(F1.coefficient([d-i-j, i, j]))
                    s = V[i, j]/(6-2*i-3*j)
                    if s < t:
                        t = s
        print(f"V = {V}")
        if p == 2:
            tests = [(0,0), (1,0), (2,0)]
        else:
            tests = [(0,0), (0,1), (1,1)]
        if all([V[i, j]/(6-2*i-3*j) > t for i, j in tests]):
            break
        else:
            prec += 5
    # this is not what should be returned. L is in general not yet the correct extension;
    # we may need to add ramification to have t in the value group. Then we also have to 
    # do another transformation.
    return v_L, F1, t


def approximate_solution(G, v_K, prec):
    r""" Approximate solution of the system defined by G.

    INPUT:

    - ``G`` -- a list of three polynomials in `K[a,b,c]`
    - ``v_K`` -- a discrete valuation on `K`
    - ``prec`` -- desired precision

    We assume that `G=[g_1, g_2, g_3]` is a Groebner basis of an ideal
    defining a zero-dimensional scheme over `K`, with respect to the lexicographic
    order with `c > b > a`. We further assume that `g_3` is monic and univariate
    in `a`, that `g_2` is monic of degree one in `b` with coefficients in `K[a]`,
    and that `g_1` is monic of degree one in `c` with coefficients in `K[a]`.

    OUTPUT:

    a tuple `(v_L, \alpha, \beta, \gamma)`, where `v_L` is an extension of `v_K`
    to a finite field extension `L/K` and `\alpha, \beta, \gamma` are elements
    of `L` with positive valuation, satisfying the system defined by `G` up to
    precision `prec`.

    """
    K = v_K.domain()
    Kh = pAdicNumberField(K, v_K)
    c, b, _ = G[0].parent().gens()
    f = G[2].univariate_polynomial()
    f = lcm(a.denominator() for a in f.coefficients()) * f
    F = approximate_factorization(Kh, f)
    print(f"approximate factorization of f={f}:")
    print(F)
    print()

    # we have to choose the *right* factor of f, so that the solutions
    # alpha, beta, gamma have positive valuations. At the moment I don't
    # have a good way to choose the right factor at this point. Therefore,
    # we try all factors until a good one is found. 
    for g in F:
        print(f"We try the factor g={g.approximate_polynomial()}")
        # L is the stem field of g
        L = Kh.simple_extension(g.approximate_polynomial())
        v_L = L.valuation()
        # we find an approximate root of g over L
        gL = g.base_change(L, ignore_linear_factors=False)
        for h in gL:
            if h.degree() == 1:
                alpha = -h.approximate_polynomial()[0]
                break
        else:
            raise ValueError("something is wrong: no root found")
        # now alpha is an approximate root of g over L
        # beta and gamma are polynomials in alpha
        # if alpha, beta, gamma is not a sufficiently good solution,
        # we have to improve the precision of alpha and repeat
        N = prec
        while True:
            print(f"We try to find alpha, beta, gamma with precision {N}." )
            beta = L.approximation(- G[1](c, b, alpha).univariate_polynomial()[0], N)
            gamma = L.approximation(- G[0](c, b, alpha).univariate_polynomial()[0], N)
            if all([v_L(G[i](gamma, beta, alpha)) > prec for i in range(3)]):
                print(f"We found a solution with precision {N}")
                # ok, alpha, beta, gamma are solutions up to the desired precision
                # but we stil have to check whether they all have positive valuation
                if v_L(alpha) > 0 and v_L(beta) > 0 and v_L(gamma) > 0:
                    return v_L, alpha, beta, gamma
                else:
                    # we assume that we got the wrong factor g of f
                    print(f"But v_L(alpha)={v_L(alpha)}, v_L(beta)={v_L(beta)}, v_L(gamma)={v_L(gamma)}")
                    print
                    break
            else:
                # we improve the precision of alpha
                N += 5
                # this is bad, because ``approximate_root`` only works (for the moment)
                # if f has coefficients in ZZ
                alpha = L.approximate_root(f, alpha, N)
    # if we get here, we are unlucky
    raise ValueError("No solution found!")
        