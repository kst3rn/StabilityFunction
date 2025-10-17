r"""
Resolution of cusps on plane models of plane curves.
====================================================

EXAMPLES:

    sage: from cusp_resolution import resolve_cusps
    sage: R.<z,x,y> = PolynomialRing(QQ, ("z","x","y"))
    sage: v_2 = QQ.valuation(2)
    sage: F = 2*z^4 + 2*x*y*z - x^3*z + y^2*z^2 + x^4 + y^4
    sage: resolve_cusp(F, v_2)

"""


from sage.all import PolynomialRing, Matrix, vector

import sys
sys.path.insert(0, '/home/stefan/code/mclf_alt')   # parent of the package dir, not .../mclf
import mclf
from mclf.padic_extensions.padic_number_fields import pAdicNumberField
from mclf.padic_extensions.approximate_factorizations import approximate_factorization


def resolve_cusp(F, v_K):
    r""" Return a base change matrix resolving the cusp.

    INPUT:

    - ``F`` -- a trivariat form of degree `3` over a field `K`
    - ``v_K`` -- a nontrivial discrete valuation on `K`

    It is assumed that `F` that
    - `F` defines a smmoth quartic curve `X`,
    - `F` is integral and primitive with respect to `v_K`, so that it defines
      an integral model `\mathcal{X}` of `X`, and
    - the special fiber has a cusp in normal form at the point `P=(1:0:0)`.

    OUTPUT:

    a pair `(v_L, T)`, where `v_L` is an extension of `v_K` to a finite field
    extension `L/K` and `T` is an upper triangular`(3,3)`-matrix over `L`,
    representing the base change to the plane model resolving the cusp `P`.

    """
    # check validity of the input (to do)
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
    S = PolynomialRing(R, ("x", "y"))
    x0, y0 = S.gens()
    f = F(K.one(), x0 + a, y0 + b + c*x0)

    if p == 2:
        A = f[0, 0]
        B = f[1, 0]
        C = f[2, 0]
    else:
        A = f[0, 0]
        B = f[0, 1]
        C = f[1, 1]
    J = R.ideal([A, B, C])
    G = J.groebner_basis()
    assert len(G) == 3, "Unexpected Groebner basis length."
    return G

    prec = 5
    while True:
        v_L, alpha, beta, gamma = approximate_solution(G, v_K, prec)
        L = v_L.domain()
        T = Matrix(L, 3, 3, [[1, alpha, beta], [0, 1, gamma], [0, 0, 1]])
        z, x, y = F.variables()
        F1 = F(T * vector((z, x, y)))

        raise NotImplementedError
