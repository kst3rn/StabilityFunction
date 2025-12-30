r""""
Parametric optimization
=======================

Let `K` be a field with a discrete valuation `v_K` and `S:=K[t]` the polynomial
ring over `K`. Let `F\in S[x_0,...,x_n]` be a homogeneous polynomial of degree
`d` in `n+1` variables. Given a finite extension `L/K` and and element `a\in L`,
we obtain a homogenous form `F_a\in L[x]` by substituting `a` for `t`. We can
then compute the minimum of the stability function `m_{F_a}` on the apartment
corresponding to `F_a`.

We can extend `a\mapsto m_{F_a}` to a function `v\mapsto m_v` on the analytic
space `(\AA_K^1)^{an}` such that `a` corresponds to the pseudovaluation
`v(f):=v_L(f(a))`. The goal of this script is to express this function in the
form

.. MATH::

    m_v = \min(h_1(v),\ldots,h_r(v)),

where each `h_i` is an explicit valuative function of the form

.. MATH::

    h_i(v) = a_0 + a_1 v(f_1) + \ldots + a_n v(f_n),

where `f_1,\ldots,f_n\in S` are irreducible ponomials in `t` and
`a_0,\ldots,a_n` are positive rational numbers.


EXAMPLES::

    sage: from parametric_optimization import minimum_as_valuative_function
    sage: S.<t> = QQ[]
    sage: v0 = GaussValuation(S, QQ.valuation(2))
    sage: R.<x0,x1,x2> = S[]
    sage: F = 1352*x2^3-x0^3-78*x0*x2^2+x1^2*x2
    sage: F1 = F(t^2*x0, t^3*x1, x2)*t^2/4
    sage: F2 = F1(x0+x2, x0+x1, x2)
    sage: h, e = minimum_as_valuative_function(F2, v0)
    sage: h
    min(-12 + 96*v(t), -18 + 96*v(t), -24 + 84*v(t) + 3*v(t^4 + 26), -24 + 84*v(t) + 2*v(t^6 + 78*t^2 - 1352))

    sage: e
    12

    sage: h.local_maxima()
    [(63, Point of type II on Berkovich line, corresponding to v(t) >= 1),
     (+Infinity, Point of type I on Berkovich line given by t = 0),
     (3, Point of type II on Berkovich line, corresponding to v(t^4 - 2) >= 5/2),
     (2, Point of type II on Berkovich line, corresponding to v(t^4 + 2) >= 2)]

"""

from sage.all import ZZ, QQ, matrix, vector, lcm, prod
from sage.geometry.polyhedron.constructor import Polyhedron
from semistable_model.stability import MinimumOfValuativeFunctions, valuative_function


def minimum_as_valuative_function(F, v0):
    r"""
    Compute the minimum of the stability function of F_a on the apartment
    corresponding to F_a.

    INPUT:

    - ``F`` -- a homogeneous polynomial in n+1 variables over a polynomial
                ring S = K[t], `where `K` is a field with a discrete valuation
    - ``v0`` -- a MacLane pseudovaluation on S.

    OUTPUT:

    A pair `(h, e)` where `h` is an admissible function defined on the discoid
    with root `v_0` and `e` is a positive integer.

    To explain what this means, write `F = \sum_{i\in I} a_i x^i`, where `I`
    is the set of exponents of `F`, `x=(x_0,\ldots,x_n)`, and the `a_i` are
    polynomials in `t`.

    Let us fix a MacLane pseudovaluation `v\geq v_0` on `S`. The
    *stability function* for `F` on the apartment corresponding to
    `x=(x_0,...,x_n)` (w.r.t. `v`) is the function which associated
    to `w=(w_0,...,w_n)` the maximum

    .. MATH::

        \phi_F(w) = \max_{i\in I} l_i(w) - v(a_i),

    where

    .. MATH::

        l_i(w) = \sum_{j=0}^n (d/(n+1) - i_j) w_j.

    If `F` defines a stable hypersurface (which we assume) then this function
    assumes a minimum `m_v`.

    The goal of this function is to express this function in the
    form `m_v=-h(v)/e`, with

    .. MATH::

        h(v) = \min(h_1(v),\ldots,h_r(v)),

    where each `h_i` is an explicit valuative function of the form

    .. MATH::

        h_i(v) = a_0 + a_1 v(f_1) + \ldots + a_n v(f_n),

    and `f_1,\ldots,f_n\in S` are irreducible ponomials in `t` and
    `a_0,\ldots,a_n` are positive rational numbers.



    ALGORITHM:

    Let `I` be the set of exponents of `F`. Let `a_i` denote the coefficient
    of the monomial `x^i` in `F`, a polynomial in `t`.

    tbc
    """
    II = F.exponents()
    monomials = [F.parent().monomial(*i) for i in II]
    a = [F.monomial_coefficient(m) for m in monomials]
    A = matrix_from_I(II)
    E = eliminate_w(A)
    e = sum(E[0])
    h = []
    for k in range(E.nrows()):
        f_k = prod(a_i**E[k, i] for i, a_i in enumerate(a))
        h_k = valuative_function(f_k, v0)
        h.append(h_k)
    return MinimumOfValuativeFunctions(h), e


# -----------------------------------------------------------------------------

#               helper functions


def matrix_from_I(II):
    """
    Return the matrix corresponding to I.

    INPUT:

    - ``I`` -- a list of triples `(i_0,i_1,i_2)` of nonegative integers
                such that `i_0 + i_1 + i_2 = d`, for a fixed `d`.

    OUTPUT:

    a matrix with `|I|` rows and 3 columns, where the `j`-th row is the linear
    form `l_j(w)` corresponding to the triple `I[j]`.
    """
    d = sum(II[0])
    N = len(II)
    A = matrix(QQ, N, 3)
    for j in range(N):
        i = II[j]
        A[j] = vector([QQ(d)/3 - i[0], QQ(d)/3 - i[1], QQ(d)/3 - i[2]])
    return A


##############################################################################
# Fourier–Motzkin for a single variable
##############################################################################


def fm_eliminate_one_variable(ineqs, var_index):
    r"""
    Eliminate the variable with index = var_index from the inequality system.

    INPUT:
      - ineqs: a list of inequalities in "cdd format":
          [b, a_0, ..., a_{p-1}] means  b + sum_j a_j v_j >= 0
        All must have the same length p+1, for some p >= 1.

      - var_index: which variable to eliminate (0-based)

    OUTPUT:
      - A new list of inequalities (again in cdd format, but for p-1 variables)
        in which `var_index` is removed.
        They describe exactly the projection of the original region onto the
        remaining p-1 variables.

    EXPLANATION:
      The function splits the inequalities into three sets:
        S_plus, S_minus, S_zero
      and forms all pairwise combinations from S_plus x S_minus to cancel out
      the chosen variable. S_zero's inequalities are carried over unchanged
      (except for dropping the eliminated coordinate from the format).
    """
    if not ineqs:
        # Nothing to eliminate
        return []
    # The total number of variables
    num_vars = len(ineqs[0]) - 1  # if each ineq is [b, a_0..a_{p-1}], p = len(ineq) - 1

    # Partition into S_plus, S_minus, S_zero
    S_plus = []
    S_minus = []
    S_zero = []
    for ineq in ineqs:
        coeff = ineq[1 + var_index]  # coefficient of the variable we are eliminating
        if coeff > 0:
            S_plus.append(ineq)
        elif coeff < 0:
            S_minus.append(ineq)
        else:
            S_zero.append(ineq)

    new_ineqs = []

    # 1) Keep the S_zero inequalities (dropping the eliminated variable's coefficient)
    for ineq in S_zero:
        new_ineq = ineq[:1]  # the constant term b
        for j in range(num_vars):
            if j != var_index:
                new_ineq.append(ineq[1 + j])
        new_ineqs.append(new_ineq)

    # 2) For each pair (p in S_plus, n in S_minus), form the linear combination
    #    that cancels the variable var_index
    for p_ineq in S_plus:
        cp = p_ineq[1 + var_index]
        for n_ineq in S_minus:
            cn = n_ineq[1 + var_index]
            # We want: alpha = cp, beta = -cn
            # so alpha*(n_ineq) + beta*(p_ineq) => coefficient_of_var_index = 0.
            alpha = cp
            beta = -cn
            combined = [alpha*n_ineq[k] + beta*p_ineq[k] for k in range(num_vars+1)]
            # Remove var_index from combined
            new_ineq = [combined[0]]  # new b
            for j in range(num_vars):
                if j != var_index:
                    new_ineq.append(combined[1 + j])

            new_ineqs.append(new_ineq)

    return new_ineqs

##############################################################################
# Main function to eliminate x from A*x <= c
##############################################################################


def eliminate_w(A):
    r""" Given an `m` x `(n+1)` matrix `A` over `\QQ`,
    eliminate `w` from the system `A\cdot w \leq t`.

    Here `w=(w_0, ..., w_{n-1})` and `t=(t_0, ..., t_{m-1})`.

    INPUT:

    - `A` -- a matrix over `\QQ`.

    OUTPUT:

    A matrix `E=(e_{k,j})` with `m` columns and nonnegative integer entries
    `e_{k,j}\geq 0` such that `E\cdot t >= 0` if and only if there exists an
    `w` such that `A\cdot w <= t`. Moreover, the sum of the entries in each
    row of `E` is the same positive integer

    .. MATH::

        e := \sum_{j=0}^{m-1} e_{k,j}, \quad\forall k.

    See Lemma 5.2 of [SternWewers]_.

    ALGORITHM:

    We use the Fourier-Motzkin elimination algorithm to eliminate the variables
    `w_{n-1}, w_{n-2}, ..., w_0` in that order. The algorithm proceeds in two
    steps:

    1) Build the initial inequalities in cdd format for `A\cdot w \leq t`:
         row `i`:  `- (A[i,0]*w_0 + ... + A[i,n-1]*w_{n-1}) + t_i >= 0`

    2) Eliminate `w_{n-1}, w_{n-2}, ..., w_0` in that order so indices remain
         stable for each step (After removing the highest index, the lower ones
         don't shift.)

    """
    m = A.nrows()
    n = A.ncols()

    # 1) Build the initial inequalities in cdd format for A*x <= c:
    #    row i:  - (A[i,0]*x_0 + ... + A[i,n-1]*x_{n-1}) + c_i >= 0
    ineqs = []
    # total number of variables = n + m
    # cdd format length = 1 + (n+m) = b plus each variable's coefficient
    for i in range(m):
        row_i = [QQ(0)]*(1 + n + m)  # [b, a_0, ..., a_{n+m-1}]
        # b=0
        row_i[0] = 0
        # x-part (the first n)
        for j in range(n):
            row_i[1+j] = -A[i, j]
        # c-part (the last m)
        row_i[1 + n + i] = 1
        ineqs.append(row_i)

    # 2) Eliminate x_{n-1}, x_{n-2}, ..., x_0 in that order
    #    so indices remain stable for each step
    #    (After removing the highest index, the lower ones don't shift.)
    for var_index in reversed(range(n)):
        ineqs = fm_eliminate_one_variable(ineqs, var_index)

    # After that, only t_0,...,t_{m-1} remain.
    # Before returning the result, we simplify the inequalities
    # by removing redundant ones. we use the build -in polyhedron
    # class to do this.
    P = Polyhedron(ieqs=ineqs)
    ineqs = [ineq[1:] for ineq in P.inequalities_list()]  # remove the constant term zero
    # Multiply each inequality with a positive integer such that all entries are integers
    # and the sum of the entries is equal for each inequality
    E = []
    for ineq in ineqs:
        denom_lcm = lcm([a.denominator() for a in ineq])
        new_ineq = [a * denom_lcm for a in ineq]
        E.append(new_ineq)

    # Compute the sums for each line
    sums = [sum(ineq) for ineq in E]
    # Compute the lcm of these sums
    lcm_sums = lcm(sums)
    # Multiply each line with the appropriate constant to make the sums equal to lcm_sums
    for i, e_i in enumerate(sums):
        factor = lcm_sums / e_i
        E[i] = [a * factor for a in E[i]]

    assert len(E) > 0, f"E is empty! A = {A}"
    return matrix(ZZ, E)
