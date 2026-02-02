# -*- coding: utf-8 -*-
r"""
Approximate solutions for zero-dimensional ideals
=================================================

Let `K` be a number field and `v_K` a `p`-adic valuation on `K`. Let

.. MATH::

    J = (F_1,\ldots,F_m) \subset K[x_1,\ldots,x_n]

be an ideal of dimension zero. Then the system

.. MATH::

    F_1(x_1,\ldots,x_n)=\cdots=F_m(x_1,\ldots,x_n)=0

has only finitely many solutions in an algebraic closure of `K`.

We are interested in *approximate solutions* with respect to `v_K`. By this we
mean a sequence

.. MATH::

    a^{(k)} = (a_1^{(k)},\ldots,a_n^{(k)}) \in L^n

for a finite extension of valued fields `(L,v_L)/(K,v_K)` such that

.. MATH::

    v_L(F_i(a^{(k)})) \to \infty \quad (k\to\infty)

for all generators `F_i` (equivalently, such that the residual valuations can
be made arbitrarily large).

In this module we currently restrict to *integral* approximate solutions, i.e.
those with

.. MATH::

    v_L(a_i^{(k)}) \ge 0 \quad \text{for all } i,k,

and optionally to *strictly positive* solutions with `v_L(a_i^{(k)}) > 0`.

Algorithmic approach
--------------------

We attempt to bring the ideal `J` into **shape position** via restricted linear
changes of coordinates of the form

.. MATH::

    x_n \leftarrow x_n + b_1 x_1 + \cdots + b_{n-1} x_{n-1},\qquad x_i \leftarrow x_i \ (i<n),

with small integers `b_i`. We use lexicographic order with

.. MATH::

    x_1 > x_2 > \cdots > x_n,

so that in strict shape position the reduced Groebner basis has the form

.. MATH::

    x_1 - r_1(x_n),\ \ldots,\ x_{n-1} - r_{n-1}(x_n),\ f(x_n).

In that case, all solutions are parametrized by roots of `f(x_n)=0`, and the
remaining coordinates are obtained by evaluation

.. MATH::

    x_i = r_i(x_n).

We then use the existing MacLane-based univariate approximation machinery to
produce approximate roots of `f` over suitable extensions `(L,v_L)/(K,v_K)`,
and evaluate the `r_i` to obtain approximate solutions of the original system.

NOTE:
Since shape position is reached via a heuristic (partially probabilistic)
coordinate change, the function may raise an error after a specified number of
unsuccessful attempts.

EXAMPLES::

    sage: v2 = QQ.valuation(2)
    sage: R.<x,y> = QQ[]
    sage: J = R.ideal([x^2 + y^2 + 2, 2*x^2 + x*y + y])
    sage: S = approximate_solutions(J, v2)
    sage: len(S)
    2
    sage: s = S[0]
    sage: s.extension()
    ...
    sage: s.approximation()
    ...
    sage: s.improve_approximation()
    ...
"""

from sage.all import SageObject, ZZ


def approximate_solutions(J, v_K, positive_valuation=True, one_solution=False,
                          max_tries=8, bound=2, seed=1, try_list=None):
    r"""
    Return approximate solutions of a zero-dimensional system.

    INPUT:

    - ``J`` -- an ideal in a polynomial ring `K[x_1,...,x_n]` over a number field `K`;
               it is assumed that we have lex term order with `x_1 > \ldots > x_n`
    - ``v_K`` -- a `p`-adic valuation on `K`
    - ``positive_valuation`` -- boolean (default: ``True``); if set, require `v_L(a_i) > 0`
                                for all coordinates of the approximation
    - ``one_solution`` -- boolean (default: ``False``); if set, return only the first
                          solution found (or ``None`` if none is found)
    - ``max_tries`` -- positive integer; number of coordinate changes attempted
    - ``bound`` -- positive integer; random integers `b_i` are sampled from `[-bound, bound]`
    - ``seed`` -- integer seed for deterministic pseudo-random tries
    - ``try_list`` -- optional list of tuples `(b_1,...,b_{n-1})` tried before random tries

    OUTPUT:

    If ``one_solution`` is ``False``:
        a list of instances of :class:`ApproximateSolution`.

    If ``one_solution`` is ``True``:
        a single instance of :class:`ApproximateSolution`, or ``None``.

    The returned objects represent approximate solutions `a^(k) in L^n` (for some finite
    extension `(L,v_L)/(K,v_K)`) whose residual valuations can be made arbitrarily large
    by repeated calls to :meth:`ApproximateSolution.improve_approximation`.

    Raises an error in any of the following cases:
    
    - ``dim(J) != 0`` 
    - the term order id not `lex`
    - if shape position is not reached after ``max_tries`` attempts.
    """
    if try_list is None:
        try_list = []

    _check_dim_zero(J)
    _check_lex_order(J)

    # Step 1: attempt shape position
    J1, phi, f, rs = _put_in_shape_position(J, max_tries=max_tries, bound=bound,
                                            seed=seed, try_list=try_list)

    # Step 2: build one ApproximateSolution per "univariate branch" for f(x_n)=0
    sols = []
    for br in _branches_for_shape_polynomial(f, v_K):
        sol = ApproximateSolution(J1, phi, f, rs, v_K, br,
                                  positive_valuation=positive_valuation)
        if sol.is_admissible():
            if one_solution:
                return sol
            sols.append(sol)

    if one_solution:
        return None
    return sols


class ApproximateSolution(SageObject):
    r"""
    An object representing an *approximate solution* to a zero-dimensional system.

    Internally, a solution is represented by a univariate approximation branch for the
    shape polynomial `f(x_n)=0` over some finite extension `(L,v_L)/(K,v_K)`. The remaining
    coordinates are obtained by evaluating `x_i = r_i(x_n)`.

    Calling :meth:`improve_approximation` refines the approximation of `x_n`, and thereby
    improves the residual valuations of the defining equations.
    """

    def __init__(self, J1, phi, f, rs, v_K, branch, positive_valuation=True):
        self._J1 = J1
        self._phi = phi
        self._f = f
        self._rs = rs
        self._v_K = v_K
        self._branch = branch
        self._positive_valuation = positive_valuation

        self._R = J1.ring()
        self._xs = tuple(self._R.gens())
        self._x_shape = self._xs[-1]

        # Branch must provide (L,v_L) and a way to refine x_n
        self._extension = _branch_field(branch)
        self._extension_valuation = _branch_valuation(branch)

        self._approximation = None
        self._residual_min = None

        # initialize first approximation
        self.improve_approximation()

    def extension(self):
        r"""Return the field extension in which this approximate solution lives."""
        return self._extension

    def extension_valuation(self):
        r"""Return the valuation on the extension field."""
        return self._extension_valuation

    def approximation(self):
        r"""Return the current approximation vector `[a_1,...,a_n]`."""
        return self._approximation

    def residual_valuation(self):
        r"""Return `min_i v_L(F_i(a))` for the chosen generators at the current approximation."""
        return self._residual_min

    def is_admissible(self):
        r"""
        Return True if the current approximation satisfies the requested valuation constraint
        on the coordinates (integral or strictly positive).
        """
        vL = self._extension_valuation
        if vL is None or self._approximation is None:
            return False

        if self._positive_valuation:
            return all(vL(ai) > 0 for ai in self._approximation)
        else:
            return all(vL(ai) >= 0 for ai in self._approximation)

    def improve_approximation(self):
        r"""
        Improve the current approximation.

        OUTPUT:

        The improved approximation vector `[a_1,...,a_n]`.

        This method refines the underlying univariate approximation for `x_n`, evaluates
        `x_i=r_i(x_n)` for `i<n`, and updates the cached residual valuation.
        """
        x_shape = _branch_next_root_approx(self._branch)
        vec = []
        for r in self._rs:
            vec.append(r(x_shape))
        vec.append(x_shape)

        self._approximation = vec
        self._residual_min = _min_residual_valuation(self._J1, self._xs, vec,
                                                     self._extension_valuation)
        return self._approximation


# --------------------------------------------------------------------
# Helper functions (technical; keep at end of file)
# --------------------------------------------------------------------

def _check_dim_zero(J):
    d = J.dimension()
    if d != 0:
        raise ValueError("Expected dim(J)=0, got dim(J)=%s." % d)

def _check_lex_order(J):
    R = J.ring()
    if R.term_order().name() != 'lex':
        raise ValueError(
            "approximate_solutions requires lexicographic term order. "
            f"Got term order: {R.term_order()}."
        )

def _put_in_shape_position(J, max_tries=8, bound=2, seed=1, try_list=None):
    r"""
    Try restricted transforms x_n <- x_n + sum_{i<n} b_i x_i until strict shape form holds.

    OUTPUT:
        (J1, phi, f, rs) as in the module documentation.

    Raises NotImplementedError if unsuccessful after max_tries.
    """
    if try_list is None:
        try_list = []

    R = J.ring()
    xs = tuple(R.gens())
    n = len(xs)

    tries = [tuple(0 for _ in range(n - 1))]
    for t in try_list:
        tries.append(tuple(int(x) for x in t))

    gen = _random_shape_transforms(n, bound=bound, seed=seed)
    while len(tries) < max_tries:
        tries.append(next(gen))

    last_err = None
    for bs in tries:
        J1, phi = _apply_shape_transform_to_ideal(J, xs, bs)
        try:
            G = J1.groebner_basis()
            f, rs = _strict_shape_data_from_groebner(G, xs)
            return J1, phi, f, rs
        except Exception as e:
            last_err = e

    raise NotImplementedError(
        "Failed to put ideal into strict shape position after %s tries. Last error: %s"
        % (max_tries, last_err)
    )


def _strict_shape_data_from_groebner(G, xs):
    n = len(xs)
    x_shape = xs[-1]

    if len(G) != n:
        raise ValueError("Expected Groebner basis of length %s, got %s." % (n, len(G)))

    f = G[-1]
    if not _is_univariate_in(f, x_shape):
        raise ValueError("Not in shape position: last Groebner element not univariate in x_n.")

    r_map = {}
    for i in range(n-1):
        xi = xs[i]
        found = None
        for p in G:
            if p.degree(xi) != 1:
                continue
            if p.coefficient({xi: 1}) != 1:
                continue
            # must not involve xj (j<n) other than xi
            bad = False
            for xj in xs[:-1]:
                if xj != xi and (xj in p.variables()):
                    bad = True
                    break
            if bad:
                continue
            r = -p.subs({xi: 0})
            if not _is_univariate_in(r, x_shape):
                continue
            found = r
            break
        if found is None:
            raise ValueError("Missing equation x_%s - r(x_n)." % (i+1))
        r_map[i] = found

    rs = [r_map[i] for i in range(n-1)]
    return f, rs


def _is_univariate_in(p, var):
    vars_ = p.variables()
    return vars_ == () or vars_ == (var,)

def _apply_shape_transform_to_ideal(J, xs, bs):
    R = J.ring()
    x_shape = xs[-1]
    new_x_shape = x_shape + sum(ZZ(bs[i]) * xs[i] for i in range(len(xs)-1))
    phi = R.hom(list(xs[:-1]) + [new_x_shape])
    return phi(J), phi

def _random_shape_transforms(n, bound, seed):
    state = seed & 0xFFFFFFFF
    while True:
        bs = []
        for _ in range(n - 1):
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            b = int(state % (2 * bound + 1)) - bound
            bs.append(b)
        yield tuple(bs)


def _branches_for_shape_polynomial(f, v_K):
    r"""
    Return a list of branch objects describing approximate root clusters of f with v(root) >= 0
    (or >0, depending on your univariate machinery).

    This is the only place that must be connected to your MacLane-based code.
    """
    # TODO: replace by your actual API from approximate_factors.py
    from semistable_model.curves.approximate_factors import approximate_factors
    return list(approximate_factors(f, v_K))


def _branch_field(branch):
    return getattr(branch, "field", None)


def _branch_valuation(branch):
    return getattr(branch, "valuation", None)


def _branch_next_root_approx(branch):
    if hasattr(branch, "next_root_approx"):
        return branch.next_root_approx()
    if hasattr(branch, "__next__"):
        return next(branch)
    raise NotImplementedError("Branch object does not provide a refinement method for x_n.")


def _min_residual_valuation(J1, xs, vec, vL):
    subs = {xs[i]: vec[i] for i in range(len(vec))}
    # default: track residual on Groebner basis generators; you can change to J1.gens()
    vals = [vL(F.subs(subs)) for F in J1.groebner_basis()]
    return min(vals) if vals else None
