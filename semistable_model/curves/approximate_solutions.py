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

- It is assumed that J is radical. For 0-dimensional non-radical ideals, shape 
  position may fail and is currently not supported.
- Since shape position is reached via a heuristic (partially probabilistic)
  coordinate change, the function may raise an error after a specified number of
  unsuccessful attempts.


EXAMPLES::

    sage: v2 = QQ.valuation(2)
    sage: R = PolynomialRing(QQ, 2, ["x", "y"], order='lex')
    sage: x, y = R.gens()
    sage: J = R.ideal([x^2 + y^2 + 2, 2*x^2 + x*y + y])
    sage: S = approximate_solutions(J, v2)
    sage: len(S)
    1
    sage: s = S[0]; s
    approximate solution in 2 variables over ... (residual valuation ≥ 10)

    sage: s.extension()
    Number Field in alpha with defining polynomial y^2 + 40*y - 80

    sage: s.approximation()
    [-19/4*alpha - 15, 1389599/1867687*alpha - 354816/1867687]

    sage: s.residual_valuation()
    10
    sage: s.improve_approximation()
    [243270472557965/196*alpha - 810932028059849/343, -839*alpha + 64/7]

    sage: s.residual_valuation()
    13

"""

from sage.all import SageObject, ZZ, QQ, PolynomialRing
from semistable_model.curves.approximate_factors import approximate_roots


def approximate_solutions(J, v_K, positive_valuation=True, one_solution=False,
                          max_tries=8, bound=2, seed=1, try_list=None, 
                          check_is_radical=True):
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
    - ``check_is_radical`` -- boolean (default: ``False``); if set, test whether the ideal
                              `J` is radical  

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
    - ``check_is_radical``is ``True`` and `J` is not radical
    - the term order id not `lex`
    - if shape position is not reached after ``max_tries`` attempts.
    """
    if try_list is None:
        try_list = []

    _check_dim_zero(J)
    # this seems to take too long in some examples
    if check_is_radical:
        _check_is_radical(J)
    _check_lex_order(J)

    # Step 1: attempt shape position
    phi, f, rs, bs = _put_in_shape_position(J, max_tries=max_tries, bound=bound,
                                            seed=seed, try_list=try_list)

    # Step 2: build one ApproximateSolution per "univariate branch" for f(x_n)=0
    sols = []
    for br in approximate_roots(f, v_K, 
                                positive_valuation=positive_valuation):
        sol = ApproximateSolution(J, phi, f, rs, v_K, br, bs,
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
    def __init__(self, J, phi, f, rs, v_K, branch, bs, positive_valuation=True):
        self._J = J
        self._phi = phi
        self._f = f
        self._rs = rs
        self._v_K = v_K
        self._branch = branch
        self._positive_valuation = positive_valuation

        self._R = J.ring()
        self._xs = tuple(self._R.gens())
        self._x_shape = self._xs[-1]

        # store shear coefficients b_1,...,b_{n-1}
        self._bs = tuple(int(b) for b in bs)
        if len(self._bs) != len(self._xs) - 1:
            raise ValueError("bs must have length n-1.")

        self._extension = branch.extension()
        self._extension_valuation = branch.extension_valuation()

        self._approximation = None
        self._residual_min = None

        self.improve_approximation()

    def __repr__(self):
        """
        String representation of an approximate solution.
        """
        n = len(self._xs)

        try:
            L = self._extension
            field_str = f"over {L}"
        except Exception:
            field_str = "over unknown extension"

        try:
            prec = self._residual_min
            if prec is None:
                prec_str = "residual valuation not computed"
            else:
                prec_str = f"residual valuation ≥ {prec}"
        except Exception:
            prec_str = "residual valuation unknown"

        return (
            f"approximate solution in {n} variables "
            f"{field_str} "
            f"({prec_str})"
        )

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

    def improve_approximation(self):
        r"""
        Improve the current approximation.

        OUTPUT:

        The improved approximation vector `[a_1,...,a_n]` in the original coordinates.

        This method refines the univariate approximation for the shape variable,
        constructs the corresponding solution in shape coordinates, maps it back
        via the inverse coordinate change, and updates the cached residual valuation.
        """
        br = self._branch  # the shape root 
        # Refine the shape root (does nothing if exact)
        br.improve_approximation()

        # Shape-coordinate for x_n' (i.e. y_n)
        y_shape = br.approximation()

        # Solution in shape coordinates (for J1)
        # y_i = r_i(y_n) for i < n, y_n = y_shape
        y = [br.eval(r, simplify=True) for r in self._rs] + [y_shape]

        # Map back to original coordinates via phi
        phi = self._phi
        subs = {self._xs[i]: y[i] for i in range(len(self._xs))}
        x = [phi(self._xs[i]).subs(subs) for i in range(len(self._xs))]

        self._approximation = x
        self._residual_min = _min_residual_valuation(
            self._J, self._xs, x, self._extension_valuation
        )

        return self._approximation

    def is_admissible(self):
        r"""
        Limit-based admissibility check using ApproximateRoot.value() and value_of_poly().
        """
        br = self._branch

        # valuations of x_1,...,x_{n-1}: x_i = r_i(y_n)
        vals = [br.value_of_poly(r) for r in self._rs]

        # valuation of x_n = y_n - sum b_i r_i(y_n)
        R = self._rs[0].parent()   # univariate ring K[x_n]
        t = R.gen()
        q = t
        for i, r in enumerate(self._rs):
            q -= ZZ(self._bs[i]) * R(r)
        vals.append(br.value_of_poly(q))

        if self._positive_valuation:
            return all(v > 0 for v in vals)
        else:
            return all(v >= 0 for v in vals)


# --------------------------------------------------------------------
# Helper functions (technical; keep at end of file)
# --------------------------------------------------------------------


def _check_dim_zero(J):
    d = J.dimension()
    if d != 0:
        raise ValueError("Expected dim(J)=0, got dim(J)=%s." % d)
    

def _check_is_radical(J):
    J_r = J.radical()
    if not all(J.reduce(g) == 0 for g in J_r.gens()):
        raise ValueError("J is not radical.")
 

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
        (phi, f, rs, bs) as in the module documentation.

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
            return phi, f, rs, bs
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
    # now we can coerce f into a univariate polynomial ring
    f = f.univariate_polynomial()

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
            found = r.univariate_polynomial()
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


def _min_residual_valuation(J, xs, vec, vL):
    subs = {xs[i]: vec[i] for i in range(len(vec))}
    vals = [vL(F.subs(subs)) for F in J.gens()]
    return min(vals) if vals else None


# -----------------------------------------------------------------------------

#          Test functions

def test_approximate_solutions(v_K=QQ.valuation(2), N=10, d=3, n=3, prec=10, n_tries=10):
    K = v_K.domain()
    R = PolynomialRing(K, n, "x", order="lex")
    for _ in range(N):
        while True:
            J = R.ideal([R.random_element(d) for _ in range(n)])
            try:
                _check_dim_zero(J)
                _check_is_radical(J)
                break
            except:
                pass
        print(J)
        S = approximate_solutions(J, v_K, positive_valuation=False)
        for s in S:
            print(f"s = {s}")
            if s.residual_valuation() < prec:
                for _ in range(n_tries):
                    t = s.residual_valuation()
                    print(f"t = {t}")
                    s.improve_approximation()
                    if t > prec:
                        break
                else:
                    raise ValueError("residual doesn't increase enough")
        print()
