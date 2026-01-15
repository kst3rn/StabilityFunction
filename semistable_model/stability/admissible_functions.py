r"""
Admissible functions
=====================

Let `K` be a field equipped with a discrete valuation `v_K` and let `S:=K[t]`
be the polynomial ring over `K`. A *MacLane pseudovaluation* on `S` is a
pseudovaluation `v:S\to \mathbb{R}\cup\{\infty\}` whose restriction to `K` is
`v_K`, satisfies `v(t)\geq 0` and has a representation as an inductive valuation,
or as a limit pseudovaluation.

It is convenient to identify a MacLane pseudovaluation with the corresponding
point on the Berkovich unit disk `X` over the completion of `K`. We shall
typically write `\xi` for a point on `X` and `v=v_\xi` for its corresponding
pseudovaluation on `S`.

The Berkovich unit disk is a partially ordered set with respect to the order

.. MATH::

    \xi\leq \xi' \quad \Leftrightarrow v_\xi(f) \leq v_{\xi'}(f) \quad \forall f\in S.

Its unique minimal element is the Gauss point `\xi_G`, corespoding to the Gauss
valuation `v_G`. Its maximal elements are the points of type I, which correspond
to the pseudovaluation which are not valuations.

For an element `\xi\in X` we set

.. MATH::

    D_\xi := \{ \xi'\in X : \xi'\geq \xi \}.

This is called a *closed diskoid*

A *valuative function* on `X` is a function `h:X\to \mathbb{R}` of the form

.. MATH::

    h(\xi) = a\cdot v_\xi(f),

with `a\in\mathbb{Q}`, `a>0`, and `f\in S`. An *admissible function* is a function
`h:X\to \mathbb{R}` which is the minimum of a finite set of valuative functions.

We note that an admissible function is increasing, i.e. `\xi\leq \xi'` implies `h(\xi)\leq h(\xi')`.

The goal of this module is to implement a class `AdmissibleFunction` which represents
admissible functions on `X`. The most important subgoal is to implement a function which
computes the subset of `X` where a given admissible function achieves its maximum.

"""

from sage.all import SageObject, QQ, Infinity, FunctionField, Polynomial, Factorization
from mclf import BerkovichLine, BerkovichTree


# -----------------------------------------------------------------------------
# Admissible functions


class AdmissibleFunction(SageObject):
    r"""
    An admissible function on the closed unit disk.

    Let `K` be a field equipped with a discrete valuation `v_K` and let `S:=K[t]`
    be the polynomial ring over `K`. A MacLane pseudovaluation on `S` is a
    pseudovaluation `v:S\to \mathbb{R}\cup\{\infty\}` whose restriction to `K` is
    `v_K`, satisfies `v(t)\geq 0` and has a representation as an inductive valuation,
    or as a limit pseudovaluation.

    Let `X` denote the set of all MacLane pseudovaluations on `S`. This is a subset
    of the closed unit disk, considered as a Berkovich annalytic space over the
    completion of `K`. The elements of `X` represent the points of the Berkovich
    unit disk which are in some sense computable.

    In [OssenWewers]_ a general theory of admissible functions on `K`-analytic curves
    is developed. This class implements a special case of this notion.
    We allow two kinds of admissible functions:

    - A single valuative function a_1*v(f_1) + ... + a_n*v(f_n), where the `f_i` are
      *polynomials* and the `a_i` are *positive* rational numbers.
    - A minimum of several such valuative functions.

    These are realized, respectively, by the following subclasses:

    - :class:`ValuativeFunction`
    - :class:`MinimumOfValuativeFunctions`

    Every admissible function is defined on a *closed diskoid* `D_\xi` in the Berkovich
    unit disk, where `\xi` is a point of type II, corresponding to a MacLane valuation.
    The point `\xi` is called the *root* of the admissible function.

    Let `h` be an admissible function with root `\xi_0`. A *skeleton* for `h` is a
    finite tree `T` in the diskoid `D_{\xi_0}` such that

    - the restriction of `h` to any edge of `T` is an affine function,
    - `h` is locally constant on the complement of `T`.

    Such a skeleton always exist. For instance, the tree with root `\xi_0`, spanned
    by the roots of the *prime factors* of `h` lying in `D_{v_0}` is a skeleton.
    Here a *prime factor* of `h` is any monic irreducible polynomial `f` dividing
    any of the polynomials `f_i` in the valuative functions defining `h`.

    """

    def base_field(self):
        r"""
        Return the base field of this admissible function.

        OUTPUT:

        The base field of this admissible function.
        """
        return self.polynomial_ring().base_ring()

    def polynomial_ring(self):
        r"""
        Return the polynomial ring of this admissible function.

        OUTPUT:

        The polynomial ring of this admissible function.
        """
        return self.root().domain()

    def function_field(self):
        r"""
        Return the function field of this admissible function.

        OUTPUT:

        The function field of this admissible function.
        """
        return self.berkovich_line().function_field

    def root(self):
        r"""
        Return the root of the domain of this admissible function.

        OUTPUT:

        A point of type II `\xi\in X` such that the closed diskoid `D_\xi`
        is the domain of this admissible function.

        """
        return self._root

    def prime_factors(self):
        r"""
        Return the prime factors of the admissible function.

        OUTPUT:

        A list of all prime factors which are active in this admissible function.
        Being *active* means that the primne factor has a root in the domain of
        this admissible function.
        """
        return self._prime_factors

    def berkovich_line(self):
        r"""
        Return the Berkovich line of the domain of this admissible function.

        OUTPUT:

        The Berkovich line of the domain of this admissible function.

        """
        return self._berkovich_line

    def skeleton(self):
        r"""
        Return the skeleton of this admissible function.

        """
        if not hasattr(self, '_skeleton'):
            self._make_skeleton()
        return self._skeleton

    def _make_skeleton(self):
        X = self.berkovich_line()
        F = X.function_field()
        xi0 = self.root()
        T = BerkovichTree(X, root=xi0)
        for f in self.prime_factors():
            for xi, _ in X.prime_divisor(F(f)):
                if xi0.is_leq(xi):
                    T, _ = T.add_point(xi)
        self._skeleton = T

    def restriction_to_subtree(self, T):
        r""" Return the restriction of this admissible function to the subtree `T`.

        INPUT:

        - ``T`` -- a subtree of the skeleton of this admissible function.

        OUTPUT:

        The restriction of this admissible function to the discoid `D_\xi` where
        `\xi` is the root of the subtree `T`.

        This method is abstract and must be implemented by the subclasses.

        """
        raise NotImplementedError()

    def maximum(self):
        r"""
        Return the maximum of this admissible function.

        OUTPUT:

        A pair `(m, M)` where `m` is the maximum of this admissible function
        and `M` is the list of  points `\xi` on the domain such that the
        maximum `m` is attained precisely on the closed diskoids `D_\xi`.

        Note that `\xi` may also be a point of type I, in which case
        `D_\xi=\{\xi\}`.

        This method is abstract and must be implemented by the subclasses.
        """
        raise NotImplementedError()


class ValuativeFunction(AdmissibleFunction):
    r""" A valuative function of the form a_0 + a_1*v(f_1) + ... + a_n*v(f_n).

    INPUT:

    - ``a0`` -- a rational number,
    - ``coefficients`` -- a tuple (a_1, ..., a_n) of positive rational numbers,
    - ``functions`` -- a list of monic, irreducible polynomials f_1, ..., f_n,
    - ``root`` -- a point of type II on the Berkovich unit disk, or ``None``,
    - ``skeleton`` -- a tree in the Berkovich unit disk, or ``None``.

    It is assumed that the irreducible polynomials `f_i` are pairwise distinct
    and have positive effective degree with respect to the root.

    If ``root`` is not given, it is assumed that it is equal to the root of the
    tree ``skeleton``. If ``skeleton`` is not given, then a skeleton is
    constructed from the prime factors of the function. Otherwise it is assumed
    that ``skeleton`` is actually a skeleton for the function.

    If neither ``root`` nor ``skeleton`` is given, an error is raised.

    We note the following:

    - a valuative function of this form is increasing,
    - every valuative function defined on a closed diskoid which is increasing
      can be brought into this form.
    - we may further assume that valuative functions `v(f_i)` are *active*,
      i.e. nonconstant on the domain of the function, the closed diskoid given
      by ``root``.

    EXAMPLES::

        sage: from sage.all import QQ
        sage: from mclf import BerkovichLine
        sage: import admissible_functions
        sage: from admissible_functions import ValuativeFunction
        sage: F.<t> = FunctionField(QQ)
        sage: X = BerkovichLine(F, QQ.valuation(2))
        sage: f = t^2 + 1
        sage: xi0 = X.gauss_point()
        sage: h = ValuativeFunction(0, [1], [f], xi0)
        sage: h
        1*v(t^2 + 1)

        sage: h.skeleton()
        Berkovich tree with 2 vertices

    """

    def __init__(self, constant, coefficients, functions, root=None, skeleton=None):
        r"""
        Initialize the valuative function.

        We check that the input already is in normal form, so that each
        entry of ``functions`` is a monic, irreducible polynomial with positive
        effective degree with respect to the root. This implies that the function
        is constant iff the list of functions is empty.
        """
        if root is None:
            assert skeleton is not None, "Either root or skeleton must be given."
            root = skeleton.root()
        X = root.berkovich_line()
        F = X.function_field()
        functions = [F(f) for f in functions]

        assert constant in QQ
        assert len(coefficients) == len(functions)
        assert all(f.denominator().is_one() and f.numerator().is_irreducible()
                   and f.numerator().is_monic() and
                   effective_degree(f, root) > 0 for f in functions
                   )  # the latter conditon checks that v(f) is nonconstant
        assert all(a in QQ and a > 0 for a in coefficients)

        self._root = root
        self._constant = constant
        self._coefficients = coefficients
        self._functions = functions
        self._prime_factors = set(functions)
        self._berkovich_line = X
        if skeleton is None:
            self._make_skeleton()
        else:
            self._skeleton = skeleton

    def __repr__(self):
        if self.is_constant():
            return str(self.constant())
        if self._constant == 0:
            return "+ ".join([f"{coeff}*v({func})" for coeff, func
                              in zip(self._coefficients, self._functions)])
        return " + ".join(
            [str(self._constant)] +
            [f"{coeff}*v({func})" for coeff, func in zip(self._coefficients, self._functions)]
        )

    def __call__(self, xi):
        r"""
        Return the value of this valuative function at the pseudovaluation `xi`.

        INPUT:

        - ``xi`` -- a point on the domain of this valuative function.

        OUTPUT:

        The value of this valuative function at ´\xi´.

        EXAMPLES::

            sage: from sage.all import QQ
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import ValuativeFunction
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: f = t^2 + 1
            sage: xi0 = X.gauss_point()
            sage: h = ValuativeFunction(0, [1], [f], xi0)
            sage: h(X.gauss_point())
            0
        """
        if self.is_constant():
            return self.constant()
        return self._constant + sum(a*xi.v(f) for a, f in zip(self._coefficients, self._functions))

    def is_constant(self):
        r"""
        Return True if this valuative function is constant.

        OUTPUT:

        True if this function is constant, False otherwise.
        """
        return len(self._functions) == 0

    def constant(self, check=True):
        r"""
        Return the constant value of this constant valuative function.

        OUTPUT:

        The constant value of this constant valuative function.

        If this function is not constant, an error is raised, unless check=False.
        """
        if check:
            assert self.is_constant(), "This valuative function is not constant."
        return self._constant

    def value_at_root(self):
        r"""
        Return the value of this valuative function at the root.

        OUTPUT:

        The value of this valuative function at the root.
        """
        return self(self.root())

    def restriction_to_subtree(self, T):
        r"""
        Return the restriction of this valuative function to the subtree `T`.

        INPUT:

        - ``T`` -- a subtree of the skeleton of this valuative function.

        OUTPUT:

        The restriction of this valuative function to the discoid `D_\xi` where
        `\xi` is the root of the subtree `T`.

        EXAMPLES::

            sage: from sage.all import QQ
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import valuative_function
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: f = t*(t^2 + 2)
            sage: xi0 = X.gauss_point()
            sage: h = valuative_function(f, xi0)
            sage: T = h.skeleton()
            sage: T1 = T.children()[0]
            sage: h.restriction_to_subtree(T1)
            1*v(t)+ 1*v(t^2 + 2)

        """
        xi = T.root()
        assert xi.type() == "II", "the root of T must be of type II"
        assert self.root().is_leq(xi), "T is not a subtree of the skeleton of this function."
        constant = self._constant
        coefficients = []
        functions = []
        for a, f in zip(self._coefficients, self._functions):
            if effective_degree(f, xi) > 0:
                coefficients.append(a)
                functions.append(f)
            else:
                constant += a*xi.v(f)
        return ValuativeFunction(constant, coefficients, functions, root=None, skeleton=T)

    def maximum(self):
        r"""
        Return the maximum of this admissible function.

        OUTPUT:

        A pair `(m, M)` where `m` is the maximum of this admissible function
        and `M` is the list of points `\xi` such that the
        maximum `m`is attained precisely on the closed diskoids `D_\xi`.

        EXAMPLES::

            sage: from sage.all import QQ
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import valuative_function
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: f = 2*(t^2 + 1)*t
            sage: xi0 = X.gauss_point()
            sage: h = valuative_function(f, xi0)
            sage: h.maximum()
            (+Infinity,
             [Point of type I on Berkovich line given by t^2 + 1 = 0,
              Point of type I on Berkovich line given by t = 0])

        """
        if self.is_constant():
            return self.constant(), [self.root()]
        # the maximum is infty, and it is attained in the roots of the prime factors
        M = []
        for f in self.prime_factors():
            for xi, _ in self.berkovich_line().prime_divisor(f):
                if self.root().is_leq(xi):
                    M.append(xi)
        return Infinity, M

    def affine_function(self, p):
        r"""
        Return the affine function a + b*t which is the restriction of this
        valuative function to the path `p`.

        INPUT:

        - ``p`` -- a path in the Berkovich line of the domain of this valuative
                    function.

        OUTPUT:

        The pair `(a, b)` representing the affine function a + b*t which is the
        restriction of this valuative function to the path `p`.

        It is assumed (and not checked) that the path `p` is contained in an edge
        of the Berkovich tree spanned by the prime factors of the valuative
        function. This implies that the restriction of the function to the path
        is an affine function.

        Suppose that the path `p` has initial point `v_0` and end point `v_1`.
        Assume first that `v_1` in inductive, determined by `v_1(\phi)=s_1`,
        where `\phi` is a key polynomial for `v_1` and `s_1` is rational or
        `\infty`. Then the parametrization of the path is given by

        .. MATH::

            s\in[s_0,s_1] \mapsto v_s, \quad v_s(\phi)=s.

        Let `f` be an irreducible polynomial. Then `v_s(f)` is a continuous,
        piecewise affine function, with right derivative at `s` equal to the
        effective degree of `f` with respect to `v_s`. By assumption, the prime
        factors of the valuative function define affine functions on the path,
        so their effective degree is constant. This allows us to easily compute
        the pair `(a, b)`.

        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import valuative_function, Path
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: f = 2*t*(t^2 + 1)
            sage: xi0 = X.gauss_point()
            sage: h = valuative_function(f, xi0)
            sage: xi1 = X.point_from_discoid(t, 1)
            sage: p = Path(xi0, xi1)
            sage: h.affine_function(p)
            (1, 1)

        """
        s0, s1 = p.boundaries()
        b = QQ(0)
        for coeff, f in zip(self._coefficients, self._functions):
            b += coeff*p.derivative(f)
        a = self.value_at_root() - b*s0
        # test the result (this test seems to slow down computations a lot)
        # s = s0
        # while s <= s1 and s < 10:
        #     xi = p.point(s)
        #     assert self(xi) == a + b*s, f"Error at s={s}: f = {self}, p = {p}, xi = {xi}"
        #     s += 1
        return a, b


class MinimumOfValuativeFunctions(AdmissibleFunction):
    r""" A minimum of several valuative functions.

    INPUT:

    - ``functions`` -- a list of :class:`ValuativeFunction` instances,
    - ``constant`` -- a rational number, or ``Infinity`` (default: ``Infinity``).
    - ``root`` -- a MacLane valuation on K[t]. (default: ``None``)
    - ``skeleton`` -- a tree in the Berkovich unit disk, or ``None``.

    OUTPUT:

    A :class:`MinimumOfValuativeFunctions` instance, representing the minimum of
    the given valuative functions and the constant function with value ``constant``.

    If any of the functions in ``functions`` is constant, the constant value is
    updated and the function is removed from the list.

    If ``root`` is not given, it is assumed that it is equal to the root of the
    tree ``skeleton``. If ``skeleton`` is not given, then a skeleton is constructed
    from the prime factors of the function. Otherwise it is assumed that ``skeleton``
    is actually a skeleton for the function.

    If neither ``root`` nor ``skeleton`` is given, the root is taken to be
    the root of the first function in the list.

    EXAMPLES::

        sage: from sage.all import QQ, FunctionField
        sage: from mclf import BerkovichLine
        sage: from admissible_functions import valuative_function, MinimumOfValuativeFunctions
        sage: F.<t> = FunctionField(QQ)
        sage: X = BerkovichLine(F, QQ.valuation(2))
        sage: f = t^2 + 1
        sage: xi0 = X.gauss_point()
        sage: h1 = valuative_function(t^2 + 1, xi0)
        sage: h2 = valuative_function(2, xi0)
        sage: h = MinimumOfValuativeFunctions([h1, h2])
        sage: h
        min(1, 1*v(t^2 + 1))

        sage: h3 = valuative_function(2*(t^2 + 1)*t, xi0)
        sage: h4 = valuative_function(1, xi0)
        sage: h = MinimumOfValuativeFunctions([h3, h4])
        sage: h
        0

    """

    def __init__(self, functions, constant=Infinity, root=None, skeleton=None):
        r"""
        Initialize the minimum of a list of valuative functions.

        We have to make sure that the function is either recognizably constant
        (i.e. the list _functions is empty) or that there is one direction in
        which the function grows.
        """
        assert all(isinstance(f, ValuativeFunction) for f in functions)
        if root is None:
            if skeleton is None:
                root = functions[0].root()
            else:
                root = skeleton.root()
        self._root = root
        # we update the constant by taking the minimum of the constants functions
        for f in functions:
            if f.is_constant():
                constant = min(constant, f.constant())
        prime_factors = set()
        nonconstant_functions = []
        for f in functions:
            if f.value_at_root() < constant:
                # since valuative functions are increasing, we can ignore
                # functions which are larger than the constant
                nonconstant_functions.append(f)
                prime_factors.update(f.prime_factors())
        self._prime_factors = prime_factors
        self._constant = constant
        self._functions = nonconstant_functions
        self._berkovich_line = root.berkovich_line()
        if skeleton is None:
            self._make_skeleton()
        else:
            self._skeleton = skeleton

    def __repr__(self):
        if self._constant == Infinity:
            if len(self._functions) == 0:
                return "Infinity"
            if len(self._functions) == 1:
                return self._functions[0]
            return "min(" + ", ".join(str(f) for f in self._functions) + ")"
        if len(self._functions) == 0:
            return str(self._constant)
        return "min(" + str(self._constant) + ", " + \
            ", ".join(str(f) for f in self._functions) + ")"

    def __call__(self, xi):
        r"""
        Return the value of this admissible function at the point `\xi`.

        INPUT:

        - ``xi`` -- a point on the domain of this admissible function.

        OUTPUT:

        The value of this admissible function at `\xi`.

        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import valuative_function, MinimumOfValuativeFunctions
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: f = 2*(t^2 + 1)*t
            sage: xi0 = X.gauss_point()
            sage: h1 = valuative_function(t^2 + 1, xi0)
            sage: h2 = valuative_function(2, xi0)
            sage: h = MinimumOfValuativeFunctions([h1, h2])
            sage: xi = X.point_from_discoid(t, 1)
            sage: h(xi)
            0
            sage: xi = X.point_from_discoid(t-1, 1)
            sage: h(xi)
            1

        """
        if self.is_constant():
            return self.constant()
        return min([self.constant()] + [f(xi) for f in self._functions])

    def is_constant(self):
        r"""
        Return True if this valuative function is constant.

        OUTPUT:

        True if this function is constant, False otherwise.
        """
        return len(self._functions) == 0

    def constant(self):
        r""" Return the constant part of this admissible function.

        This admissible function is of the form `\min(a, h_1, ..., h_n)`, where
        `a` is a constant and `h_1,\ldots,h_n` are nonconstant valuative functions.

        This method returns `a`.

        """
        return self._constant

    def nonconstant_functions(self):
        r""" Return the list of nonconstant valuative functions in this admissible function.

        This admissible function is of the form `\min(a, h_1, ..., h_n)`, where
        `a` is a constant and `h_1,\ldots,h_n` are nonconstant valuative functions.

        This method returns `[h_1, ..., h_n]`.

        """
        return self._functions

    def value_at_root(self):
        r"""
        Return the value of this valuative function at the root.

        """
        return self(self.root())

    def restriction_to_subtree(self, T):
        r"""
        Return the restriction of this admissible function to the subtree `T`.

        INPUT:

        - ``T`` -- a subtree of the skeleton of this admissible function.

        OUTPUT:

        The restriction of this admissible function to the discoid `D_\xi` where
        `\xi` is the root of the subtree `T`.

        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import valuative_function, MinimumOfValuativeFunctions
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: xi0 = X.gauss_point()
            sage: h1 = valuative_function(t*(t^2 + 2)*(t+1), xi0)
            sage: h2 = valuative_function(2, xi0)
            sage: h = MinimumOfValuativeFunctions([h1, h2], root=xi0)
            sage: h
            min(1, 1*v(t)+ 1*v(t + 1)+ 1*v(t^2 + 2))

            sage: h.skeleton()
            Berkovich tree with 5 vertices

            sage: T1 = h.skeleton().children()[0]
            sage: h.restriction_to_subtree(T1)
            1

            sage: T2 = h.skeleton().children()[1]
            sage: h.restriction_to_subtree(T2)
            Traceback (most recent call last):
            ...
            AssertionError: the root of T must be of type II

        """
        c = self._constant
        nonconstant_functions = []
        for f in self._functions:
            f_T = f.restriction_to_subtree(T)
            if f_T.is_constant():
                c = min(c, f_T.constant())
            else:
                nonconstant_functions.append(f_T)
        return MinimumOfValuativeFunctions(nonconstant_functions, c, skeleton=T)

    def maximum(self):
        r"""
        Return the maximum of this admissible function.

        OUTPUT:

        A pair `(m, M)` where `m` is the maximum of this admissible function
        and `M` is the list of points `\xi` such that the
        maximum `m`is attained precisely on the closed diskoids `D_\xi`.
        """
        if self.is_constant():
            return self.value_at_root(), [self.root()]
        local_maxima = self.local_maxima()
        m = Infinity
        M = []
        for m1, xi1 in local_maxima:
            if m1 > m:
                m = m1
                M = [xi1]
            elif m1 == m:
                M.append(xi1)
        return m, M

    def local_maxima(self):
        r"""
        Return the local maxima of this admissible function.

        OUTPUT:

        A list of pairs `(m, xi)` where `m` is a local maximum of this
        admissible function and `xi` is the point such
        that the local maximum `m` is attained precisely on the closed diskoid
        `D_\xi`.

        ALGORITHM:

        If this admissible function is constant, the unique local maximum is
        attained on the root.

        Let `T` denote the skeleton of the
        admissible function `h`. Then the function `h` is constant outside
        of `T` and piecewise affine and concave on the edges of `T`.

        The local maxima of the admissible function are
        precisely the leafs of the subtree of `T` where the function is
        nonconstant.

        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import valuative_function, MinimumOfValuativeFunctions
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: xi0 = X.gauss_point()
            sage: h1 = valuative_function(2*(t^2 +1)*t, xi0)
            sage: h2 = valuative_function(4, xi0)
            sage: h = MinimumOfValuativeFunctions([h1, h2], root=xi0)
            sage: h.local_maxima()
            [(2, Point of type II on Berkovich line, corresponding to v(t + 1) >= 1/2),
             (2, Point of type II on Berkovich line, corresponding to v(t) >= 1)]

        """
        if self.is_constant():
            return [(self.constant(), self.root())]
        # if this admissible function is not constant, the maximum can't be
        # attained on the root. But the local maxima all lie on the
        # berkovich tree spanned by the prime factors.
        T = self.skeleton()
        ret = []
        for T1 in T.children():
            xi1 = T1.root()
            if xi1.type() == "I":
                ret.append(self.local_maximum_on_branch(xi1))
            else:
                h1 = self.restriction_to_subtree(T1)
                if not h1.is_constant():
                    ret += h1.local_maxima()
                else:
                    # there is a unique local maximum on the path
                    # from the root to the child v1
                    ret.append(self.local_maximum_on_branch(xi1))
        return ret

    def local_maximum_on_branch(self, xi1):
        r"""
        Return the local maximum of this admissible function on the branch
        spanned by the root and the point `xi1`.

        INPUT:

        - ``xi`` -- a point on the domain of this admissible function,
                    strictly larger than the root.

        OUTPUT:

        A pair `(m, xi)` where `m` is the local maximum of this admissible
        function on the branch spanned by `\xi_1`, i.e. the path from the
        root to `\xi_1`.

        It is assume that `\xi_1` is a child node of the root `\xi_0` in the
        skeleton of this admissible function. This implies that the
        restriction of the function to the branch is the minimum of affine
        functions.

        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import valuative_function, MinimumOfValuativeFunctions
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: xi0 = X.gauss_point()
            sage: h1 = valuative_function(t^2 - 4, xi0)
            sage: h2 = valuative_function(2, xi0)
            sage: h = MinimumOfValuativeFunctions([h1, h2], root=xi0)
            sage: xi1 = h.skeleton().children()[0].root()
            sage: h.local_maximum_on_branch(xi1)
            (1, Point of type II on Berkovich line, corresponding to v(t) >= 1/2)

        """
        if self.is_constant():
            return self.constant(), self.root()
        if self(xi1) == Infinity:
            return Infinity, xi1

        path = Path(self.root(), xi1)
        s0, s1 = path.boundaries()
        affine_functions = [(self.constant(), QQ.zero())]
        for f in self.nonconstant_functions():
            affine_functions.append(f.affine_function(path))
        m, s = max_of_min(affine_functions, s0, s1)
        xi = path.point(s)
        assert self(xi) == m, (f"something is wrong: path = {path}, s0 = {s0}, "
                               + "s1 = {s1}, m = {m}, s = {s}, xi = {xi}, "
                               + "affine_functions = {affine_functions}")
        return m, xi


class Path(SageObject):
    r"""
    A path in the Berkovich line.

    INPUT:

    - ``xi0`` -- a point of type II on the Berkovich unit disk.
    - ``xi1`` -- a point on the Berkovich unit disk which is strictly larger than `\xi_0`.

    OUTPUT:

    A :class:`Path` instance representing the path from `\xi_0` to `\xi_1`.

    The main purpose of this class is to be able to have an explicit parametrization
    of this path, and to be able to compute the derivative of a valuative function
    with respect to this parametrization.

    Let `v_i` be the MacLane pseudovaluation corresponding to `\xi_i`, `i=0,1`.
    Suppose first that `v_1` is inductive, determined by `v_1(\phi)=s_1`, where
    `\phi` is a key polynomial for `v_1` and `s_1` is rational or `\infty`. Then
    we use the parametrization

    .. MATH::

        s\in[s_0,s_1] \mapsto v_s, \quad v_s(\phi)=s,

    where `s_0 = v_0(\phi)`.

    If `v_1` is a limit valuation, we need to find a sufficiently good approximation
    of `v_1` whose key polynomial `\phi` has maximal degree. We can then use the
    parametrization given by `\phi` up to the chosen approximation of `v_1`. In
    order to go beyond this approximation, we first have to improve this and adapt
    the choice of `\phi`.

    EXAMPLES::

        sage: from sage.all import QQ, FunctionField
        sage: from mclf import BerkovichLine
        sage: from admissible_functions import Path
        sage: F.<t> = FunctionField(QQ)
        sage: X = BerkovichLine(F, QQ.valuation(2))
        sage: xi0 = X.gauss_point()
        sage: xi1 = X.point_from_discoid(t, 1)
        sage: p = Path(xi0, xi1)
        sage: p
        Path from Point of type II on Berkovich line, corresponding to v(t) >= 0 to Point of type II on Berkovich line, corresponding to v(t) >= 1

    Here is an example that caused trouble previously:

    """

    def __init__(self, xi0, xi1):
        r"""
        Initialize the path.
        """
        assert xi0.type() == "II" and xi0.is_leq(xi1), "The input is not valid."
        # compute a key polynomial to use for the parametrization
        if xi1.is_inductive():
            # xi1 is inductive, so there is a natural choice for the key polynomial
            phi, s1 = xi1.discoid()
            s0 = xi0.v(phi)
            self._xi1_is_inductive = True
        else:
            # xi1 is a limit point; we need a sufficiently good approximation
            # whose key polynomial phi has the maximal degree
            # this can be checked by computing the effective degree of phi
            phi, _ = xi1.improved_approximation().discoid()
            while effective_degree(phi, xi1.approximation()) > 1:
                phi, _ = xi1.improved_approximation().discoid()
            s0 = xi0.v(phi)
            s1 = Infinity
            self._s1_approx = xi1.v(phi)
            self._v1_is_inductive = False
        self._xi0 = xi0
        self._xi1 = xi1
        self._s0 = s0
        self._s1 = s1
        self._phi = phi

    def __repr__(self):
        return f"Path from {self._xi0} to {self._xi1}"

    def initial_point(self):
        r"""
        Return the initial point of this path.

        """
        return self._xi0

    def boundaries(self):
        r"""
        Return the boundaries of the interval parametrizing this path.

        OUTPUT:

        A pair `(s_0, s_1)` of rational numbers such that the path is
        parametrized by the interval `[s_0, s_1]`. The value of `s_1`
        may be `Infinity`.

        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import Path
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: xi0 = X.gauss_point()
            sage: xi1 = X.point_from_discoid(t, 1)
            sage: p = Path(xi0, xi1)
            sage: p.boundaries()
            (0, 1)
        """
        return self._s0, self._s1

    def derivative(self, f, s=None):
        r"""
        Return the derivative of the function `f` along this path.

        INPUT:

        - ``f`` -- a polynomial in `K[t]`.
        - ``s`` -- a rational number in the interval `[s_0, s_1]`.

        OUTPUT:

        The right derivative of the function `f` at the point `s`.
        If `s` is not given, the derivative is computed at the initial point.


        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import Path
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: xi0 = X.gauss_point()
            sage: xi1 = X.point_from_discoid(t, 1)
            sage: p = Path(xi0, xi1)
            sage: p.derivative(t^2 + 2)
            2
            sage: p.derivative(t^2 + 2, 1/4)
            2
            sage: p.derivative(t^2 + 2, 1/2)
            0

        """
        assert f.denominator().degree() == 0, "f must be a polynomial."
        f = f.numerator()
        if s is None:
            s = self._s0
        xi = self.point(s)
        v = xi.pseudovaluation_on_polynomial_ring()
        phi = self._phi.numerator()
        psi, e = v.equivalence_decomposition(phi)[0]
        for g, m in v.equivalence_decomposition(f):
            if v(g-psi) > v(psi):
                return QQ(m/e)
        return QQ.zero()

    def point(self, s):
        r"""
        Return the point on this path corresponding to the parameter `s`.

        INPUT:

        - ``s`` -- a rational number in the interval `[s_0, s_1]`.

        OUTPUT:

        The point on the path corresponding to the parameter `s`.

        Since we use the parametrization given by a fixed key polynomial `\phi`,
        we have to find the MacLane valuation `v_s` such that `v_s(\phi) = s`.
        Here we use the functionality of MCLF to find `v_s`.

        EXAMPLES::

            sage: from sage.all import QQ, FunctionField
            sage: from mclf import BerkovichLine
            sage: from admissible_functions import Path
            sage: F.<t> = FunctionField(QQ)
            sage: X = BerkovichLine(F, QQ.valuation(2))
            sage: xi0 = X.gauss_point()
            sage: xi1 = X.point_from_discoid(t, 1)
            sage: p = Path(xi0, xi1)
            sage: p.point(1/2)
            Point of type II on Berkovich line, corresponding to v(t) >= 1/2

        """
        s0, s1 = self.boundaries()
        xi1 = self._xi1
        assert s0 <= s <= s1, "The parameter must be in the interval [s0, s1]."
        if s == s0:
            return self._xi0
        if s == s1:
            return xi1
        if not xi1.is_inductive():
            phi, s1_approx = xi1.approximation().discoid()
            while s1_approx <= s:
                phi, s1_approx = xi1.improved_approximation().discoid()
                self._s1_approx = s1_approx
                self._phi = phi
        else:
            phi, _ = xi1.discoid()
        return xi1.berkovich_line().point_from_discoid(phi, s)

# -----------------------------------------------------------------------------
# Helper functions


def effective_degree(f, xi):
    r""" Return the effective degree of a polynomial with respect to a point of type II.

    INPUT:

    - ``f`` -- a polynomial in K[t],
    - ``xi`` -- a point of type II on the Berkovich unit disk.

    OUTPUT:

    The effective degree of the polynomial `f` with respect to the point `xi`.

    EXAMPLES::

        sage: from sage.all import QQ
        sage: from mclf import BerkovichLine
        sage: from admissible_functions import effective_degree
        sage: F.<t> = FunctionField(QQ)
        sage: X = BerkovichLine(F, QQ.valuation(2))
        sage: f = t^2 + 1
        sage: xi = X.gauss_point()
        sage: effective_degree(f, xi)
        2

    """
    F = xi.berkovich_line().function_field()
    f = F(f)
    assert xi.type() == "II", "The point must be of type II."
    assert xi.is_in_unit_disk(), "The point must be on the unit disk."
    assert xi.parameter() == F.gen(), "The valuation must be induced by the standard parameter."
    v = xi.pseudovaluation_on_polynomial_ring()
    assert f.denominator().degree() == 0, "f must be a polynomial."
    return v.effective_degree(f.numerator())


def max_of_min(functions, s0, s1):
    r""" Return the maximum of the minimum of linear functions on an interval.

    INPUT:

    - ``functions`` -- a list of pairs `(a, b)` of rational numbers,
                       with `b\geq 0`, and at least one `b > 0`,
    - ``s0`` -- a rational number,
    - ``s1`` -- a rational number `>s_0`, or `Infinity`.

    OUTPUT:

    The pair `(m, s)` where `m` is the maximum of the minimum of the linear
    functions `a + b*s`, on the interval `[s0, s1]` and `s` is the lowest
    value of `s` where this maximum is attained.

    Special cases:
      1. If s1 == infinity and all b_i > 0, then f(t) tends to infinity as
         t goes to infinity and no finite t attains the supremum. In that
         case the function returns (Infinity, None).
      2. If s1 == infinity and some b_i are 0, then for large t the minimum
         is given by one of the constant functions (with b_i = 0). In that
         case, if m = min{a_i : b_i = 0},
         then f(t)=m for all sufficiently large t, and the smallest such t
         is computed.

    EXAMPLES::

        sage: from admissible_functions import max_of_min
        sage: max_of_min([(1, 1), (2, 1)], 0, 1)
        (2, 1)
        sage: max_of_min([(1, 1), (2, 1)], 0, 2)
        (3, 2)
        sage: max_of_min([(1, 1), (2, 1)], 0, Infinity)
        (+Infinity, None)
        sage: max_of_min([(1, 1), (2, 0)], 0, Infinity)
        (2, 1)

    """
    # Case 1: s1 is infinity (oo)
    if s1 == Infinity:
        # Check if any function is constant (i.e. has b == 0)
        constant_functions = [(a, b) for (a, b) in functions if b == 0]
        if not constant_functions:
            # All functions are strictly increasing so f(t) -> oo as t -> oo.
            return (Infinity, None)
        else:
            # The maximum is determined by the constant functions.
            m = min(a for (a, b) in constant_functions)
            t_candidates = [s0]
            for (a, b) in functions:
                if b > 0:
                    if a < m:
                        # Solve a + b*t = m  =>  t = (m - a)/b.
                        t_candidates.append((m - a) / b)
                    else:
                        # a is already at least m at t = s0.
                        t_candidates.append(s0)
            t_star = max(t_candidates)
            assert min(a + b * t_star for (a, b) in functions) == m, "Something went wrong."
            return (m, t_star)

    # Case 2: s1 is finite.
    # The candidate set will include the endpoints and all intersection points.
    candidates = [s0, s1]
    n = len(functions)
    for i in range(n):
        for j in range(i + 1, n):
            a_i, b_i = functions[i]
            a_j, b_j = functions[j]
            if b_i != b_j:
                # Find the intersection of f_i and f_j:
                # a_i + b_i*t = a_j + b_j*t   =>   t = (a_j - a_i) / (b_i - b_j)
                t_int = (a_j - a_i) / (b_i - b_j)
                if s0 <= t_int <= s1:
                    candidates.append(t_int)
    # Remove duplicates and sort the candidate t-values.
    candidates = sorted(set(candidates))

    # Now, evaluate f(t) = min_i (a_i + b_i*t) at each candidate,
    # and pick the candidate with the highest value.
    best_val = -Infinity
    best_t = None
    for t in candidates:
        current_val = min(a + b * t for (a, b) in functions)
        # Choose the candidate with a larger value;
        # in case of a tie, pick the smaller t.
        if current_val > best_val or (current_val == best_val and (best_t is None or t < best_t)):
            best_val = current_val
            best_t = t
    return (best_val, best_t)


def valuative_function(f, root):
    r""" Return the valuative function of a polynomial f.

    INPUT:

    - ``f`` -- a polynomial in K[t], or the factorization of such a polynomial,
    - ``root`` -- a point of type II on the Berkovich unit disk,
                  or a MacLane valuation on K[t]

    OUTPUT:

    A :class:`ValuativeFunction` instance, representing the valuative function
    of `f` on the closed diskoid `D` with root ``root``.

    EXAMPLES::

        sage: from sage.all import QQ
        sage: from mclf import BerkovichLine
        sage: from admissible_functions import valuative_function
        sage: F.<t> = FunctionField(QQ)
        sage: X = BerkovichLine(F, QQ.valuation(2))
        sage: xi0 = X.gauss_point()
        sage: f = 2*(t^2 + 2)*t
        sage: valuative_function(f, xi0)
        1 + 1*v(t) + 1*v(t^2 + 2)

        sage: xi = X.point_from_discoid(t-1, 1)
        sage: valuative_function(f, xi)
        1

    """
    if hasattr(root, "berkovich_line"):
        xi0 = root
        assert xi0.type() == "II", "the root must be a point of type II"
        assert xi0.is_in_unit_disk(), "the root must lie in the unit disk"
        F = xi0.berkovich_line().function_field()
    else:
        v0 = root
        assert v0.is_discrete_valuation(), ("the root is neither a point on "
                                            + "the Berkovich line nor a MacLane valuation")
        K = v0.domain().base_ring()
        F = FunctionField(K, v0.domain().variable_name())
        v_K = v0.restriction(K)
        X = BerkovichLine(F, v_K)
        xi0 = X.point_from_pseudovaluation_on_polynomial_ring(v0, F.gen())
    if isinstance(f, Polynomial):
        f = F(f)
        factorization = f.factor()
    elif f in F:
        factorization = f.factor()
    elif isinstance(f, Factorization):
        factorization = Factorization([(F(g), e) for g, e in f], unit=f.unit())
    else:
        raise ValueError(f"f= {f} is neither a polynomial nor a factorization")
    constant = xi0.v(factorization.unit())
    coefficients = []
    functions = []
    for g, e in factorization:
        if effective_degree(g, xi0) > 0:
            coefficients.append(e)
            functions.append(g)
        else:
            constant += e*xi0.v(g)
    return ValuativeFunction(constant, coefficients, functions, root=xi0)


def linear_sum_of_valuative_functions(coefficients, functions):
    r"""
    Return the linear form a_0 + a_1*v(f_1) + ... + a_n*v(f_n).

    INPUT:

    - ``coefficients`` -- a list `[a_0, a_1, ..., a_n]` of positive rational numbers,
    - ``functions`` -- a list of valuative functions, defined on the same closed diskoid.

    OUTPUT:

    A :class:`ValuativeFunction` instance
    """
    raise NotImplementedError()
