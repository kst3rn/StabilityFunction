r"""
Approximate factors of polynomials over (fake) p-adic number fields
===================================================================

Let `K` be a number field and `v_K` a nontrivial discrete valuation on `K`.
Let `\hat{K}` denote the completion of `K` with respect to `v_K`.

A polynomial `f\in K[x]` is called *strongly irreducible* if it is
irreducible over `\hat{K}`. Given `f\in K[x]`, a *strong prime factor* of `f`
is a monic irreducible factor `g\in\hat{K}[x]` of `f`.

This module provides functionality for computing arbitrarily precise
approximations of all strong prime factors of a given polynomial
`f\in K[x]`. It uses, in an essential way, the theory of MacLane on
inductive valuations which is already implemented in Sage.


A *MacLane* valuation is a discrete pseudovaluation `v` on `K[x]` extending `v_K`
and such that `v(x)\geq 0`. Recall that the set `V_K` of all MacLane valuations 
forms a partially orderd set, and that this poset is a *rooted tree*. Its least
element of the *Gauss valuation*.

Given a monic, integral (with respect to `v_K`) and strongly irreducible 
polynomial `g\in K[x]` there exists a unique MacLane pseudovaluation `v_g`
such that `v_g(g)=\infty`. This is a maximal element of `V_K`.
    
A MacLane valuation `v` is called an *approximate factor* of `f` if there exists
a strong prime factor `g` of `f` such that `v\leq v_g`. It is called an 
*approximate prime factor* if there is a unique such factor `g`, and then it
is called an *approximation* of `g`.

In our implementation, a strong prime factor `g` of `f` is represented by one of
its approximations. Such an approximation can be improved, using MacLane's method,
and any sequence of improvements will converge towards the 'true' prime factor `g`.


EXAMPLES:

    sage: R.<x> = QQ[]
    sage: v_2 = QQ.valuation(2)
    sage: f = x^6 + 4*x + 1
    sage: F = approximate_factorization(f, v_2); F
    [approximate prime factor of x^6 + 4*x + 1 of degree 2,
     approximate prime factor of x^6 + 4*x + 1 of degree 4]
    
    sage: g = F[0]
    sage: g.approximate_factor()
    x + 1

We see that the degree of the prime factor is not necessary equal
to its first approximation. But it is after one improvement:

    sage: g.improve_approximation()
    sage: g.approximate_factor()
    x^2 +1 

You can also force the approximation to have a guaranteed 
precision:

    sage: g.approximate_factor(10)
    x^2 - 28/11*x + 7/23

The *precision* of this approximation beeing `\geq N` means that
the valuation of the difference between a root of the approximation
and the nearest root of the true prime factor is `\geq N`. This is
of the same rough magnitude as the valuation of the approximation,
but may not be equal to it: 

    sage: g.precision()
    21/2

    sage: g.valuation()(g.approximate_factor())
    23/2

"""

from sage.all import SageObject, GaussValuation, Infinity, PolynomialRing, QQ, ZZ
from sage.geometry.newton_polygon import NewtonPolygon
from sage.rings.valuation.limit_valuation import LimitValuation 


def approximate_factorization(f, v_K, g0=None, assume_squarefree=False, 
                               assume_irreducible=False):
    r""" Return the factorization of a polynomial over a p-adic number field.
    
    INPUT:

    - ``f`` -- a nonconstant polynomial over a number field `K`
    - ``v_K`` - a nontrivial discrete valuation on `K`
    - `g0` -- an approximate factor of `f`, or ``None``
    - `assume_squarefree` --  a boolean (default: ``False``)
    - `assume_irreducible` -- a boolean (default: ``False``)

    OUTPUT:

    A list of the approximate prime factors of `f` with respect to `v_K`.
    These are objects of the class :class:`ApproximatePrimeFactor`.

    If `g_0` is given, then only the factors approximated by `g_0` are returned. 

    NOTE: for the moment, only irreducible factors which are integral with respect
    to `v_K` are returned.

    
    TODO:

    - write examples

    """
    f = f.change_ring(v_K.domain())
    assert f.degree() > 0, "f must be nonconstant"
    if not assume_squarefree:
        f = f.radical() 
    if not assume_irreducible:
        ret = []
        for g, _ in f.factor():
            ret += approximate_factorization(g, v_K, g0=g0, 
                                              assume_squarefree=True, 
                                              assume_irreducible=True)
        return ret
    # make f integral and primitive
    m = min(v_K(c) for c in f.coefficients())
    f = v_K.element_with_valuation(-m)*f
    if g0 is None:
        v0 = GaussValuation(f.parent(), v_K)
        g0 = ApproximateFactor(f, v0)
    if g0.is_irreducible():
        return [ApproximatePrimeFactor(f, g0.valuation())]
    ret = []
    for g in g0.mac_lane_step():
        ret += approximate_factorization(f, v_K, g0=g, 
                                          assume_squarefree=True, 
                                          assume_irreducible=True)
    return ret


def approximate_roots(f, v_K, positive_valuation=True):
    r"""
    Return approximate roots of an univariate polynomial.

    INPUT:

    - ``f`` -- a nonconstant polynomial over a number field `K`
    - ``v_K`` -- a `p`-adic valuation on `K`
    - ``positive_valuation`` -- boolean (default: ``True``); if set, 
                                require `v_L(a) > 0` for all approximations

    OUTPUT:

    a list of instances of :class:`ApproximateRoot`.

    The returned objects represent approximate roots `a in L` (for some finite
    extension `(L,v_L)/(K,v_K)`) whose residual valuations can be made arbitrarily large
    by repeated calls to :meth:`ApproximateRoot.improve_approximation`.

    TODO:

    - write examples

    """
    factors = approximate_factorization(f, v_K)
    if positive_valuation:
        return [ApproximateRoot(g) for g in factors if g.value() > 0]
    else:
        # the test g.value() >= 0 is superfluous because the current code
        # only produces factors with this property; included for clarity 
        # and because this may change in the future
        return [ApproximateRoot(g) for g in factors if g.value()>=0]


class ApproximateFactor(SageObject):
    r""" An approximate factor of a polynomial over a p-adic number field.

    INPUT:

    - ``f`` -- a nonconstant polynomial over a number field `K`
    - ``v`` -- a discrete valuation on the polynomial ring
               to which `f` belongs

    It is assumed that `v` is an *approximate factor* of `f`, i.e. that
    `v` is a MacLane valuation on `K[x]` and that there exists a strong
    factor `g` of `f` such that `v\leq v_g`.
    
    OUTPUT:

    an object representing this approximate factor.

    This is the base class for :class:`ApproximatePrimeFactor`. It is only used
    in the process of finding the factorization of an irreducible polynomial
    into its strong prime factors, as approximations of factors which may not
    be prime.

    The method :meth:`mac_lane_step`, produces an improved approximation, which
    may consists of several approximate factors. 
      
    """
    def __init__(self, f, v):
        v_K = v._base_valuation
        assert v.domain() == f.parent(), "the domain of v must be the parent of f"
        self._polynomial = f
        self._valuation = v
        self._base_valuation = v_K
        F = v.equivalence_decomposition(f, compute_unit=False)
        self._degree = sum(e*phi.degree() for phi, e in F)
        self._is_irreducible = (len(F) == 1 and F[0][1] == 1)
        self._equivalence_decomposition = F

    def __repr__(self):
        return f"approximate factor of {self.polynomial()} of degree {self.degree()}"
    
    def base_valuation(self):
        return self._base_valuation
    
    def base_field(self):
        return self.base_valuation().domain()
    
    def polynomial(self):
        return self._polynomial
    
    def valuation(self):
        return self._valuation
    
    def degree(self):
        r""" Return the degree of this approximate factor.
        
        """
        return self._degree

    def is_irreducible(self):
        r""" Return whether this approximate factor is irreducible.
        
        """
        return self._is_irreducible

    def equivalence_decomposition(self):
        return self._equivalence_decomposition
    
    def mac_lane_step(self):
        r""" Return a list of approximate factors which refine this factor.
        
        """
        assert not self.is_irreducible(), "The MacLane step only makes sense if this factor is not irreducible"
        v0 = self.valuation()
        f = self.polynomial()
        ret = []
        for phi, _ in self.equivalence_decomposition():
            t0 = v0(phi)
            v00 = v0.augmentation(phi, t0, check=False)
            # we use v00 only for convenience; it is important not to augment it
            valuations = list(v00.valuations(f))
            # the values of v_00 on the terms of the phi-expansion of f
            a = min(valuations)
            n = min(i for i in range(len(valuations)) if valuations[i] == a)
            # n is the "degree" of f with respect to v00; this means that n*deg(phi)
            # is the number of roots of f inside the residue class corresponding to phi
            assert n > 0, "something is wrong!"
            # we find the maximal value t1 > t0 such that the discoid corresponding to
            # v1=[v0, v(phi)=t1] still contains all the prime factors of f in
            # the residue class corresponding to phi:
            np = NewtonPolygon(enumerate(valuations))
            slopes = [mu for mu in np.slopes(repetition=False) if mu < 0]
            if len(slopes) == 0:
                t1 = Infinity
            else:
                t1 = t0 - max(slopes)
            v1 = v0.augmentation(phi, t1)
            ret.append(ApproximateFactor(self.polynomial(), v1))
        return ret
    

class ApproximatePrimeFactor(ApproximateFactor):
    r""" An approximate prime factor of a polynomial over a p-adic number field.

    INPUT:

    - ``f`` -- a nonconstant polynomial over a number field `K`
    - ``v`` -- a discrete valuation on the polynomial ring
               to which `f` belongs

    It is assumed that `v` is an *approximate prime factor* of `f`, i.e. that
    `v` is a MacLane valuation on `K[x]` and that there exists a *unique* strong
    factor `g` of `f` such that `v\leq v_g`.
    
    OUTPUT:

    an object representing this approximate prime factor.
      
    """
    def __init__(self, f, v):
        K = f.base_ring()
        v_K = v.restriction(K)
        assert v.domain() == f.parent(), "the domain of v must be the parent of f"
        self._polynomial = f
        self._valuation = v
        self._base_valuation = v_K
        self._value = v(f.parent().gen())

        # we check whether this factor is really irreducible
        F = v.equivalence_decomposition(f, compute_unit=False)
        assert len(F) == 1 and F[0][1] == 1, "this factor is not irreducible"
        phi = F[0][0]
        self._degree = phi.degree()
        if phi.degree() == f.degree():
            self._valuation = v.augmentation(f, Infinity)
            self._precision = Infinity
        else:
            R = f.parent()
            S = PolynomialRing(R, "T")
            x = R.gen()
            T = S.gen()
            F = f(x+T)
            self._F = [F[i] for i in range(phi.degree() + 1)]
            # these two attributes have to updated after every improvement step,
            # in this order:
            self._v_g = v.augmentation(v.phi(), Infinity)
            self._compute_precision()

    def __repr__(self):
        return f"approximate prime factor of {self.polynomial()} of degree {self.degree()}"
    
    def _compute_precision(self):
        r""" Compute and store the current precision of this approximate prime factor.
        
        The *precision* of this 
        """
        v_g = self._v_g
        F = self._F
        self._precision = max((v_g(F[0]) - v_g(F[i]))/i for i in range(1, self.degree() + 1))

    def polynomial(self):
        r""" Return the irreducible polynomial of which this is an approximate prime factor.
        """
        return self._polynomial
    
    def valuation(self):
        r""" Return the inductive valuation underlying this approximate prime factor.
        """
        return self._valuation
    
    def degree(self):
        r""" Return the degree of this approximate factor.
        
        """
        return self._degree
    
    def value(self):
        r""" Return the value of this approximate prime factor.

        The *value* of this approximate prime factor is defined as
        `v(x)`, where `v` is the inductive valuation represeting it.
        It is equal to the valuation of a root `\alpha` of this factor,
        in a suitable extension of the completed base field.
 
        """
        return self._value
    
    def precision(self):
        r""" Return the precision of this approximate prime factor.
        
        The *precision* of this approximate prime factor is the valuation
        `v_L(\alpha-\alpha_0)`, where `\alpha` is a root of this factor,
        and `\alpha_0` is a root of the current approximation, closest to
        `\alpha`. So it measures the accurary with which the root `\alpha`
        is 'known' by the current approximation. 
        
        """
        return self._precision
    
    def improve_approximation(self):
        r""" Improve the approximation of the approximate prime factor.

        This function has no output, but it replaces the underlying inductive valuation
        by the next better approximation given by one MacLane step.
        
        """
        if self.precision() < Infinity: 
            v0 = self.valuation()
            f = self.polynomial()
            g = v0.equivalence_decomposition(f)[0][0]
            f1, f0 = f.quo_rem(g)
            t = v0(f0) - v0(f1)
            v1 = v0.augmentation(g, t)
            v_g = v0.augmentation(g, Infinity)
            self._valuation = v1
            self._v_g = v_g
            self._compute_precision()

    def approximate_factor(self, prec=None):
        r""" Return the current approximation of this approximate prime factor.

        INPUT:

        - ``prec`` -- a nonnegative rational number (default: ``None``)

        OUTPUT:

        An approximation of this approximate prime factor with precision
        at least ``prec``. If ``prec`` is not given, the current approximation
        is returned.
          
        """
        if prec is None:
            return self.valuation().phi()
        else:
            while self.precision() < prec:
                self.improve_approximation()
            return self.valuation().phi()
        

class ApproximateRoot(SageObject):
    r"""
    An object representing an *approximate root* of an univariate polynomial.

    INPUT:

    - ``g`` -- an instance of :class:`ApproximatePrimeFactor`, representing
               an irreducible factor of a polynomial `f` over the completion
               `\hat{K}` of its base field `K`
    - ``name`` -- an alphanumerical string (default: "alpha")

    OUTPUT:

    an object representing the tautological root `a\in \hat{L}`, where 

    .. MATH::

        \hat{L} := \hat{K}[x]/(g)

    Internally, an approximate root `a` is represented by 

    - an approximate prime factor `g` of `f`
    - a fixed approximation `g_0` of `g`.

    The approximation `g_0` has the property that the stem field 

    .. MATH::

        L := K(a_0) = K[x]/(g_0)

    of `g_0` is stable, after completion and up to isomorphism, for further 
    approximations of `g`. By Krasner's Lemma this condition is satisfied 
    for every sufficiently precise approximation `g_0`. 

    The element `a_0\in L`, the tautological root of `g_0`, is then the
    first approximation of the approximate root `a`. 

    Calling :meth:`improve_approximation` computes and returns an improved 
    approximation.

    """
    def __init__(self, g, name="alpha"):
        self._is_exact = (g.precision() == Infinity)
        self._g = g
        f = g.polynomial()
        self._f = f
        self._df = f.derivative()
        # 1. improve approximation until Hensel's lemma applies
        self._force_hensel()
        # 2. construct L
        K = self.base_field()
        g0 = self.prime_factor().approximate_factor()
        L = K.extension(g0, name)
        v_L = self.base_valuation().extension(L)
        self._extension = L
        self._extension_valuation = v_L
        # 3. initialize first approximation
        self._approximation = L.gen()
        self._precision = g.precision()
        # 4. construct limit valuation
        self._init_limit_valuation()

    def __repr__(self):
        return f"approximate root of {self.polynomial()}"

    def base_field(self):
        r""" Return the base field of this approximate root.
        """
        return self.prime_factor().base_field()
    
    def base_valuation(self):
        r""" Return the valuation of the base field of this approximate root.
        """
        return self.prime_factor().base_valuation()

    def prime_factor(self):
        r""" Return the strong prime factor of which ``self`` is a root. 
        """
        return self._g

    def extension(self):
        r"""Return the field extension in which this approximate root lives."""
        return self._extension

    def extension_valuation(self):
        r"""Return the valuation on the extension field."""
        return self._extension_valuation
    
    def polynomial_ring(self):
        return self.polynomial().parent()
    
    def polynomial(self):
        r""" Return the polynomial of which ``self`` is a root.
        """
        return self._f
    
    def derivative(self):
        return self._df
    
    def value(self):
        r""" Return the valuation of this root.
        """
        return self.prime_factor().value()

    def approximation(self):
        r"""Return the current approximation of this root."""
        return self._approximation
    
    def eval(self, r, simplify=True):
        r""" Evaluate a polynomial on this approximate root.
        
        INPUT:

        - ``r`` -- a polynomial over the base field `K`
        
        OUTPUT:

        the value `r(a_0)` where `a_0` is the current approximation 
        of this approximate root.

        If ``simplify`` is ``True`` then the value is simplified to the 
        current precision before it is returned.

        """
        v_L = self.extension_valuation()
        a0 = self.approximation()
        if simplify:
            return v_L.simplify(r(a0), self.precision())
        else:
            return r(a0)
    
    def precision(self):
        r""" Return a lower bound for the precision of the current approximation.
        
        If `t` is return, then `v_L(a-a_0) >= t`, where `a_0` is the 
        current approximation and `a` is the exact root. 
        """
        if self.is_exact():
            return Infinity
        return self._precision
    
    def is_exact(self):
        r""" Return whether this approximate root is exact.
        
        """
        return self._is_exact
    
    def value_of_poly(self, r):
        r"""Return the valuation of a polynomial in this root.
        
        INPUT:

        - ``r`` -- a univariate polynomial over the base field

        OUTPUT:

        the valuation `v_L(r(a))`, where `a\in \hat{L}` is the root
        of `f` represented by ``self``.

        ALGORITHM:

        the *actual root* `a\in\hat{L}` can be realized as a 
        *limit valuation* on `K[x]`, i.e. the discrete pseudovaluation
        `v` defined by `v(r) = v_L(r(a))`. Note that `v(f)=\infty`,
        where `f` is the minimal polynomial of `a` over `K`. Therefore,
        `v` is only a pseudovaluation. 

        There is a native implementation :class:`LimitValuation` in 
        Sage; however, it expects the equation `f` to be monic and integral
        with respect to `v_K`. For this reason we have first rescale
        `v` before it can be implemented, i.e. replace `f` by

        .. MATH::

            f_1 := c^{-1}\pi^{ad}\cdot f(\pi^{-a}x),

        where `d` is the degree of `f` and `c` its leading coefficient.

        """
        try:
            r = self.polynomial_ring()(r)
        except (TypeError, ValueError):
            raise ValueError(
                f"r must be coercible into {self.polynomial_ring()} (got {r.parent()})")
        return self._limit_valuation(r)

    def improve_approximation(self):
        r"""
        Improve the current approximation.

        OUTPUT:

        The improved approximation of this root.

        """
        if self.is_exact():
            # there is nothing to improve
            return
        # we can assume that the approximation can be improve using Hensel's Lemma, 
        # i.e. Newton approximation
        f = self.polynomial()
        df = self.derivative()
        v_L = self.extension_valuation()
        a0 = self._approximation
        a = a0 - f(a0)/df(a0)
        # compute the new precision
        self._precision = v_L(f(a)) - v_L(df(a))
        # simplify a to new precision
        a = self.extension_valuation().simplify(a, self.precision()+1)
        self._approximation = a
        self._count = self._count + 1
        return a
    
    def _force_hensel(self):
        r""" Improve the approximation of the prime factor `g` so that
        approximations of the root can be improved using Hensel's Lemma.
        
        This guarantees that the stem field `\hat{L}=\hat{K}[x]/(g_0)`, 
        where `g_0` is the current approximation, contains the root of `f`
        represented by ``self``.
        """
        
        if self.is_exact():
            self._precision = Infinity
            return
        g = self.prime_factor()
        g.improve_approximation()
        f = self.polynomial()
        df = self.derivative()
        m = g.valuation()(f)
        s = g.valuation()(df)
        while m <= 2*s + 1:
            g.improve_approximation()
            m = g.valuation()(f)
            s = g.valuation()(df)
        self._count = 0
        self._precision = m - s
    
    def _init_limit_valuation(self):
        r""" Construct the limit valuation corresponding to this approximate root,
        i.e. the discrete valuation v on `K[x]` such that
        
        .. MATH::

            v(r) = v_L(r(a)),

        for all polynomials `r\in K[x]`, and where `a\in\hat{L}` is the exact root
        represented by this approximate root. 
        """
        # in principal we can use LimitValuation(v, f), where v is the MacLane valuation
        # corresponding to the approximate factor.
        # the problem is that this only works if f is monic and integral
        # which we can't assume. 
        # We have to force this by scaling, i.e. replace f by
        # f(c^(-1)*x).monic() and v by the corresponding `scaled valuation' 
        f = self.polynomial()
        v_K = self.base_valuation()
        # we compute the maximal slope of the NP of f
        d = f.degree()
        a_d = v_K(f[d])
        m = max((a_d-v_K(f[i])) for i in range(d))
        # if m <= 0, f is already integral
        if m > 0:
            c = v_K.element_with_valuation(-m)
        else:
            m = QQ.zero()
            c = v_K.domain().one()
        self._scaling_factor = c
        f1 = f(c*f.parent().gen()).monic()
        v1 = _scale_inductive_valuation(self.prime_factor().valuation(), v_K, m)
        self._scaled_limit_valuation = LimitValuation(v1, f1)
    
    def _limit_valuation(self, r):
        x = r.parent().gen()
        c = self._scaling_factor
        v = self._scaled_limit_valuation
        return v(r(c*x))


def _scale_inductive_valuation(v, v_K, m):
    r""" Scale a MacLane valuation. 
    
    INPUT:

    - ``v`` -- an inductive valuation on a polynomial ring `K[x]`
    - ``v_K`` -- the base valuation, i.e. the restriction of v to `K`
    - ``m`` -- a nonnegative value of `v_K`

    OUTPUT:

    the scaling of `v` by `m`; this is the inductive valuation `v_1:=v\circ\tau`,
    where `\tau:K[x]\to K[x]` is the ring automorphism defined by
    `\tau(x) = c\cdot x` for an element `c\in K` with `v_K(c)=m`.

    EXAMPLES:

    sage: R.<x> = QQ[]
    sage: v_K = QQ.valuation(5)
    sage: f = x^9 + 2*x^8 + x^7 + 4*x^6 + 2*x^5 + x^4 + 3*x^3 + 2
    sage: v0 = GaussValuation(R, v_K)
    sage: v = v0.augmentation(f, 1)
    sage: _scale_inductive_valuation(v, v_K, 2)
    [ Gauss valuation induced by 5-adic valuation, v(x) = 2, v(x^9 + ...) = 19 ]

    """
    ci = v_K.element_with_valuation(-m)  # inverse of the scaling element
    x = v.domain().gen()
    if v.is_gauss_valuation():
        if m == 0:
            return v
        else:
            return v.augmentation(x, m) 
    v0 = _scale_inductive_valuation(v.augmentation_chain()[1], v_K, m)
    phi = v.phi()
    l = v(phi)
    # this may give an error, because phi(ci*x).monic may not be equivalence irreducible!
    #return v0.augmentation(phi(ci*x).monic(), l + phi.degree()*m)
    return _augmentation(v0, phi(ci*x).monic(), l + phi.degree()*m)


def _augmentation(v0, phi, s):
    r""" Helper function for _scale_inductive_valution.
    
    INPUT:

    - ``v_0`` -- an inductive valuation on a polynomial ring `K[x]`
    - ``phi`` -- a monic, integral and absolutely irreducible element of `K[x]`
    - ``s`` -- a rational number 

    It is assumed that `\phi` is not an equivalence unit for `v_0`, and that
    `v_0(\phi) < s`.

    OUTPUT:

    the unique inductive valuation of the form

    .. MATH::

         v = [v_0,\ldots, v(\phi)=s].

    """
    F = v0.equivalence_decomposition(phi)
    if len(F) != 1:
        raise ValueError(f"phi = {phi} is either not absolutely irreducible"+
                         f"or an equivalence unit for v0 = {v0}")
    while not v0.is_key(phi):
        v1 = v0.mac_lane_step(phi)[0]
        assert v1(phi) < s, "v0= {v0}, v1 = {v1}, phi = {phi}"
        v0 = v1
    return v0.augmentation(phi, s)


# ---------------------------------------------------------------------------------------

#     Tests


def test_precision(g):
    print(f"value = {g.value()}")
    for _ in range(10):
        g0 = g.approximate_factor()
        print(f"g0 = {g0}")
        print(f"v(g_0) = {g.valuation()(g0)}")
        print(f"prec = {g.precision()}")
        print()
        g.improve_approximation()


def test_approximate_roots(R, v_K, N=10, d = 12):
    for _ in range(N):
        f = R.random_element(d)
        print()
        print(f"f = {f}")
        roots = approximate_roots(f, v_K)
        for a in roots:
            print(f"   a_0 = {a.approximation()}")
            for _ in range(5):
                a.improve_approximation() 
            print(f"a_5 = {a.approximation()}")           
            print(f"precision = {a.precision()}")    
            print(f"v_L(f(a_5)) = {a.extension_valuation()(f(a.approximation()))}")
            h = R.random_element()
            print(f"h= {h}")
            assert a.extension_valuation()(h(a.approximation())) <= a.value_of_poly(h)


def test_scale_valuation(N=10, p=2):
    R = PolynomialRing(QQ, "x")
    x = R.gen()
    for _ in range(N):
         f = R.random_element(10)
         print(f"f = {f}")
         G = approximate_factorization(f, QQ.valuation(p))
         for g in G:
             v = g.valuation()
             print(f"v= {v}")
             s = ZZ.random_element().abs()
             print(f"s = {s}")
             v_s = _scale_inductive_valuation(v, QQ.valuation(p), s)
             print(f"v_s = {v_s}")
             h = R.random_element()
             assert v(h) == v_s(h(x/p**s)), "h = {h}, v(h) = {v(h)}, v_s(h_1) = {v_s(h(x/p**s))}"