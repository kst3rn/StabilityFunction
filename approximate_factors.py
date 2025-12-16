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

from sage.all import SageObject, GaussValuation, Infinity, PolynomialRing
from sage.geometry.newton_polygon import NewtonPolygon


def approximate_factorization(f, v_K, g0=None, assume_squarefree=False, 
                               assume_irreducible=False):
    r""" Return the factorization of a polynomial over a p-adic number field.
    
    INPUT:

    - ``f`` -- a nonconstant polynomial over a number field `k`
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
        return f"approximate factor of {self._polynomial()} of degree {self.degree()}"
    
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
        v_K = v._base_valuation
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
        return f"approximate prime factor of {self._polynomial()} of degree {self.degree()}"
    
    def _compute_precision(self):
        r""" Compute and store the current precision of this approximate prime factor.
        
        The *precision* of this 
        """
        v_g = self._v_g
        F = self._F
        self._precision = max((v_g(F[0]) - v_g(F[i]))/i for i in range(1, self.degree() + 1))

    def polynomial(self):
        r""" Return the irreducible polynomial of which this is aan approximate prime factor.
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