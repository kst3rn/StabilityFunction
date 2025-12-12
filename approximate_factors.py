r"""
Approximate factors of polynomials over (fake) p-adic number fields
===================================================================

Let `K` be a number field and `v_K` a nontrivial discrete valuation on `K`.
Let `hat{K}` denote the completion of `K` with respect to `v_K`.

A polynomial `f\in K[x]` is called *strongly irreducible* if it is
irreducible over `\hat{K}`. Given `f\in K[x]`, a *strong factor* of `f`
is a monic irreducible factor `g\in\hat{K}[x]` of `f`.

This module provides functionality for computing arbitrarily precise
approximations of all strongly irreducible factors of a given polynomial
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
a strong factor `g` of `f` such that `v\leq v_g`

EXAMPLES:

    sage: R.<x> = QQ[]
    sage: v_2 = QQ.valuation(2)
    sage: f = x^6 + 4*x + 1
    sage: F = approximate_factorization(f, v_2); F
    [approximate prime factor of x^6 + 4*x + 1 of degree 2,
     approximate prime factor of x^6 + 4*x + 1 of degree 4]
    
    sage: g1 = F[0]
    sage: g1.approximate_factor()
    x + 1

    sage: g1.approximate_factor(5)
    x^2 + 12*x - 15

"""

from sage.all import SageObject, GaussValuation, Infinity
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

    a list of the irreducible factors of `f` over the completion of `K` at `v_K`.
    The irreducible factors are objects of the class :class:`ApproximateFactor`.

    If `g_0` is given, then only the factors approximated by `g_0` are returned. 

    NOTE: for the moment, only irreducible factors which are integral with respect
    to `v_K` are returned.

    """
    assert v_K.domain() == f.base_ring(), "the domain of v_K must be the base ring of f"
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
        # we check whether this factor is really irreducible
        # as a side effect, we compute the next improved approximation
        # note that phi has the correct degree, but v.phi() may not 
        F = v.equivalence_decomposition(f, compute_unit=False)
        assert len(F) == 1 and F[0][1] == 1, "this factor is not irreducible"
        phi = F[0][0]
        self._phi = phi
        self._degree = phi.degree()
        if phi.degree() == f.degree():
            self._prec = Infinity
        else:
            self._prec = v(v.phi())

    def __repr__(self):
        return f"approximate prime factor of {self._polynomial()} of degree {self.degree()}"
    
    def polynomial(self):
        return self._polynomial
    
    def valuation(self):
        return self._valuation
    
    def degree(self):
        r""" Return the degree of this approximate factor.
        
        """
        return self._degree
    
    def prec(self):
        return self._prec
    
    def improve_approximation(self):
        r""" Improve the approximation of the approximate prime factor.

        This function has no output, but it replaces the underlying inductive valuation
        by the next better approximation given by one MacLane step.
        
        """
        if self.prec() < Infinity: 
            v0 = self.valuation()
            f = self.polynomial()
            phi = v0.equivalence_decomposition(f)[0][0]
            f1, f0 = f.quo_rem(phi)
            t = v0(f0) - v0(f1)
            v1 = v0.augmentation(phi, t)
            F = v1.equivalence_decomposition(f)
            assert len(F) == 1 and F[0][1] == 1, "something is wrong"
            self._valuation = v1
            self._prec = t

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
            while self._prec < prec:
                self.improve_approximation()
            return self.valuation().phi()