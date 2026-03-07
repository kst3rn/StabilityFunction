r"""
Stable reduction of plane quartics
==================================

This module implements the computation of the *stable reduction* of a
smooth plane quartic over a p-adic number field, restricted to the case of
*non-hyperelliptic reduction*.

Given a homogeneous quartic form `F` over a number field `K` and a p-adic
valuation `v_K` on `K`, the main entry point
:func:`stable_reduction_of_quartic` carries out the following steps:

1. Construct a GIT-semistable model of the plane quartic with respect to `v_K`,
   using the machinery provided by
   :class:`~semistable_model.curves.plane_curves_valued.PlaneCurveOverValuedField`.

2. Compute the special fibre of this model over a finite extension of the residue
   field, and test whether the resulting curve is stable.  If the special fibre
   is not stable, the reduction is hyperelliptic and the computation stops.

3. Resolve each cusp of the special fibre and determine whether the corresponding
   exceptional divisor is an elliptic tail or a pigtail.

4. Construct the component graph of the stable curve using
   :func:`component_graph_of_GIT_stable_quartic`, and classify its combinatorial
   type via canonical signatures.

The output is packaged in a lightweight container class
:class:`StableReductionResult`, which records the most important intermediate
objects (GIT model, special fibre, component graph) together with the final
reduction type.  The class is designed to be easy to inspect interactively and
to allow later serialization for large-scale experimental computations.

At present, only the non-hyperelliptic case is handled; if hyperelliptic
reduction is detected, this is recorded in the result object.

.. todo::

    - fix the error with message "there is no unique extension of 3-adic valuation 
      from Rational Field ...", see the examples below
    - Improve the performance of :func:`semistable_model.curves.approximate_factors.ApproximateRoot'
      by simplifying the coefficients of a polynomial before evaluating it on the
      approximate root.

"""

import hashlib
from sage.all import QQ, GF
from semistable_model.curves.plane_curves import ProjectivePlaneCurve
from semistable_model.curves.plane_curves_valued import PlaneCurveOverValuedField
from semistable_model.curves.cusp_resolution import resolve_cusp
from semistable_model.curves.component_graphs_of_plane_curves import component_graph_of_GIT_stable_quartic
from semistable_model.curves.genus3_reduction_types import classify_genus3_type


def stable_reduction_of_quartic(F, v_K, compute_matrix=False):
    r""" Return the stable reduction of a smooth quartic.
    
    INPUT:

    - ``F`` -- a form of degree `4` in three variables over a number field `K`
    - ``v_K`` -- a p-adic valuation on `K`

    It is assume that `F` defines a smooth plane quartic `X`.

    OUTPUT:

    An instance of :class:`StableReductionResult` describing the geometric stable
    reduction of ``X`` with respect to ``v_K``.

    The returned object has one of the following statuses:

    - ``status == "ok"``:
      The special fibre of a GIT-semistable model of ``X`` is stable and
      non-hyperelliptic.  In this case the result object contains, among other
      data:
      
        * the GIT-semistable model and its special fibre,
        * the component graph of the stable curve,
        * the canonical signature of this graph,
        * the reduction type as a string (one of the 42 genus-3 types).

    - ``status == "hyperelliptic"``:
      The special fibre is not stable, which in genus 3 corresponds to
      hyperelliptic reduction.  No component graph is produced in this case.

    - ``status == "fail"``:
      An exception occurred during the computation.  The exception message is
      recorded in the ``warnings`` attribute of the result object.

    If ``compute_matrix`` is ``True`` then for each cusp/tail we record a triple
    `(v_L, t, E, T)`, where `v_L` is the valuation on the field extension over
    the cusp could be resolved, `t` is the thickness of the node, `E` is the tail 
    (a semistable plane cubic) and `T` is the base change matrix representing 
    the resolution.

    If ``compute_matrix` is ``False`` (the default) then instead a triple 
    `(v_L, t, E, e)` is recorded, where `t`, `E` are as before, `v_L` is the 
    valuation on a subextension and `e` is a positive integer indicating the 
    ramification index necessary to obtain the actual extension.

    ALGORITHM (overview):

    1. Construct a GIT-semistable model of ``X`` with respect to ``v_K`` and
       compute its special fibre over a finite extension of the residue field.

    2. Test whether the special fibre is stable.  If not, record hyperelliptic
       reduction and stop.

    3. For each cusp of the special fibre, compute the resolution and determine
       whether the corresponding exceptional divisor is an elliptic tail or a
       pigtail.

    4. Build the component graph of the resulting stable curve and classify its
       combinatorial type via canonical signatures.

    NOTES:

    - Only the non-hyperelliptic case is handled at present.
    - The function is designed for expensive computations; the result object is
      intended to be inspected interactively or stored for later analysis.

    data describing the stable reduction of `X` with respect to `v_K`, 
    if `X` does not have hyperelliptic reduction. If it has hyperelliptic 
    reduction then ``None`` is returned.

    EXAMPLES:

    The following example produced an error in a previous version:

        sage: from semistable_model.curves.stable_reduction_of_quartics import stable_reduction_of_quartic
        sage: R.<x,y,z> = QQ[]
        sage: F = -x^3*y - 8*x*y^3 - 7*x^3*z - 7*x^2*y*z + 5*x*y^2*z + 6*x^2*z^2 - 6*y*z^3
        sage: SR = stable_reduction_of_quartic(F, QQ.valuation(2)); SR  # long time
        StableReductionResult(ok, type=0---0e)

        sage: SR.git_extension
        Number Field in a1 with defining polynomial x^4 + 2*x^2 + 4

        sage: SR.tail_data
        {(z2 : z2 : 1): (2-adic valuation,
         Projective Plane Curve over Finite Field in u2 of size 2^2 defined by x^3 + u2*y^2*z + (u2 + 1)*y*z^2,
         1)}
        
        sage: SR.to_json_record()
        {'id': 'c2dac3579f4dd19d',
         'K_defpoly': 'QQ',
         'K_gen': None,
         'p': 2,
         'F': '-x^3*y - 8*x*y^3 - 7*x^3*z - 7*x^2*y*z + 5*x*y^2*z + 6*x^2*z^2 - 6*y*z^3',
         'reduction_type': '0---0e',
         'git_dfe': [4, 2, 2],
         'cusp_dfes': [[8, 1, 8]]}         
    
    This example also caused an error, see issue #67

        sage: F = 4*x^4 - 2*x*y^3 + 6*x^3*z + x*y*z^2 + 10*y^2*z^2 + 7*x*z^3 - 7*y*z^3 + 9*z^4
        sage: SR = stable_reduction_of_quartic(F, QQ.valuation(3)); SR
        StableReductionResult(ok, type=0---0e)

    The following example shows that issue #65 has been fixed:

        sage: F = 3*x^4 - 12*x^2*y^2 - 5*y^4 - 5*x^3*z - 8*x^2*y*z + 3*x*y^2*z + 18*y*z^3
        sage: SR = stable_reduction_of_quartic(F, QQ.valuation(3)); SR
        StableReductionResult(ok, type=0---0e)

    This example shows that the problem described in commit d3eca6a has been fixed:

        sage: F = -8*x^4 - 2*y^4 + 4*x^3*z + y^3*z + 5*x^2*z^2 - 4*y^2*z^2 + 5*z^4
        sage: SR = stable_reduction_of_quartic(F, QQ.valuation(2)); SR
        StableReductionResult(ok, type=1ee)

    """
    K = F.base_ring()
    if not v_K.domain() == K:
        raise ValueError("base field of F must be the domain of v_K")
    res = StableReductionResult(status="fail", F=F, v_K=v_K)
    try:
        # 1.) Find GIT-semistable model
        X = PlaneCurveOverValuedField(F, v_K)
        XX = X.git_semistable_model_with_rational_cusps()
        v_L = XX.base_ring_valuation()
        L = v_L.domain()
        # define the special fiber, with reduction and lifting maps
        # this should later be moved to `plane_curve_valued`
        k = v_L.residue_field()
        k1 = GF(k.cardinality())
        phi = k.an_embedding(k1)
        phi_inv = phi.inverse()
        red0 = v_L.reduce
        red = lambda a: phi(red0(a))
        lift0 = v_L.lift
        lift = lambda a: lift0(phi_inv(a))
        f = XX.defining_polynomial().map_coefficients(red, k1)
        Xs = ProjectivePlaneCurve(f)
        # Xs = XX.special_fiber()   # old code
        res.git_model = XX
        res.git_special_fiber = Xs
        res.git_extension = L
        res.git_extension_valuation = v_L
        
        # 2) detect hyperelliptic reduction early
        if not Xs.is_git_stable():
            res.status = "hyperelliptic"
            return res

        # 3) resolve cusps -> tail types
        cusps = Xs.rational_cusps()
        # this is a list of `flags`
        cusp_data = []  # list of pairs (P, "e"/"m")
        # must be a list of pairs (P, tail_type), where P is a *rational*
        # point and ``tail_type`` is "e" or "m"     
        for C in cusps:
            T = C.move_to_e0_x2()
            # we want to lift T to a matrix in L; for this we first have
            # to change the base ring of T to the residue field of v_L
            M = T.map_coefficients(lift, L)
            cusp_model = XX.apply_matrix(M)
            tail = resolve_cusp(cusp_model.defining_polynomial(), v_L, 
                                   compute_matrix=compute_matrix)
            E = tail[2]
            P = Xs.point(C.point)
            res.tail_data[P] = tail
            if E.is_smooth():
                cusp_data.append((P, "e"))
            else:
                cusp_data.append((P, "m"))
        
        # 4) build component graph + classify
        G = component_graph_of_GIT_stable_quartic(Xs, cusp_data=cusp_data)
        res.component_graph = G
        res.canonical_signature = G.canonical_signature()
        res.reduction_type = classify_genus3_type(G)
        res.status = "ok"
        return res

    except Exception as e:
        res.warnings.append(f"Exception: {e}")
        return res


class StableReductionResult:
    """
    Lightweight container for the stable reduction computation of a plane quartic.

    The object is meant to be:
      - convenient to inspect interactively,
      - easy to serialize (at least in 'light' form) for databases,
      - extensible as the pipeline grows.
    """

    def __init__(self, status, F, v_K):
        # status: "ok", "hyperelliptic", "fail"
        self.status = status

        # input provenance
        self.K = v_K.domain()               
        self.v_K = v_K
        self.F = F

        # GIT model data 
        self.git_model = None                  # a GIT-semistable model of X 
        self.git_special_fiber = None          # its special fiber
        self.git_extension = None              # the field extension over which a git ss model exists
        self.git_extension_valuation = None    # its p-adic valuation
        
        # tail data
        self.tail_data = {}              # keys: points P, values:list of triples (v_L, E, T) or (v_L, E, e)

        # final combinatorics
        self.component_graph = None      # ComponentGraph
        self.canonical_signature = None  # the canonical signature of the graph
        self.reduction_type = None       # string among 42 types

        # logging / diagnostics
        self.warnings = []
        self.debug = {}                  # anything else you want to stash

    # ---------------------------
    # Convenience helpers
    # ---------------------------

    def is_ok(self):
        return self.status == "ok"

    def summary(self):
        if self.status == "ok":
            return (f"StableReductionResult(ok, type={self.reduction_type})")
        return f"StableReductionResult({self.status})"

    def __repr__(self):
        return self.summary()
    
    def to_record(self):
        """
        Return a lightweight record (dict) for the example database.

        Only defined for status == "ok". Otherwise returns None.

        The returned dictionary contains the following entries:

        - ``"K_defpoly"``: the defining polynomial of the number field `K`
        (as a string). If `K = QQ`, this is `"QQ"`.

        - ``"K_gen"``: the name of the chosen generator of `K`
        (if applicable).

        - ``"p"``: the rational prime defining the p-adic valuation.

        - ``"F"``: the homogeneous quartic form defining the curve,
        given as a string in `K[x,y,z]`.

        - ``"reduction_type"``: the reduction type of the stable
        reduction (one of the 42 genus-3 types).

        - ``"git_dfe"``: a triple ``(d, f, e)`` describing the extension
        `L_0/K` over which a GIT-semistable model exists, where
            * `d = [L_0 : K]` is the field degree,
            * `f` is the residue field degree, and
            * `e` is the ramification index (so `d = e * f`).

        - ``"cusp_dfes"``: a list of triples ``(d, f, e)``, one for each
        cusp, describing the field extension `L/L_0` required to resolve
        that cusp.  Again, `d = e * f`.

        These data are sufficient to reconstruct the example and to
        filter examples by reduction type and ramification behaviour.

        """
        if self.status != "ok":
            return None

        K = self.K
        if K is QQ:
            K_defpoly = "QQ"
            K_gen = None
        else:
            K_defpoly = str(K.defining_polynomial())
            K_gen = str(K.gen())
        p = int(self.v_K.p())  # assuming this works for your valuation
        F_str = str(self.F)

        # (d,f,e) for L0/K
        v0 = self.git_extension_valuation
        git_dfe = extension_triple(self.v_K, v0)

        # cusp extension data:
        # tail_data entries are either (v_L, t, E, T) or (v_L, t, E, e)
        cusp_dfes = []
        for tail in self.tail_data.values():
            vL = tail[0]
            third = tail[3]

            # If resolve_cusp returned a matrix, then the cusp resolution was done over vL (e=1)
            if hasattr(third, "nrows"):  # crude: matrix has nrows()
                # extension is whatever vL is over v0
                dfe = extension_triple(v0, vL)
                cusp_dfes.append(dfe)
                continue

            # Otherwise third is "e_needed" (positive int)
            e_needed = int(third)

            # Case 1: vL already extends v0 (possibly nontrivially)
            if vL.domain() != v0.domain():
                dfe0 = extension_triple(v0, vL)
            else:
                dfe0 = (1, 1, 1)

            # Case 2: additional totally ramified extension of index e_needed over vL
            # (this is exactly how resolve_cusp constructs it when compute_matrix=True)
            dfe1 = (e_needed, 1, e_needed) if e_needed != 1 else (1, 1, 1)

            # Store the combined extension over L0 as (d,f,e) if you want,
            # or store both levels. You asked for L/L0, so combine:
            d = dfe0[0] * dfe1[0]
            f = dfe0[1] * dfe1[1]
            e = dfe0[2] * dfe1[2]
            cusp_dfes.append((d, f, e))
        cusp_dfes.sort()  # to have a consistent output
        payload = f"{K_defpoly}|{p}|{F_str}".encode("utf-8")
        id = hashlib.sha256(payload).hexdigest()[:16]


        rec = {
            "id": id,
            "K_defpoly": K_defpoly,
            "K_gen": K_gen,
            "p": p,
            "F": F_str,
            "reduction_type": self.reduction_type,
            "git_dfe": git_dfe,              # (d,f,e) for L0/K
            "cusp_dfes": cusp_dfes,          # list of (d,f,e) for each cusp extension L/L0
        }
        return rec

    def to_json_record(self):
        r""" Return the record in json format.
        
        For the moment, this just means that tuples are turned into lists.
        """
        rec = self.to_record()
        if rec is None:
            return None
        rec = dict(rec)  # shallow copy
        rec["git_dfe"] = list(rec["git_dfe"])
        rec["cusp_dfes"] = [list(t) for t in rec["cusp_dfes"]]
        return rec


# ---------------------------------------------------------------------------

# Helper functions

def _residue_degree(v):
    k = v.residue_field()
    # the following gives wrong result if k is a quotient ring
    # Finite field: degree over prime field
    # if hasattr(k, "degree"):
    #     try:
    #         return int(k.degree())
    #     except Exception:
    #         pass
    # Fallback for finite fields where degree() may not exist
    if hasattr(k, "cardinality") and hasattr(k, "characteristic"):
        q = int(k.cardinality())
        p = int(k.characteristic())
        # q = p^f
        f = 0
        t = 1
        while t < q:
            t *= p
            f += 1
        if t == q:
            return f
    raise ValueError("Cannot determine residue field degree from valuation.")

def extension_triple(v_base, v_ext):
    """
    Return (d,f,e) for the extension of valued fields induced by v_ext over v_base.
    Assumes v_ext extends v_base.
    """
    K = v_base.domain()
    L = v_ext.domain()
    # field degree d = [L:K]
    d = int(L.absolute_degree() / K.absolute_degree()) 
    # residue degree f
    fK = _residue_degree(v_base)
    fL = _residue_degree(v_ext)
    if fL % fK != 0:
        raise ValueError(f"Residue degree not divisible; inconsistent valuations? f_K = {fK}, f_L = {fL}")
    f = int(fL // fK)
    if d % f != 0:
        raise ValueError("d not divisible by f; inconsistent (d,f) data?")
    e = int(d // f)
    return (d, f, e)
