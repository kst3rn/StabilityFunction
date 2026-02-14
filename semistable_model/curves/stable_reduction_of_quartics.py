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
"""

from sage.all import Curve
from semistable_model.curves.plane_curves_valued import PlaneCurveOverValuedField
from semistable_model.curves.cusp_resolution import resolve_cusp
from semistable_model.curves.component_graphs_of_plane_curves import component_graph_of_GIT_stable_quartic
from semistable_model.curves.genus3_reduction_types import classify_genus3_type


def stable_reduction_of_quartic(F, v_K):
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

    The following example produces an error:

        sage: R.<x,y,z> = QQ[]
        sage: F = -x^3*y - 8*x*y^3 - 7*x^3*z - 7*x^2*y*z + 5*x*y^2*z + 6*x^2*z^2 - 6*y*z^3
        sage: SR = stable_reduction_of_quartic(F, QQ.valuation(2)); SR  # long time
        StableReductionResult(fail)

        sage: SR.warnings
        ["Exception: unsupported operand parent(s) for +: ... over Finite Field in z2 of size 2^2'"]

    The error is caused by a bug in sage, see issue https://github.com/sagemath/sage/issues/41643

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
        Xs = XX.special_fiber()
        res.git_model = XX
        res.special_fiber = Xs
        res.v_L = v_L
        # the attribute "git_wild" registers whether the the git-semistable model
        # is defined over a wild or tame extension of K
        if v_K.p().divides(v_L(v_K.uniformizer())/v_L(v_L.uniformizer())):
            res.git_wild = True

        else:
            res.git_wild = False

        # 2) detect hyperelliptic reduction early
        if not Xs.is_git_stable():
            res.status = "hyperelliptic"
            return res

        # 3) resolve cusps -> tail types
        cusps = Xs.rational_cusps()
        # this is a list of `flags`
        # below we need Xs to be an object of `Curve`
        Xs = Curve(Xs.defining_polynomial()) 
        cusp_data = [] 
        # must be a list of pairs (P, tail_type), where P is a *rational*
        # point and ``tail_type`` is "e" or "m"     
        for C in cusps:
            T = C.move_to_e0_x2()
            # we want to lift T to a matrix in L; for this we first have
            # to change the base ring of T to the residue field of v_L
            phi = T.base_ring().an_embedding(v_L.residue_field())
            M = T.map_coefficients(phi, v_L.residue_field()).map_coefficients(v_L.lift, L)
            cusp_model = XX.apply_matrix(M)
            _, E, _ = resolve_cusp(cusp_model.defining_polynomial(), v_L)
            P = Xs.point(C.point)
            if E.is_smooth():
                cusp_data.append((P, "e"))
            else:
                cusp_data.append((P, "m"))

        # res.tail_types = {P:"e"/"m", ...}

        # 4) build component graph + classify
        G = component_graph_of_GIT_stable_quartic(Xs, cusp_data=cusp_data)

        # the above gives an error in rare examples:
        # F = -x^3*y - 8*x*y^3 - 7*x^3*z - 7*x^2*y*z + 5*x*y^2*z + 6*x^2*z^2 - 6*y*z^3
        # v_K = QQ.valuation(2)

        res.component_graph = G
        res.signature = G.canonical_signature()
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
        self.F_input = F

        # normalization / transformations
        self.F_normalized = None
        self.change_of_coordinates = None   # e.g. matrix in GL3(K)

        # model data (GIT model etc.)
        self.git_model = None               # your model object, if you keep it
        self.special_fiber = None           # plane curve over residue field (or splitting field)
        self.residue_field = None
        self.splitting_field = None
        self.field_extension_data = {}      # e.g. degrees, defining polynomials

        # cusp / tail data
        self.cusps = []                     # list of points on special_fiber
        self.tail_types = {}                # cusp -> "e"/"m"
        self.cusp_resolution_data = {}      # optional: per cusp diagnostics

        # final combinatorics
        self.component_graph = None         # ComponentGraph
        self.signature = None               # canonical signature (tuple)
        self.reduction_type = None          # string among 42 types

        # logging / diagnostics
        self.warnings = []
        self.debug = {}                     # anything else you want to stash

    # ---------------------------
    # Convenience helpers
    # ---------------------------

    def is_ok(self):
        return self.status == "ok"

    def summary(self):
        if self.status == "ok":
            return (f"StableReductionResult(ok, type={self.reduction_type}, "
                    f"sig={self.signature})")
        return f"StableReductionResult({self.status})"

    def __repr__(self):
        return self.summary()

    # ---------------------------
    # Serialization
    # ---------------------------

    def to_dict(self, light=True):
        """
        Convert to a dict suitable for storage.

        If light=True, omit heavy Sage objects and keep only strings/primitive data.
        """
        d = {
            "status": self.status,
            "reduction_type": self.reduction_type,
            "signature": self.signature,
            "warnings": list(self.warnings),
            "field_extension_data": dict(self.field_extension_data),
        }

        # Store polynomials as strings if possible
        if self.F_input is not None:
            d["F_input"] = str(self.F_input)
        if self.F_normalized is not None:
            d["F_normalized"] = str(self.F_normalized)

        # Cusps / tails: store as strings (points can be non-JSON)
        if self.cusps:
            d["cusps"] = [str(P) for P in self.cusps]
        if self.tail_types:
            d["tail_types"] = {str(P): t for P, t in self.tail_types.items()}

        if not light:
            # keep heavy objects for pickling / in-session storage
            d.update({
                "K": self.K,
                "v_K": self.v_K,
                "git_model": self.git_model,
                "special_fiber": self.special_fiber,
                "residue_field": self.residue_field,
                "splitting_field": self.splitting_field,
                "component_graph": self.component_graph,
                "debug": dict(self.debug),
            })

        return d
