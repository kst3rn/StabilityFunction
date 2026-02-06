# Test that GIT-stable models over ZZ_p, for p > 3, 
# are defined over a tamely ramified extension

from sage.all import QQ, ZZ, PolynomialRing, Curve
from semistable_model.curves.stable_reduction_of_quartics import stable_reduction_of_quartic

def random_smooth_quartic(coeff_bound=1000):
    R.<x,y,z> = QQ[]
    while True:
        F = R.zero()
        for m in R.monomials_of_degree(4):
            if random() > 0.7:
                F += ZZ.random_element(coeff_bound)*m
        if not F.is_zero() and Curve(F).is_smooth():
            return F

def test_tame_assumption(p, N=1000):
    for i in range(N):
        F = random_smooth_quartic()
        X = PlaneCurveOverValuedField(F, QQ.valuation(p))
        XX = X.semistable_model()
        v_L = XX.base_ring_valuation()
        if XX.special_fiber().is_stable() and p.divides(1/v_L(v_L.uniformizer())):
            print(f"X = {X} has semistable model over {XX.base_ring_valuation()}")
