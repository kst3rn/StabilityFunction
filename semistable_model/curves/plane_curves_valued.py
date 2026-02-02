
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import matrix, identity_matrix, PolynomialRing, GF
from semistable_model.curves import ProjectivePlaneCurve
from semistable_model.valuations import LinearValuation


class PlaneCurveOverValuedField(ProjectivePlaneCurve):
  r"""
  Construct...

  EXAMPLES::
    sage: R.<x,y,z> = QQ[]
    sage: F = z*y^2 - x^3 - x*z^2
    sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2)); Y
    Projective Plane Curve with defining polynomial -x^3 + y^2*z - x*z^2 over Rational Field with 2-adic valuation
  """

  def __init__(self, polynomial, base_ring_valuation):
    r"""
    Construct...

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = z*y^2 - x^3 - x*z^2
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      Projective Plane Curve with defining polynomial -x^3 + y^2*z - x*z^2 over Rational Field with 2-adic valuation
    """
    super().__init__(polynomial)
    self._base_ring_valuation = base_ring_valuation


  def __repr__(self):
    return f"Projective Plane Curve with defining polynomial {self.defining_polynomial()} over {self.base_ring()} with {self.base_ring_valuation()}"


  def base_ring_valuation(self):
    return self._base_ring_valuation


  def base_ring(self):
    return self.base_ring_valuation().domain()


  def base_change(self, valuation_extension):
    r"""
    Return the base change of `self`.
    """
    PolRin0 = self.defining_polynomial().parent()
    PolRin1 = PolRin0.change_ring(valuation_extension.domain())
    phi = PolRin1.coerce_map_from(PolRin0)
    if phi is None:
      raise NotImplementedError(f"No coercion from the polynomial ring over {self.base_ring()} to the polynomial ring over {R}")
    new_poly = phi(self.defining_polynomial())
    return PlaneCurveOverValuedField(new_poly, valuation_extension)


  def degree(self):
    return self.defining_polynomial().degree()


  def git_semistable_model(self, ramification_index=None):
    r"""
    Return a semistable model of `self`.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = y^4 + 2*x^3*z + x*y^2*z + 2*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.git_semistable_model()
      sage: X.base_ring()
      Number Field in piK with defining polynomial x^2 + 2
      sage:
      sage: F = 16*x^4 + y^4 + 8*y^3*z + 16*x*y*z^2 + 4*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.git_semistable_model()
      sage: X.base_ring()
      Number Field in piL with defining polynomial x^12 + 2*x^6 + 2
      sage: X.generic_fiber().base_ring() == X.base_ring()
      True
      sage: X.has_semistable_reduction()
      True
      sage: X.special_fiber()
      Projective Plane Curve with defining polynomial x^4 + x^2*y^2 + y*z^3 over Finite Field of size 2
      sage:
      sage: F = 4*x^4 + 4*x*y^3 + y^4 + 2*x*z^3 + 4*y*z^3 + z^4
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.git_semistable_model()
      sage: X.base_ring()
      Number Field in piK with defining polynomial x^4 + 8*x + 4
      sage: X.generic_fiber().base_ring() == X.base_ring()
      True
      sage: X.has_semistable_reduction()
      True
    """
    from semistable_model.stability import semistable_reduction_field
    from semistable_model.stability import StabilityFunction

    F = self.defining_polynomial()
    v_K = self.base_ring_valuation()
    L = semistable_reduction_field(F, v_K, ramification_index)
    v_L = v_K.extension(L)
    X_L = self.base_change(v_L)
    phiL = StabilityFunction(X_L.defining_polynomial(), v_L)
    a, b = phiL.global_minimum()
    T = b.move_to_origin().base_change_matrix()
    return PlaneModel(X_L, T)


  def git_semistable_model_with_rational_cusps(self, ramification_index=None):
    r"""
    Return a semistable models of `self` such that all cusps of its
    reduction are rational.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = y^4 + 2*x^3*z + x*y^2*z + 2*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.git_semistable_model_with_rational_cusps()
      sage: X.base_ring()
      Number Field in a1 with defining polynomial x^4 - 2*x^3 + x^2 - 6*x + 9
      sage: Xs = X.special_fiber()
      sage: Xs.rational_cusps()
      [Projective flag given by [u1, 1, 1] and u1*x + u1*y + z,
      Projective flag given by [u1 + 1, u1, 1] and u1*x + z]
    """
    X = self.git_semistable_model(ramification_index)
    Xs = X.special_fiber()
    L_tr = X.base_ring()
    d = Xs.splitting_field_of_singular_points().degree()
    Rk = PolynomialRing(GF(2), 'x')
    g_bar = Rk.irreducible_element(d)
    g = g_bar.change_ring(L_tr)
    L_mixed_relative = L_tr.extension(g, names='b')
    L_mixed_absolute = L_mixed_relative.absolute_field(names='a')
    L_mixed_absolute = L_mixed_absolute.optimized_representation()[0]
    v_K = self.base_ring_valuation()
    v_L_mixed = v_K.extension(L_mixed_absolute)
    from semistable_model.stability import StabilityFunction
    Y_L_mixed = self.base_change(v_L_mixed)
    phi = StabilityFunction(Y_L_mixed.defining_polynomial(), v_L_mixed)
    a, b = phi.global_minimum()
    T = b.move_to_origin().base_change_matrix()
    return PlaneModel(Y_L_mixed, T)


  def git_semistable_models_with_e2_x0_cusps(self, ramification_index=None):
    r"""
    Return a list of semistable models such that all cusps of their
    reductions are rational and at least one cusp is in canonical form.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = y^4 + 2*x^3*z + x*y^2*z + 2*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X1, X2 = Y.git_semistable_models_with_e2_x0_cusps()
      sage: X1.special_fiber().rational_cusps()
      [Projective flag given by [0, 0, 1] and x,
      Projective flag given by [u1, u1, 1] and x + u1*y + z]
      sage: X2.special_fiber().rational_cusps()
      [Projective flag given by [0, 0, 1] and x,
      Projective flag given by [1, u1, 1] and (u1 + 1)*x + y + z]
    """
    X = self.git_semistable_model_with_rational_cusps(ramification_index)
    L = X.base_ring()
    Xs = X.special_fiber()
    cusps = Xs.rational_cusps()
    v = X.base_ring_valuation()
    models = []
    for C in cusps:
      T = C.move_to_e2_x0()
      M = [[0,0,0],[0,0,0],[0,0,0]]
      for i in range(3):
        for j in range(3):
          M[i][j] = v.lift(T[i][j])
      M = matrix(L, M)
      models.append(X.apply_matrix(M))
    return models


  def git_semistable_models_with_e0_x2_cusps(self, ramification_index=None):
    r"""
    Return a list of semistable models such that all cusps of their
    reductions are rational and at least one cusp is in canonical form.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = y^4 + 2*x^3*z + x*y^2*z + 2*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X1, X2 = Y.git_semistable_models_with_e0_x2_cusps()
      sage: X1.special_fiber().rational_cusps()
      [Projective flag given by [1, 0, 0] and z,
      Projective flag given by [u1 + 1, 1, 1] and x + u1*y + z]
      sage: X2.special_fiber().rational_cusps()
      [Projective flag given by [1, 0, 0] and z,
      Projective flag given by [1, u1, 1] and u1*x + u1*y + z]
    """
    X = self.git_semistable_model_with_rational_cusps(ramification_index)
    L = X.base_ring()
    Xs = X.special_fiber()
    cusps = Xs.rational_cusps()
    v = X.base_ring_valuation()
    models = []
    for C in cusps:
      T = C.move_to_e0_x2()
      M = [[0,0,0],[0,0,0],[0,0,0]]
      for i in range(3):
        for j in range(3):
          M[i][j] = v.lift(T[i][j])
      M = matrix(L, M)
      models.append(X.apply_matrix(M))
    return models



class PlaneModel(ProjectivePlaneCurve):
  r"""
  Construct a plane model of a projective plane curve
  over a valued field to the following conditions.

  INPUT:
  - ``generic_fiber`` -- a curve over a valued field.
  - ``base_change_matrix`` -- an invertible matrix.
  """

  def __init__(self, generic_fiber, base_change_matrix):
    r"""
    Construct a plane model of a projective plane curve
    over a valued field.
    """
    if not base_change_matrix.is_invertible():
      raise ValueError("The base change matrix must be invertible.")
    from semistable_model.stability import BTB_Point
    b = BTB_Point(generic_fiber.base_ring_valuation(),
                  base_change_matrix,
                  [0,0,0])
    F = b.hypersurface_model(generic_fiber.defining_polynomial())
    self._bruhat_tits_building_point = b
    super().__init__(F)
    self._generic_fiber = generic_fiber


  def __repr__(self):
    return f"Plane Model of {self.generic_fiber()}"


  def generic_fiber(self):
    return self._generic_fiber


  def base_ring_valuation(self):
    return self.generic_fiber().base_ring_valuation()


  def as_point_on_BTB(self):
    return self._bruhat_tits_building_point


  def base_change_matrix(self):
    return self.as_point_on_BTB().base_change_matrix()


  def adapted_basis(self):
    return self.as_point_on_BTB().linear_valuation().adapted_basis()


  def apply_matrix(self, T):
    r"""
    Return the plane model of the generic fiber of `self` with
    base change matrix given by `T * self.base_change_matrix()`.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = 16*x^4 + y^4 + 8*y^3*z + 16*x*y*z^2 + 4*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = PlaneModel(Y, identity_matrix(QQ, 3))
      sage: X.base_change_matrix()
      [1 0 0]
      [0 1 0]
      [0 0 1]
      sage: T = matrix(QQ, [[1,0,0],[1,1,0],[0,0,1]])
      sage: T_X = X.apply_matrix(T)
      sage: T_X.base_change_matrix()
      [1 0 0]
      [1 1 0]
      [0 0 1]
      sage: M = matrix(QQ, [[1,0,0],[0,1,0],[0,1,1]])
      sage: MT_X = T_X.apply_matrix(M)
      sage: MT_X.base_change_matrix()
      [1 0 0]
      [1 1 0]
      [1 1 1]
    This method defines a left action.
      sage: MT_X.base_change_matrix() == M*T
      True
      sage: MT_X.base_change_matrix() == T*M
      False
    """
    return PlaneModel(self.generic_fiber(), T * self.base_change_matrix())


  def special_fiber(self):
    r"""
    Return the special fiber of `self`.
    """
    F = self.defining_polynomial()
    R = F.parent()
    E = identity_matrix(R.base_ring(), R.ngens())
    v_K = self.as_point_on_BTB().base_ring_valuation()
    v = LinearValuation(R, v_K, E, [0]*R.ngens())
    f_wrong = v.reduction(self.defining_polynomial())
    k_wrong = f_wrong.base_ring()
    p = k_wrong.characteristic()
    d = k_wrong.degree()
    k_right = GF(p**d)
    phi = k_wrong.an_embedding(k_right)
    f_right = f_wrong.change_ring(phi)
    return ProjectivePlaneCurve(f_right)


  def has_semistable_reduction(self):
    r"""
    Return `True` if the special fiber of `self` is
    semistable and `False` otherwise.
    """
    return self.special_fiber().is_semistable()


  def resolve_cusp(self):
    r"""
    If `self.special_fiber()` has no rational cusp given by
    (1:0:0) and V_+(x_2) an error is raised.

    EXAMPLES::
    First, we compute two models with rational cusps in the right position.
      sage: R.<x,y,z> = QQ[]
      sage: F = y^4 + 2*x^3*z + x*y^2*z + 2*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X1, X2 = Y.git_semistable_models_with_e0_x2_cusps()
      sage: X1.base_ring()
      Number Field in a1 with defining polynomial x^4 - 2*x^3 + x^2 - 6*x + 9

    Now we can resolve the cusps of X1 and X2.
      sage: v_L1, T1, F1b = X1.resolve_cusp()
      sage: v_L1
      2-adic valuation
      sage: v_L1.domain()
      Number Field in alpha with defining polynomial a^8 + (8*a1^3 + 8*a1)*a^7 + (8*a1^3 + 8*a1^2 + 8)*a^6 + (16*a1 + 16)*a^5 + (16*a1^3 - 8*a1^2 + 16*a1 - 12)*a^4 + (48*a1^3 + 16*a1)*a^3 + (16*a1^3 + 32*a1^2 + 40*a1 - 4)*a^2 + (8*a1^3 + 16*a1^2 + 8*a1 + 64)*a + 24*a1^3 + 32*a1^2 + 16*a1 - 4 over its base field
      sage: F1b
      y^3 + x^2*z + x*z^2
      sage: F1b.base_ring()
      Finite Field in u1 of size 2^2
      sage:
      sage: v_L2, T2, F2b = X2.resolve_cusp()
      sage: v_L2.domain()
      Number Field in alpha with defining polynomial a^8 + 8*a^7 + (8*a1 + 8)*a^6 + 16*a1*a^5 + (-8*a1^3 + 4*a1^2 - 8*a1 + 8)*a^4 + (32*a1^3 + 16*a1^2 + 32*a1)*a^3 + (16*a1^3 + 40*a1 + 12)*a^2 + (32*a1^3 + 48*a1 - 8)*a + 24*a1^3 - 12*a1^2 + 64*a1 + 52 over its base field
      sage: F2b
      (u1 + 1)*y^3 + u1*x^2*z + x*z^2
      sage: F2b.base_ring()
      Finite Field in u1 of size 2^2
      sage: F1b.base_ring()
      Finite Field in u1 of size 2^2
    """
    from semistable_model.curves import resolve_cusp
    return resolve_cusp(self.defining_polynomial(), self.base_ring_valuation())

