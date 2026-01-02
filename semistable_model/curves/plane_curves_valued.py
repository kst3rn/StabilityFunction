
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import identity_matrix
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


  def semistable_model(self):
    r"""
    Return a semistable model of `self`.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = 16*x^4 + y^4 + 8*y^3*z + 16*x*y*z^2 + 4*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.semistable_model()
      sage: X.base_ring()
      Number Field in piL with defining polynomial x^12 + 2*x^6 + 2
      sage: X.has_semistable_reduction()
      True
      sage: X.special_fiber()
      Projective Plane Curve with defining polynomial x^4 + x^2*y^2 + y*z^3 over Finite Field of size 2
      sage: X.generic_fiber().base_ring()
      Number Field in piL with defining polynomial x^12 + 2*x^6 + 2
      sage:
      sage: F = 4*x^4 + 4*x*y^3 + y^4 + 2*x*z^3 + 4*y*z^3 + z^4
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.semistable_model()
      sage: X.base_ring()
      Number Field in piK with defining polynomial x^4 + 2*x^3 + 2*x^2 + 2
      sage: X.has_semistable_reduction()
      True
      sage: X.special_fiber()
      Projective Plane Curve with defining polynomial x^4 + x^2*y^2 + x*y^3 + y^3*z + y^2*z^2 + y*z^3 over Finite Field of size 2
      sage: X.generic_fiber().base_ring()
      Number Field in piK with defining polynomial x^4 + 2*x^3 + 2*x^2 + 2
    """
    from semistable_model.stability import semistable_reduction_field
    from semistable_model.stability import StabilityFunction

    F = self.defining_polynomial()
    v_K = self.base_ring_valuation()
    L = semistable_reduction_field(F, v_K)
    v_L = v_K.extension(L)
    X_L = self.base_change(v_L)
    phiL = StabilityFunction(X_L.defining_polynomial(), v_L)
    a, b = phiL.global_minimum()
    T = b.move_to_origin().base_change_matrix()
    return PlaneModel(X_L, T)


    def semistable_models_with_canonical_cusps(self):
      r"""
      Return a list with semistable models of `self` with cusps in canonical form.
      """
      raise NotImplementedError



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
    return ProjectivePlaneCurve(v.reduction(self.defining_polynomial()))


  def has_semistable_reduction(self):
    r"""
    Return `True` if the special fiber of `self` is
    semistable and `False` otherwise.
    """
    return self.special_fiber().is_semistable()

