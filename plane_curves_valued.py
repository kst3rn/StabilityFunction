
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from plane_curves import ProjectivePlaneCurve
from extension_search import find_semistable_model


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


  def degree(self):
    return self.defining_polynomial().degree()


  def semistable_model(self):
    r"""
    Return a semistable model of `self`.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: F = 16*x^4 + y^4 + 8*y^3*z + 16*x*y*z^2 + 4*x*z^3
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.semistable_model(); X
      Plane Model of Projective Plane Curve with defining polynomial 16*x^4 + y^4 + 8*y^3*z + 16*x*y*z^2 + 4*x*z^3 over Rational Field with 2-adic valuation
      sage: X.base_ring()
      Number Field in piL with defining polynomial x^12 + 2*x^6 + 2
      sage:
      sage: F = 4*x^4 + 4*x*y^3 + y^4 + 2*x*z^3 + 4*y*z^3 + z^4
      sage: Y = PlaneCurveOverValuedField(F, QQ.valuation(2))
      sage: X = Y.semistable_model(); X
      Plane Model of Projective Plane Curve with defining polynomial 4*x^4 + 4*x*y^3 + y^4 + 2*x*z^3 + 4*y*z^3 + z^4 over Rational Field with 2-adic valuation
      sage: X.base_ring()
      Number Field in piK with defining polynomial x^4 + 2*x^3 + 2*x^2 + 2
    """
    btb_point, F = find_semistable_model(self.defining_polynomial(), self.base_ring_valuation())
    return PlaneModel(F, self, btb_point)



class PlaneModel(ProjectivePlaneCurve):
  r"""
  Construct...
  """

  def __init__(self, polynomial, generic_fiber, btb_point):
    super().__init__(polynomial)
    self._generic_fiber = generic_fiber
    self._bruhat_tits_building_point = btb_point


  def __repr__(self):
    return f"Plane Model of {self.generic_fiber()}"


  def generic_fiber(self):
    return self._generic_fiber


  def point_on_BruhatTitsBuilding(self):
    return self._bruhat_tits_building_point


  def special_fiber(self):
    r"""
    Return the special fiber of `self`.
    """
    v = self.point_on_BruhatTitsBuilding().linear_valuation()
    return ProjectivePlaneCurve(v.reduction(self.defining_polynomial()))


  def has_semistable_reduction(self):
    r"""
    Return `True` if the special fiber of `self` is
    semistable and `False` otherwise.
    """
    return self.special_fiber().is_semistable()

