
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


  def semistable_model(self):
    raise NotImplementedError("work in progress")


class PlaneModel(ProjectivePlaneCurve):
  r"""
  Construct...
  """

  def __init__(self, polynomial, generic_fiber, base_change_matrix):
    super().__init__(polynomial)
    self._generic_fiber = generic_fiber
    self._base_change_matrix = base_change_matrix

