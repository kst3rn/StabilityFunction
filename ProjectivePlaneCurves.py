
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from functools import cached_property
from functools import cache
from sage.all import *


class ProjectivePlaneCurve:
  r"""
  Construct a projective plane curve to the following conditions.

  INPUT:
  - ``polynomial`` -- homogeneous polynomial in K[x_0, x_1, x_2].
  """

  def __init__(self, polynomial):
    r"""
    Construct a projective plane curve to the following conditions.

    INPUT:
    - ``polynomial`` -- homogeneous polynomial in K[x_0, x_1, x_2].
    """

    if not polynomial.is_homogeneous():
      raise TypeError
    if not len(polynomial.parent().gens()) == 3:
      raise ValueError

    if not polynomial.base_ring().is_field():
      raise ValueError

    self.polynomial = polynomial
    self._degree = self.polynomial.degree()
    self.polynomial_ring = self.polynomial.parent()
    self.base_ring = self.polynomial.base_ring()
    self.projective_plane = ProjectiveSpace(self.polynomial_ring)
    self.plane_curve = self.projective_plane.curve(self.polynomial)
    self.standard_basis = self.polynomial_ring.gens()


  def __repr__(self):
    return f"Projective Plane Curve with defining polynomial {self.polynomial}"


  def get_base_ring(self):
    return self.base_ring


  def get_polynomial(self):
    return self.polynomial


  def get_standard_basis(self):
    return self.standard_basis


  def degree(self):
    return self._degree


  def tangent_cone_at(self, P):
    r"""
    Return the tangent cone of self at the point `P`.

    INPUT:
    - ``P`` -- a point on `self`.

    EXAMPLES:
      sage: R.<x0,x1,x2> = GF(3)[]
      sage: f = x0^2*x2 - x1^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x1^3 + x0^2*x2
      sage: P = [0,0,1]
      sage: C = X.tangent_cone_at(P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x1^3 + x0^2*x2 at [0, 0, 1]
      sage: C.get_polynomial()
      x^2

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0-3*x2)^2*x2 - (x1+5*x2)^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x1^3 + x0^2*x2 - 15*x1^2*x2 - 6*x0*x2^2 - 75*x1*x2^2 - 116*x2^3
      sage: P = [3,-5,1]
      sage: C = X.tangent_cone_at(P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x1^3 + x0^2*x2 - 15*x1^2*x2 - 6*x0*x2^2 - 75*x1*x2^2 - 116*x2^3 at [3, -5, 1]
      sage: C.get_polynomial()
      x^2

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x1^2*x2 - x0^3 - x0^2*x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x0^3 - x0^2*x2 + x1^2*x2
      sage: P = [0,0,1]
      sage: C = X.tangent_cone_at(P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x0^3 - x0^2*x2 + x1^2*x2 at [0, 0, 1]
      sage: C.get_polynomial()
      -x^2 + y^2

      sage: R.<x0,x1,x2> = GF(7)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: P = [1,2,4]
      sage: C = X.tangent_cone_at(P); C
      Tangent cone of Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4 at [2, 4, 1]
      sage: C.get_polynomial()
      -3*x - 3*y
    """

    return PPC_TangentCone(self, P)


  @cache
  def is_smooth(self):
    r"""
    Return `True` if `self` is smooth and `False` otherwise.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x1^2*x2 - x0^3 - x0^2*x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x0^3 - x0^2*x2 + x1^2*x2
      sage: X.is_smooth()
      False

      sage: R.<x0,x1,x2> = GF(7)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_smooth()
      True
    """

    return self.plane_curve.is_smooth()


  def is_reduced(self):
    r"""
    Return `True` if `self` is reduced and `False` otherwise.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_reduced()
      True

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_reduced()
      False

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0 + x1)^2*(x1 + x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x1 + 2*x0*x1^2 + x1^3 + x0^2*x2 + 2*x0*x1*x2 + x1^2*x2
      sage: X.is_reduced()
      False
    """

    return not any(multiplicity > 1 for factor, multiplicity in self._decompose)


  def is_irreducible(self):
    r"""
    Return `True` if `self` is irreducible and `False` otherwise.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_irreducible()
      True

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_irreducible()
      False

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0 + x1)*(x1 + x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1 + x1^2 + x0*x2 + x1*x2
      sage: X.is_irreducible()
      False
    """

    if len(self._decompose) > 1:
      return False
    return self.is_reduced()


  def is_semistable(self):
    r"""
    Return `True` if `self` is semistable and `False` otherwise.

    EXAMPLES:
    A nodal cubic is semistable.
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x1^2*x2 + x0^3 + x0^2*x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x0^2*x2 + x1^2*x2
      sage: X.is_semistable()
      True

    A cuspidal cubic is unstable.
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x1^2*x2 + x0^3 + x0^2*x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x0^2*x2 + x1^2*x2
      sage: X.is_semistable()
      False
      sage:
      sage: R.<x0,x1,x2> = GF(3)[]
      sage: f = x0^3 + x1^2 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x1^2*x2
      sage: X.is_semistable()
      False

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_semistable()
      True
      sage:
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_semistable()
      False

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0^2 + x1*x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + 2*x0^2*x1*x2 + x1^2*x2^2
      sage: X.is_semistable()
      True
    """

    if self.is_smooth():
      return True
    elif self.instability() is not None:
      return False
    return True


  def is_stable(self):
    r"""
    Return `True` if `self` is stable and `False` otherwise.

    EXAMPLES:
    A smooth curve is stable.
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.is_stable()
      True

    A singular but stable curve.
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x0*x1^2*x2 + x0*x1*x2^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x0*x1^2*x2 + x0*x1*x2^2
      sage: X.is_stable()
      True

    A nonreduced but stable curve.
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0^4 + x1^3*x2 + x2^4)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^8 + 2*x0^4*x1^3*x2 + x1^6*x2^2 + 2*x0^4*x2^4 + 2*x1^3*x2^5 + x2^8
      sage: X.is_stable()
      True

    A properly semistable curve.
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0^2 + x1*x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + 2*x0^2*x1*x2 + x1^2*x2^2
      sage: X.is_stable()
      False
    """

    if self.is_smooth():
      return True

    if not self.is_semistable():
      return False

    # X_red is a conic.
    if self.degree() % 2 == 0:
      G, m = self._decompose[0]
      if m == self.degree() / 2 and G.degree() == 2:
        return False

    # Search for a line of multiplicity d/3.
    if self.degree() % 3 == 0:
      for Y, m in self._decompose:
        if Y.degree() == 1 and m == self.degree() / 3:
          return False

    # Search for point of multiplicity 2d/3 or a point
    # of multiplicity d/3 < m <= 2d/3 and a line in the
    # tangent cone of multiplicity >= m/2.
    X_red_sing = self._reduced_singular_points
    for P in X_red_sing:
      m = self.multiplicity(P)
      if m == 2 * self.degree() / 3:
        return False
      elif m > self.degree() / 3:
        for L, L_mult in PPC_TangentCone(self, P).embedded_lines():
          if L_mult >= m / 2:
            if ProjectiveFlag(P, L).is_semiunstable(self):
              return False

    return True


  def instability(self):
    r"""
    Return an instability of `self` or `None` if `self` is
    semistable.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x1^2*x2 - x0^3 - x0^2*x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x0^3 - x0^2*x2 + x1^2*x2
      sage: X.instability()
      sage:

      sage: R.<x0,x1,x2> = GF(3)[]
      sage: f = x0^3 + x1^2 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x1^2*x2
      sage: X.instability()
      Flag attached to Projective Plane Curve with defining polynomial x0^3 + x1^2*x2 given by [0, 0, 1] and x1

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.instability()
      Flag attached to Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4 given by x0 + x1 + x2
    """

    # Search for a line of multiplicity > d/3.
    for Y, m in self._decompose:
      if Y.degree() == 1 and m > self.degree() / 3:
        return ProjectiveFlag(None, Y)

    # Search for a point of multiplicity > 2d/3 or a point
    # of multiplicity d/3 < m <= 2d/3 and a line in the
    # tangent cone of multiplicity >= m/2.
    X_red_sing = self._reduced_singular_points
    for P in X_red_sing:
      m = self.multiplicity(P)
      if m > 2 * self.degree() / 3:
        return ProjectiveFlag(P, None)
      elif m > self.degree() / 3:
        for L, L_mult in PPC_TangentCone(self, P).embedded_lines():
          if L_mult > m / 2:
            P_on_L_flag = ProjectiveFlag(P, L)
            if P_on_L_flag.is_unstable(self):
              return P_on_L_flag

    return None


  def elementary_instability_direction(self, shape):
    r"""
    Return the element `lambda` of the base ring of `self` such
    that `self` has an instability diagonalized by the basis
    corresponding to the elementary matrix `T_{ij}(lambda)` where
    i, j = shape.

    INPUT:
    - ``shape`` -- a pair `(i, j)` of distinct integers in {0, 1, 2}.

    OUTPUT:
    - ``lambda`` -- an element of self.get_base_ring() such that there
                    is an instability diagonalized by the basis
                    (x_0, x_1, x_2) * T_{ij}(lambda),
                    where
                    i = shape[0],
                    j = shape[1],
                    (x_0, x_1, x_2) = self.get_standard_basis(),
                    and `T_{ij}(lambda)` is the elementary matrix with
                    `lambda` at the (i,j)-th position.

    EXAMPLES:
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x1^2*x2 + x0^3 + x0^2*x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x0^2*x2 + x1^2*x2
      sage: X.elementary_instability_direction((1,0))
      1
      sage: X.elementary_instability_direction((2,0))
      None
      sage: X.elementary_instability_direction((2,1))
      None

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x1^2*x2 - x0^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x0^3 + x1^2*x2
      sage: X.elementary_instability_direction((1,0))
      None
      sage: X.elementary_instability_direction((2,0))
      0
      sage: X.elementary_instability_direction((2,1))
      None

    REMARK:
    This method does not search for instabilities that are
    diagonalized by self.get_standard_basis().
    """

    if shape[0] == shape[1]:
      raise ValueError(f"The entries of {shape} must be distinct.")

    i, j = shape
    k = 3 - i - j
    x_i = self.get_standard_basis()[i]
    x_j = self.get_standard_basis()[j]
    x_k = self.get_standard_basis()[k]

    # Search for a line of multiplicity > d/3.
    for G, m in self._decompose:
      if G.degree() == 1 and m > self.degree() / 3:
        G_vars = set(G.variables())
        if G_vars == {x_j, x_i}:
          return G[x_i] / G[x_j]
        if m > 2 * self.degree() / 3 and G_vars.issupset({x_j, x_i}):
          return G[x_i] / G[x_j]

    # Search for a point of multiplicity > 2d/3 or a point
    # of multiplicity d/3 < m <= 2d/3 and a line in the
    # tangent cone of multiplicity >= m/2.
    X_red_sing = self._reduced_singular_points
    for P in X_red_sing:
      m = self.multiplicity(P)
      if m > 2 * self.degree() / 3 and P[j] != 0 and P[i] != 0 and P[k] == 0:
        return -P[j] / P[i]
      elif m > self.degree() / 3:
        for L, L_mult in PPC_TangentCone(self, P).embedded_lines():
          if L_mult > m / 2 and ProjectiveFlag(P, L).is_unstable(self):
            L_vars = set(L.variables())
            if L_vars == {x_j, x_i}:
              lambda_L = L[x_i] / L[x_j]
              if P[j] == 0 and P[i] == 0:
                return lambda_L
              elif P[k] == 0 and P[i] != 0:
                lambda_P = -P[j] / P[i]
                if lambda_L == lambda_P:
                  return lambda_L
            elif L_vars == {x_k}:
              if P[i] != 0 and P[k] == 0:
                return -P[j] / P[i]

    return None


  def reduced_subscheme(self):
    r"""
    Return the reduced subscheme of `self` as a projective plane curve.

    EXAMPLES:
      sage: R.<x0,x1,x2> = GF(2^3)[]
      sage: f = (x0^2 + x1*x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^2*x2^2
      sage: X.reduced_subscheme()
      Projective Plane Curve with defining polynomial x0^2 + x1*x2

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.reduced_subscheme()
      Projective Plane Curve with defining polynomial x0 + x1 + x2

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.reduced_subscheme()
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
    """

    f = prod(factor for factor, multiplicity in self._decompose)

    return ProjectivePlaneCurve(f)


  def irreducible_components(self):
    r"""
    Return the list of irreducible components of `self`.

    EXAMPLES:
      sage: R.<x0,x1,x2> = GF(2^3)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.irreducible_components()
      [Projective Plane Curve with defining polynomial x0 + x1 + x2]

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.irreducible_components()
      [Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4]

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0 + x1)*(x1 - x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1 + x1^2 - x0*x2 - x1*x2
      sage: X.irreducible_components()
      [Projective Plane Curve with defining polynomial -x1 + x2,
       Projective Plane Curve with defining polynomial x0 + x1]
    """

    return [ProjectivePlaneCurve(factor)
            for factor, multiplicity in self._decompose]


  def nonreduced_components(self):
    r"""
    Return the list of nonreduced components with corresponding
    multiplicities of `self`.

    OUTPUT:
    A list of tuples `(m, Y)` where `Y` is an irreducible component
    contained in `self` with multiplicity `m` such that m > 1.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0 + x1)*x2^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x2^2 + x1*x2^2
      sage: X.nonreduced_components()
      [(2, Projective Plane Curve with defining polynomial x2)]

      sage: R.<x0,x1,x2> = GF(2^3)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.nonreduced_components()
      [(4, Projective Plane Curve with defining polynomial x0 + x1 + x2)]

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^3 * (x1 + x2)^2 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3*x1^2*x2 + 2*x0^3*x1*x2^2 + x0^3*x2^3
      sage: X.nonreduced_components()
      [(2, Projective Plane Curve with defining polynomial x1 + x2),
       (3, Projective Plane Curve with defining polynomial x0)]

      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0 + x1)*(x1 - x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1 + x1^2 - x0*x2 - x1*x2
      sage: X.nonreduced_components()
      []
    """

    return [(multiplicity, ProjectivePlaneCurve(factor))
            for factor, multiplicity in self._decompose
            if multiplicity > 1]


  def singular_points(self):
    r"""
    Return the list of singular points of `self`.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.singular_points()
      []

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.singular_points()
      [(0 : 1 : 1), (1 : 0 : 1), (1 : 1 : 0)]
    """

    return self.plane_curve.singular_points()


  def multiplicity(self, P):
    r"""
    Return the multiplicity of `self` at the point `P`, i.e. the
    degree of the tangent cone at `P`.

    INPUT:
    - ``P`` -- a point on the projective plane.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: P = [0,0,1]
      sage: f = x0*x2^2 + x0^2*x1
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x1 + x0*x2^2
      sage: X.multiplicity(P)
      1
      sage: 
      sage: f = x0^2*x2^2 + x0^2*x1^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x1^2 + x0^2*x2^2
      sage: X.multiplicity(P)
      2
      sage: 
      sage: f = x0^3*x2 + x0^2*x1^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x1^2 + x0^3*x2
      sage: X.multiplicity(P)
      3
    """

    return self.plane_curve.multiplicity(P)


  def maximal_multiplicity(self):
    r"""
    Return the maximum of all multiplicities of rational points
    on `self`.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0 * (x1 + x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1^2 + 2*x0*x1*x2 + x0*x2^2
      sage: X.maximal_multiplicity()
      3
      sage:
      sage: X.multiplicity([1,1,-1])
      2
      sage: X.multiplicity([0,1,-1])
      3

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^3 * (x1 + x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2
      sage: X.maximal_multiplicity()
      4
      sage:
      sage: max(X.multiplicity(P) for P in X.singular_points())
      4

    MATHEMATICAL INTERPRETATION:
    The maximal multiplicity of the scheme defined by self.polynomial at
    a rational point P is sought. This occurs either:
    (1) At a singular point P of the support (reduced subscheme).
    (2) At a generic (smooth) point P of an irreducible component of the support.
        If self.polynomial = ... * factor_i^e_i * ..., the multiplicity of self
        at a generic point of the component defined by factor_i is e_i.
    The method computes the maximum over all such values.
    """

    X_red_sing = self._reduced_singular_points
    max_sing_mult = max((self.multiplicity(P) for P in X_red_sing), default=1)

    component_mults = [mult for mult, comp in self.nonreduced_components()]
    max_comp_mult = max(component_mults, default=1)

    return max(max_sing_mult, max_comp_mult)



  def singular_locus_dimension(self):
    r"""
    Return the dimension of the singular locus of `self`.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0*(x1 + x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1^2 + 2*x0*x1*x2 + x0*x2^2
      sage: X.singular_locus_dimension()
      1
      sage:
      sage: f = x0^2*x2 - x1^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x1^3 + x0^2*x2
      sage: X.singular_locus_dimension()
      0
      sage:
      sage: f = x0^4 + x1^4 + x2^4
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^4 + x2^4
      sage: X.singular_locus_dimension()
      -1
    """

    if self.is_smooth():
      return -1
    if self.is_reduced():
      return 0
    return 1


  def maximal_multiplicity_points(self): # upgrade to infinite fields: add nonreduced components to the list
    r"""
    Return the list of points of maximal multiplicity.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0*(x1 + x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1^2 + 2*x0*x1*x2 + x0*x2^2
      sage: X.maximal_multiplicity_points()
      [(0 : -1 : 1)]

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^3 * (x1 + x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2
      sage: X.maximal_multiplicity_points()
      [(0 : 1 : 1)]
      sage:
      sage: f = (x0^2 + x1*x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^2*x2^2
      sage: X.maximal_multiplicity_points()
      [(0 : 0 : 1), (0 : 1 : 0), (1 : 1 : 1)]
    """

    max_mult = self.maximal_multiplicity()
    points_with_max_multiplicity = []

    if self.base_ring.is_finite():
      for P in self.singular_points():
        if self.multiplicity(P) == max_mult:
          points_with_max_multiplicity.append(P)
      return points_with_max_multiplicity

    if any(multiplicity == max_mult for multiplicity, component in self.nonreduced_components()):
        raise NotImplementedError

    for P in self.reduced_subscheme().singular_points():
      if self.multiplicity(P) == max_mult:
        points_with_max_multiplicity.append(P)

    return points_with_max_multiplicity


  def points_with_high_multiplicity(self): # upgrade to infinite fields: add nonreduced components to the list
    r"""
    Return a list of points on `self` with multiplicity strictly
    greater than self.degree() / 2.

    OUTPUT:
    A list of tuples `(P, m)` where `P` is a point contained in `self` with
    multiplicity `m` such that
    m > self.degree() / 2.

    EXAMPLES:
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^2 * (x1 + x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x1^2 + x0^2*x2^2
      sage: X.points_with_high_multiplicity()
      [((0 : 1 : 1), 4)]
      sage: 
      sage: f = (x0^2 + x1*x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^2*x2^2
      sage: X.points_with_high_multiplicity()
      []
      sage:
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0 * (x1 + x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1^2 + x0*x2^2
      sage: X.points_with_high_multiplicity()
      [((0 : 1 : 1), 3), ((1 : 0 : 0), 2), ((1 : 1 : 1), 2)]
    """

    L = []
    for P in self.plane_curve.singular_points():
      m = self.plane_curve.multiplicity(P)
      if m > self.degree() / Integer(2):
        L.append((P, m))

    return L


  def lines_with_high_multiplicity(self):
    r"""
    Return a list of lines contained in `self`.

    OUTPUT:
    A list of triples `(L, m, G)` where `L` is a line contained
    in `self` with multiplicity `m` such that either
    0 < m <= self.degree() / 3 and self.polynomial = L^m * G
    or
    G = None.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0 * (x1 + x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1^2 + 2*x0*x1*x2 + x0*x2^2
      sage: X.lines_with_high_multiplicity()
      [(x0, 1, x1^2 + 2*x1*x2 + x2^2), (x1 + x2, 2, None)]

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^2 * x1^3 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x1^3*x2
      sage: X.lines_with_high_multiplicity()
      [(x2, 1, x0^2*x1^3), (x0, 2, x1^3*x2), (x1, 3, None)]
      sage:
      sage: f = x0
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0
      sage: X.lines_with_high_multiplicity()
      [(x0, 1, None)]
    """

    L = []
    polynomial_factors = list(self.polynomial.factor())
    for i, (factor, factor_multiplicity) in enumerate(polynomial_factors):
      if factor.degree() > 1:
        continue

      if factor_multiplicity > self.degree() / Integer(3):
        L.append((factor, factor_multiplicity, None))
      else:
        G = Integer(1)
        for j, (Gfactor, Gfactor_multiplicity) in enumerate(polynomial_factors):
          if j != i:
            G = G * Gfactor**Gfactor_multiplicity
        L.append((factor, factor_multiplicity, G))

    return L


  def _pseudo_instabilities(self): # Adjust nonreduced line
    r"""
    Return a list of pseudo-instabilities of `self`.

    EXAMPLES:
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = (x0^2 + x1*x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x1^2*x2^2
      sage: X.pseudo_instabilities()
      []

      sage: R.<x0,x1,x2> = GF(3)[]
      sage: f = x0^3 + x1^3 + x0^2*x2 - x0*x2^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x1^3 + x0^2*x2 - x0*x2^2
      sage: X.pseudo_instabilities()
      [Pseudo-instability of Projective Plane Curve with defining polynomial x0^3 + x1^3 + x0^2*x2 - x0*x2^2 given by [2, 2, 1] and x0 + x2]

      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^3 * (x1 + x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2
      sage: X.pseudo_instabilities()
      [Pseudo-instability of Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by [0, 0, 1],
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by [0, 1, 0],
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by [0, 1, 1],
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by [0, 1, 1] and x1 + x2,
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by x0]

      sage: R.<x0,x1,x2> = GF(3)[]
      sage: f = (x0 - x1)^2 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2
      sage: X.pseudo_instabilities()
      [Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [0, 0, 1] and x0 - x1,
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [1, 1, 0],
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [1, 1, 1] and x0 - x1,
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [2, 2, 1] and x0 - x1,
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [1, 1, 0] and x2,
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by x0 - x1]

    MATHEMATICAL INTERPRETATION:
    Give reference.
    """

    list_of_pseu_inst = []
    # CASES (b) and (d)
    for P, m in self.points_with_high_multiplicity():
      if m > 2 * self.degree() / 3:  # CASE (b)
        list_of_pseu_inst.append(ProjectiveFlag(P, None))
      elif m > self.degree() / 2:  # CASE (d) 1/2
        for L, L_multiplicity in PPC_TangentCone(self, P).embedded_lines():
          if L_multiplicity > m / 2:
            list_of_pseu_inst.append(ProjectiveFlag(P, L))

    # CASES (a) and (c)
    for L, L_multiplicity, G in self.lines_with_high_multiplicity():
      if L_multiplicity > self.degree() / Integer(3):  # CASE (a)
        list_of_pseu_inst.append(ProjectiveFlag(None, L))
      else: # CASE (c) 1/2
        L_curve = self.projective_plane.curve(L)
        G_curve = self.projective_plane.curve(G)
        for P in L_curve.intersection_points(G_curve):
          if L_curve.intersection_multiplicity(G_curve, P) > (self.degree() - L_multiplicity) / Integer(2):
            list_of_pseu_inst.append(ProjectiveFlag(P, L))

    return list_of_pseu_inst


  def flags(self):
    r"""
    Return a generator object of flags attached to `self`.

    EXAMPLES:
    A properly semistable quartic.
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0^2 + x1*x2)^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + 2*x0^2*x1*x2 + x1^2*x2^2
      sage: list(X.flags())
      []

    A properly semistable cubic.
      sage: R.<x0,x1,x2> = GF(3)[]
      sage: f = x0^3 + x1^3 + x0^2*x2 - x0*x2^2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x1^3 + x0^2*x2 - x0*x2^2
      sage: list(X.flags())
      [Flag attached to Projective Plane Curve with defining polynomial x0^3 + x1^3 + x0^2*x2 - x0*x2^2 given by [2, 2, 1] and x0 + x2]

    An unstable quartic.
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^3 * (x1 + x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2
      sage: list(X.flags())
      [Flag attached to Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by [0, 1, 1],
       Flag attached to Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by [0, 1, 1] and x1 + x2,
       Flag attached to Projective Plane Curve with defining polynomial x0^3*x1 + x0^3*x2 given by x0]

    An unstable cubic.
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = (x0 - x1)^2 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x2 - 2*x0*x1*x2 + x1^2*x2
      sage: list(X.flags())
      [Flag attached to Projective Plane Curve with defining polynomial x0^2*x2 - 2*x0*x1*x2 + x1^2*x2 given by [1, 1, 0],
       Flag attached to Projective Plane Curve with defining polynomial x0^2*x2 - 2*x0*x1*x2 + x1^2*x2 given by [1, 1, 0] and x2,
       Flag attached to Projective Plane Curve with defining polynomial x0^2*x2 - 2*x0*x1*x2 + x1^2*x2 given by x0 - x1]

    MATHEMATICAL INTERPRETATION:
    Give reference.
    """

    X_red_sing = self._reduced_singular_points
    # CASES (b) and (d)
    for P in X_red_sing:
      m = self.multiplicity(P)
      if m > 2 * self.degree() / 3:  # CASE (b)
        yield ProjectiveFlag(P, None)
      elif m > self.degree() / 2:  # CASE (d) 1/2
        for L, L_multiplicity in PPC_TangentCone(self, P).embedded_lines():
          if L_multiplicity > m / 2:
            yield ProjectiveFlag(P, L)

    # CASES (a) and (c)
    for L, L_multiplicity, G in self.lines_with_high_multiplicity():
      if L_multiplicity > self.degree() / 3:  # CASE (a)
        yield ProjectiveFlag(None, L)
      else: # CASE (c) 1/2
        L_curve = self.projective_plane.curve(L)
        G_curve = self.projective_plane.curve(G)
        for P in L_curve.intersection_points(G_curve):
          if L_curve.intersection_multiplicity(G_curve, P) > (self.degree() - L_multiplicity) / 2:
            yield ProjectiveFlag(P, L)


  def instabilities(self):
    r"""
    Return a list of pseudo-instabilities of `self` which are
    instabilities.

    EXAMPLES:
      sage: R.<x0,x1,x2> = GF(2)[]
      sage: f = x0^3 + x1^2 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^3 + x1^2*x2
      sage: X.instabilities()
      [Pseudo-instability of Projective Plane Curve with defining polynomial x0^3 + x1^2*x2 given by [0, 0, 1] and x1]
      sage: len(X.pseudo_instabilities())
      1

    There can be less instabilities than pseudo-instabilities as the following example shows.
      sage: R.<x0,x1,x2> = GF(3)[]
      sage: f = (x0 - x1)^2 * x2
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2
      sage: X.instabilities()
      [Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [0, 0, 1] and x0 - x1,
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [1, 1, 0],
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by [1, 1, 0] and x2,
       Pseudo-instability of Projective Plane Curve with defining polynomial x0^2*x2 + x0*x1*x2 + x1^2*x2 given by x0 - x1]
      sage: len(X.pseudo_instabilities())
      6

    It can happen that none of the pseudo-instabilities is an instability.
      sage: R.<x0,x1,x2> = GF(2^8)[]
      sage: f = x0 * (x0^3 + x1^2*x2 + x0*x2^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^4 + x0*x1^2*x2 + x0^2*x2^2
      sage: X.instabilities()
      []
      sage: len(X.pseudo_instabilities())
      1
      sage:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0 * x2 * (x1^2 + x0*x2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0*x1^2*x2 + x0^2*x2^2
      sage: X.instabilities()
      []
      sage: len(X.pseudo_instabilities())
      2
    """

    return [I for I in self.flags() if I.is_unstable()]


  @cached_property
  def _decompose(self):
    r"""
    Return the factored form of self.polynomial.

    EXAMPLES:
      sage: R.<x0,x1,x2> = QQ[]
      sage: f = x0 * x1^2 * (x0 * x1 + x2^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial x0^2*x1^3 + x0*x1^2*x2^2
      sage: X._decompose
      [(x0, 1), (x1, 2), (x0*x1 + x2^2, 1)]
    """

    return list(self.polynomial.factor())


  @cached_property
  def _reduced_singular_points(self):
    r"""
    Return the list of singular points of the reduced subscheme
    of `self`.
    """

    return self.reduced_subscheme().singular_points()



class PPC_TangentCone:
  r"""
  Construct the tangent cone of a projective plane curve at a point
  to the following conditions.

  INPUT:
  - ``projective_plane_curve`` -- a projective plane curve.
  - ``P`` -- a point in the projective plane.

    EXAMPLES:
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3 - (x-2*z)^2*z
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3
      sage: P = [2,1,1]
      sage: PPC_TangentCone(X, P)
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3 at [2, 1, 1]

      sage: R.<x,y,z> = QQ[]
      sage: f = (y*z - x^2)*(y*z + x^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^4 + y^2*z^2
      sage: P = [0,0,1]
      sage: PPC_TangentCone(X, P)
      Tangent cone of Projective Plane Curve with defining polynomial -x^4 + y^2*z^2 at [0, 0, 1]
  """

  def __init__(self, projective_plane_curve, P):
    r"""
    Construct the tangent cone of a projective plane curve at a point
    to the following conditions.

    INPUT:
    - ``projective_plane_curve`` -- a projective plane curve.
    - ``P`` -- a point in the projective plane.

    EXAMPLES:
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3 - (x-2*z)^2*z
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3
      sage: P = [2,1,1]
      sage: PPC_TangentCone(X, P)
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3 at [2, 1, 1]

      sage: R.<x,y,z> = QQ[]
      sage: f = (y*z - x^2)*(y*z + x^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^4 + y^2*z^2
      sage: P = [0,0,1]
      sage: PPC_TangentCone(X, P)
      Tangent cone of Projective Plane Curve with defining polynomial -x^4 + y^2*z^2 at [0, 0, 1]
    """

    # Convert to list
    P = list(P)
    if projective_plane_curve.get_polynomial()(P) != 0:
      raise ValueError

    self.projective_plane_curve = projective_plane_curve
    self.base_ring = projective_plane_curve.get_base_ring()
    self.polynomial_ring = PolynomialRing(self.base_ring, ['x', 'y'])
    self.gen1, self.gen2 = self.polynomial_ring.gens()

    # Coerce coordinates to self.base_ring
    for i in range(3):
      P[i] = self.base_ring(P[i])

    # Find the maximal index, i_max, with P[i_max] != 0 and normalize by P[i_max]
    self.affine_patch, self.normalized_point = _normalize_by_last_nonzero_entry(P)


  def __repr__(self):
    return f"Tangent cone of {self.projective_plane_curve} at {self.normalized_point}"


  def get_polynomial(self):
    r"""
    Return the defining polynomial of `self`.

    EXAMPLES:
      A nodal cubic with node at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3 - (x-2*z)^2*z
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3 at [2, 1, 1]
      sage: h = C.get_polynomial(); h
      -x^2 + y^2
      sage: h.parent()
      Multivariate Polynomial Ring in x, y over Rational Field
      sage: h.factor()
      (-1) * (x - y) * (x + y)

    A cuspidal cubic with cusp at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3 at [2, 1, 1]
      sage: h = C.get_polynomial(); h
      y^2
      sage: h.parent()
      Multivariate Polynomial Ring in x, y over Rational Field

    A quartic with tacnode at (0:0:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y*z - x^2)*(y*z + x^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^4 + y^2*z^2
      sage: P = [0,0,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^4 + y^2*z^2 at [0, 0, 1]
      sage: h = C.get_polynomial(); h
      y^2
      sage: h.parent()
      Multivariate Polynomial Ring in x, y over Rational Field
    """

    PPC_equation = self.projective_plane_curve.get_polynomial()
    dehomogenization = [self.gen1, self.gen2]
    dehomogenization.insert(self.affine_patch, self.polynomial_ring(0))
    dehomogenization_translated = []

    for i in range(3):
      dehomogenization_translated.append(dehomogenization[i] + self.normalized_point[i])

    f = PPC_equation(dehomogenization_translated)
    f_homo_comp_dict = f.homogeneous_components()
    minimal_degree = min(f_homo_comp_dict.keys())
    tangent_cone_polynomial = f_homo_comp_dict[minimal_degree]

    return tangent_cone_polynomial


  def embedded_polynomial(self):
    r"""
    Return the defining polynomial of the embedding of `self` into
    the projective plane at the point `self.normalized_point`.

    EXAMPLES:
    A nodal cubic with node at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3 - (x-2*z)^2*z
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3 at [2, 1, 1]
      sage: h = C.embedded_polynomial(); h
      -x^2 + y^2 + 4*x*z - 2*y*z - 3*z^2
      sage: h(P)
      0
      sage: h.factor()
      (-x + y + z) * (x + y - 3*z)
      sage: h(x + 2, y + 1, 1)
      -x^2 + y^2
      sage: C.get_polynomial()
      -x^2 + y^2

    A cuspidal cubic with cusp at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3 at [2, 1, 1]
      sage: h = C.embedded_polynomial(); h
      y^2 - 2*y*z + z^2
      sage: h(P)
      0
      sage: h.factor()
      (y - z)^2
      sage: h(x + 2, y + 1, 1)
      y^2
      sage: C.get_polynomial()
      y^2

    A quartic with tacnode at (0:0:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y*z - x^2)*(y*z + x^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^4 + y^2*z^2
      sage: P = [0,0,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^4 + y^2*z^2 at [0, 0, 1]
      sage: h = C.embedded_polynomial(); h
      y^2
      sage: h(P)
      0
      sage: C.get_polynomial()
      y^2

    MATHEMATICAL INTERPRETATION:
    First, let
      F = self.projective_plane_curve.get_polynomial(),
      K = self.projective_plane_curve.get_base_ring(),
      P = self.normalized_point,
      j = self.affine_patch,
      x_0, x_1, x_2 = self.projective_plane_curve.get_standard_basis().
    Then P[j] = 1 and for all j < i <= 2 we have
      P[i] = 0.
    Further, we set y_0 = x_0, y_1 = x_1, y_2 = x_2, y_j = 1 and 
      f = F(y_0 + P[0]*y_j, y_1 + P[1]*y_j, y_2 + P[2]*y_j).
    This is the dehomogenization of F to the affine patch j and the
    subsequent translation of the point P to the origin in this affine
    patch. Thus,
      f = _apply_matrix(T, F, affine_patch = j)
    with
      T = _ult_line_transformation(K, P).
    Now let TanCon be the homogeneous component of f of the lowest degree,
    i.e. TanCon is the homogeneous polynomial defining the tangent cone of
    self.projective_plane_curve at the point P. We view TanCon as a polynomial
    in K[x_0, x_1, x_2], although it does not depend on x_j. Let e_j be the
    j-th standard basis vector. Then we have
      TanCon(e_j) = 0 and e_j*T = P.
    Now we define
      TanCon_P = TanCon(x_0 - P[0]*x_j, x_1 - P[1]*x_j, x_2 - P[2]*x_j),
    i.e.
      TanCon_P = _apply_matrix(T.inverse(), TanCan).
    In particular,
      TanCon_P(P) = TanCon(e_j) = 0.
    Moreover, the dehomogenization of TanCon_P to the affine patch j is the
    translation of TanCon from the origin to the point P.            
    """

    T = _ult_line_transformation(self.base_ring, self.normalized_point)
    T_inverse = T.inverse()
    F = self.projective_plane_curve.get_polynomial()
    f = _apply_matrix(T, F, self.affine_patch)
    f_homo_comp_dict = f.homogeneous_components()
    minimal_degree = min(f_homo_comp_dict.keys())
    tangent_cone_polynomial = f_homo_comp_dict[minimal_degree]

    return _apply_matrix(T.inverse(), tangent_cone_polynomial)


  def get_lines(self):
    r"""
    Return linear factors of the defining polynomial of `self`.

    OUTPUT:
    A list of tuples `(L, m)` where `L` is a line contained
    in `self` with multiplicity `m`.

    EXAMPLES:
    A nodal cubic with node at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3 - (x-2*z)^2*z
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3 at [2, 1, 1]
      sage: C.get_lines()
      [(x - y, 1), (x + y, 1)]
      sage: h = C.get_polynomial(); h
      -x^2 + y^2
      sage: h.factor()
      (-1) * (x - y) * (x + y)

    A cuspidal cubic with cusp at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3 at [2, 1, 1]
      sage: C.get_lines()
      [(y, 2)]
      sage: h = C.get_polynomial(); h
      y^2

    A quartic with tacnode at (0:0:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y*z - x^2)*(y*z + x^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^4 + y^2*z^2
      sage: P = [0,0,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^4 + y^2*z^2 at [0, 0, 1]
      sage: C.get_lines()
      [(y, 2)]
      sage: h = C.get_polynomial(); h
      y^2
    """

    L = []
    tangent_cone_polynomial = self.get_polynomial()
    for factor, factor_multiplicity in list(tangent_cone_polynomial.factor()):
      if factor.degree() == 1:
        L.append((factor, factor_multiplicity))

    return L


  def embedded_lines(self):
    r"""
    Return linear factors of the polynomial self.embedded_polynomial().

    OUTPUT:
    A list of tuples `(L, m)` where `L` is a line contained with
    multiplicity `m`in the embedding of `self` into the projective
    plane at the point `self.normalized_point`.

    EXAMPLES:
    A nodal cubic with node at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3 - (x-2*z)^2*z
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 5*x^2*z + y^2*z - 8*x*z^2 - 2*y*z^2 + 5*z^3 at [2, 1, 1]
      sage: C.embedded_lines()
      [(-x + y + z, 1), (x + y - 3*z, 1)]
      sage: C.get_lines()
      [(x - y, 1), (x + y, 1)]

    A cuspidal cubic with cusp at (2:1:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y-z)^2*z - (x-2*z)^3
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3
      sage: P = [2,1,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^3 + 6*x^2*z + y^2*z - 12*x*z^2 - 2*y*z^2 + 9*z^3 at [2, 1, 1]
      sage: C.embedded_lines()
      [(y - z, 2)]
      sage: C.get_lines()
      [(y, 2)]

    A quartic with tacnode at (0:0:1).
      sage: R.<x,y,z> = QQ[]
      sage: f = (y*z - x^2)*(y*z + x^2)
      sage: X = ProjectivePlaneCurve(f); X
      Projective Plane Curve with defining polynomial -x^4 + y^2*z^2
      sage: P = [0,0,1]
      sage: C = PPC_TangentCone(X, P); C
      Tangent cone of Projective Plane Curve with defining polynomial -x^4 + y^2*z^2 at [0, 0, 1]
      sage: C.embedded_lines()
      [(y, 2)]
      sage: C.get_lines()
      [(y, 2)]
    """

    L = []
    tangent_cone_polynomial = self.embedded_polynomial()
    for factor, factor_multiplicity in list(tangent_cone_polynomial.factor()):
      if factor.degree() == 1:
        L.append((factor, factor_multiplicity))

    return L


class ProjectiveFlag:
  r"""
  Construct a projective flag to the following conditions.

  INPUT:
  - ``projective_point`` -- a point in the projective plane.
  - ``linear_form`` -- a linear form defining a line in the projective plane.
  """

  def __init__(self, projective_point=None, linear_form=None):
    r"""
    Construct a projective flag to the following conditions.

    INPUT:
    - ``projective_point`` -- a point in the projective plane.
    - ``linear_form`` -- a linear form defining a line in the projective plane.

    EXAMPLES:
      sage: ProjectiveFlag([1,2,3])
      Projective flag given by [1, 2, 3]

      sage: R.<x0,x1,x2> = QQ[]
      sage: ProjectiveFlag(linear_form = x0 + 2*x1 - x2)
      Projective flag given by x0 + 2*x1 - x2
      sage:
      sage: ProjectiveFlag([1,1,3], x0 + 2*x1 - x2)
      Projective flag given by [1, 1, 3] and x0 + 2*x1 - x2
    """

    if projective_point is None and linear_form is None:
      raise ValueError("Both arguments are None")

    if projective_point is not None:
      proj_P_list = list(projective_point)
      if linear_form is not None and linear_form(proj_P_list) != 0:
        raise ValueError(
          f"{projective_point} is not a point on the projective line given by {linear_form}")

    if projective_point is not None:
      self.point = proj_P_list
    else:
      self.point = None
    self.line = linear_form


  def __repr__(self):
    if self.point == None:
      return f"Projective flag given by {self.line}"
    elif self.line == None:
      return f"Projective flag given by {self.point}"
    else:
      return f"Projective flag given by {self.point} and {self.line}"


  def base_change_matrix(self, matrix_form = 'uut'):
    r"""
    Return a unipotent matrix transforming a flag given by some
    standard basis vector e_j and some line x_i = 0 to `self`.

    INPUT:
    - ``matrix_form`` -- a list of rational numbers or one of the strings 'ult', 'uut'.

    OUTPUT:
    A unipotent matrix.

    EXAMPLES:
      sage: K.<a,b,c,A,B,C> = QQ[]
      sage: K = K.fraction_field()
      sage: R.<x0,x1,x2> = K[]
      sage: P = [a,b,c]
      sage: L = A*x0 + B*x1 - (a*A/c + b*B/c)*x2
      sage: L(P)
      0
      sage: F = ProjectiveFlag(P, L)
      sage: T = F.base_change_matrix('ult')
      sage: T = F.base_change_matrix('ult'); T
      [     1      0      0]
      [(-B)/A      1      0]
      [   a/c    b/c      1]
      sage: c * vector([0,0,1]) * T
      (a, b, c)
      sage: L(list(vector([x0,x1,x2]) * T))
      A*x0
      sage:
      sage: P = [a,b,0]
      sage: L = A*x0 - (a*A/b)*x1 + C*x2
      sage: L(P)
      0
      sage: F = ProjectiveFlag(P, L)
      sage: T = F.base_change_matrix('ult'); T
      [     1      0      0]
      [   a/b      1      0]
      [(-C)/A      0      1]
      sage: b * vector([0,1,0]) * T
      (a, b, 0)
      sage: L(list(vector([x0,x1,x2]) * T))
      A*x0
      sage:
      sage: P = [a,0,0]
      sage: L = B*x1 + C*x2
      sage: L(P)
      0
      sage: F = ProjectiveFlag(P, L)
      sage: T = F.base_change_matrix('ult'); T
      [     1      0      0]
      [     0      1      0]
      [     0 (-C)/B      1]
      sage: a * vector([1,0,0]) * T
      (a, 0, 0)
      sage: L(list(vector([x0,x1,x2]) * T))
      B*x1
      sage:
      sage: P = [a,b,c]
      sage: L = -(b*B/a + c*C/a)*x0 + B*x1 + C*x2
      sage: L(P)
      0
      sage: F = ProjectiveFlag(P, L)
      sage: T = F.base_change_matrix('uut'); T
      [     1    b/a    c/a]
      [     0      1 (-B)/C]
      [     0      0      1]
      sage: a * vector([1,0,0]) * T
      (a, b, c)
      sage: L(list(vector([x0,x1,x2]) * T))
      C*x2
      sage:
      sage: P = [0,b,c]
      sage: L = A*x0 - (c*C/b)*x1 + C*x2
      sage: L(P)
      0
      sage: F = ProjectiveFlag(P, L)
      sage: T = F.base_change_matrix('uut'); T
      [     1      0 (-A)/C]
      [     0      1    c/b]
      [     0      0      1]
      sage: b * vector([0,1,0]) * T
      (0, b, c)
      sage: L(list(vector([x0,x1,x2]) * T))
      C*x2
      sage:
      sage: P = [0,0,c]
      sage: L = A*x0 + B*x1
      sage: L(P)
      0
      sage: F = ProjectiveFlag(P, L)
      sage: T = F.base_change_matrix('uut'); T
      [     1 (-A)/B      0]
      [     0      1      0]
      [     0      0      1]
      sage: c * vector([0,0,1]) * T
      (0, 0, c)
      sage: L(list(vector([x0,x1,x2]) * T))
      B*x1
    """

    if self.point is None:
      if matrix_form == 'uut':
        return _uut_plane_transformation(self.line)
      elif matrix_form == 'ult':
        return _ult_plane_transformation(self.line)
      elif isinstance(matrix_form, list):
        return _integral_plane_transformation(self.line, matrix_form)
      else:
        raise ValueError
    elif self.line is None:
      base_ring = self.line.base_ring()
      if matrix_form == 'uut':
        return _uut_line_transformation(base_ring, self.point)
      elif matrix_form == 'ult':
        return _ult_line_transformation(base_ring, self.point)
      elif isinstance(matrix_form, list):
        return _integral_line_transformation(base_ring, self.point, matrix_form)
      else:
        raise ValueError
    else:
      if matrix_form == 'uut':
        return _uut_flag_transformation(self.point, self.line)
      elif matrix_form == 'ult':
        return _ult_flag_transformation(self.point, self.line)
      elif isinstance(matrix_form, list):
        return _integral_flag_transformation(self.point, self.line, matrix_form)
      else:
        raise ValueError


  def is_unstable(self, projective_plane_curve):
    r"""
    Return `True` or `False` depending on whether the base
    change matrix `self.base_change_matrix()` gives rise to
    an instability of `projective_plane_curve`.

    INPUT:
    - ``projective_plane_curve`` -- a projective plane curve.

    MATHEMATICAL INTERPRETATION:
    First, let
      K = self.base_ring,
      T = self.flag_transformation(),
      F = self.proj_plane_curve.get_polynomial().
    Furthermore, let
      (x0, x1, x2) = self.proj_plane_curve.get_standard_basis()
    and
      G = F((x0,x1,x2)*T),
    i.e.
      G = _apply_matrix(T, F).
    For a multi index set I subset NN^3 we can write
      G = sum_{i in I} a_i x0^i0 * x1^i1 * x2^i2
    with a_i != 0 for all i in I. Note that with respect to the new
    basis (x0,x1,x2)*T the flag (self.point, self.line) is given by
    (e_j, x_i). Thus, it yields an instability, if there exists a
    balanced weight vector
      (w0, w1, w2) in QQ^3, i.e. w0 + w1 + w2 = 0,
    such that
      i0*w0 + i1*w1 + i2*w2 > 0
    for all i in I. 
    REMARK. Any nonzero multiple of a balanced weight vector is again
    a balanced weight vector. Thus, it suffices to consider
      (w0, w1, w2) in QQ^3
    wtih
      -1 <= w0, w1, w2 <= 1.
    Thus, we only have to maximize the function
      min(i0*w0 + i1*w1 + i2*w2 : i in I)
    under the constraints -1 <= w0, w1, w2 <= 1 and to check whether
    the maximum is > 0 or not.
    """

    if not isinstance(projective_plane_curve, ProjectivePlaneCurve):
      raise TypeError

    T = self.base_change_matrix()
    F = projective_plane_curve.get_polynomial()
    G = _apply_matrix(T, F)

    MILP = MixedIntegerLinearProgram(solver='PPL')

    v = MILP.new_variable()
    t = v['minimum']
    w0 = v['w0']
    w1 = v['w1']
    w2 = v['w2']

    MILP.set_objective(t)
    MILP.add_constraint(-1 <= w0 <= 1)
    MILP.add_constraint(-1 <= w1 <= 1)
    MILP.add_constraint(-1 <= w2 <= 1)
    MILP.add_constraint(w0 + w1 + w2 == 0)

    for i in G.dict():
      MILP.add_constraint(t <= i[0] * w0 + i[1] * w1 + i[2] * w2)

    MILP.solve()
    values = MILP.get_values(v)

    return values['minimum'] > 0


  def is_semiunstable(self, projective_plane_curve):
    r"""
    Return `True` or `False` depending on whether the base
    change matrix `self.base_change_matrix()` gives rise to
    an semiinstability of `projective_plane_curve`.

    INPUT:
    - ``projective_plane_curve`` -- a projective plane curve.

    MATHEMATICAL INTERPRETATION:
    First, let
      K = self.base_ring,
      T = self.flag_transformation(),
      F = self.proj_plane_curve.get_polynomial().
    Furthermore, let
      (x0, x1, x2) = self.proj_plane_curve.get_standard_basis()
    and
      G = F((x0,x1,x2)*T),
    i.e.
      G = _apply_matrix(T, F).
    For a multi index set I subset NN^3 we can write
      G = sum_{i in I} a_i x0^i0 * x1^i1 * x2^i2
    with a_i != 0 for all i in I. Note that with respect to the new
    basis (x0,x1,x2)*T the flag (self.point, self.line) is given by
    (e_j, x_i). Thus, it yields an semiinstability, if there exists a
    balanced weight vector
      (w0, w1, w2) in QQ^3, i.e. w0 + w1 + w2 = 0 and w != (0,0,0),
    such that
      min(i0*w0 + i1*w1 + i2*w2 : i in I) = 0.
    REMARK. Any nonzero multiple of a balanced weight vector is again
    a balanced weight vector. Thus, it suffices to consider
      (w0, w1, w2) in QQ^3
    wtih
      ||w||_1 = 1.
    Thus, we only have to maximize the function
      min(i0*w0 + i1*w1 + i2*w2 : i in I)
    on the boundary of the cube [-1, 1]^3 and to check whether the
    maximum is 0 or not.
    """

    if not isinstance(projective_plane_curve, ProjectivePlaneCurve):
      raise TypeError

    T = self.base_change_matrix()
    F = projective_plane_curve.get_polynomial()
    G = _apply_matrix(T, F)

    maximum_is_zero = False
    # positive faces
    for position in range(3):
      for plus_minus in [Integer(-1), Integer(1)]:
        MILP = MixedIntegerLinearProgram(solver='PPL')
        v = MILP.new_variable()
        t = v['minimum']
        MILP.set_objective(t)
        # Conditions to be on
        # {1}x[-1,1]^2, [-1,1]x{1}x[-1,1], [-1,1]^2x{1},
        # {-1}x[-1,1]^2, [-1,1]x{-1}x[-1,1], [-1,1]^2x{-1}
        for i in range(3):
          if i == position:
            MILP.add_constraint(v[i] == plus_minus)
          MILP.add_constraint(-1 <= v[i] <= 1)
        # Condition to be in the zero space
        MILP.add_constraint(sum(v[i] for i in range(3)) == 0)
        # All linear functions are bounded by minimum.
        for exponent in G.exponents():
          lin_func = sum(Integer(i_j) * v[j] for j, i_j in enumerate(exponent))
          MILP.add_constraint(t <= lin_func)
        max_value = MILP.solve()
        if max_value > 0:
          return False
        elif max_value == 0:
          maximum_is_zero = True

    return maximum_is_zero



class PPC_Instability:

  def __init__(self, base_change_matrix, weight_vector, geometric_type = 'not specified'):

    self.base_change_matrix = base_change_matrix
    self.weight_vector = weight_vector
    self.geometric_type = geometric_type


  def get_base_change_matrix(self):
    return self.base_change_matrix


  def get_weight_vector(self):
    return self.weight_vector


  def get_geometric_type(self):
    return self.geometric_type



# ================== helper functions ==================

def _min_index_of_nonzero_entry(L):
  r"""
  Return minimal integer i between 0 and len(L) with L[i] != 0

  INPUT:
  L - list of rational numbers

  OUTPUT:
  i - minimal integer between 0 and len(L) with L[i] != 0
  """

  return next((i for i in range(len(L)) if L[i] != 0), None)


def _max_index_of_nonzero_entry(L):
  r"""
  Return maximal integer i between 0 and len(L) with L[i] != 0

  INPUT:
  L - list of rational numbers

  OUTPUT:
  i - maximal integer between 0 and len(L) with L[i] != 0
  """

  return next((i for i in reversed(range(len(L))) if L[i] != 0), None)


def _normalize_by_first_nonzero_entry(L):
  r"""
  Return the pair (i, [L[0] / L[i], ..., L[n] / L[i]]),
  where n = len(L) and i = _min_index_of_nonzero_entry(L).
  """

  i_min = _min_index_of_nonzero_entry(L)
  first_nonzero_entry = L[i_min]

  return (i_min, [x / first_nonzero_entry for x in L])


def _normalize_by_last_nonzero_entry(L):
  r"""
  Return the pair (i, [L[0] / L[i], ..., L[n] / L[i]]),
  where n = len(L) and i = _max_index_of_nonzero_entry(L).
  """

  i_max = _max_index_of_nonzero_entry(L)
  last_nonzero_entry = L[i_max]

  return (i_max, [x / last_nonzero_entry for x in L])


def _apply_matrix(T, F, i = None):
  r"""
  Return F((x_0,...,x_n) * T) or its dehomogenization
  at `i`, i.e. x_{i} = 1 if `i` is not `None`.

  INPUT:
  - ``T`` -- matrix over K.
  - ``F`` -- polynomial in K[x_0,...,x_n].
  - ``i`` -- an integer between 0 and n.

  OUTPUT:
  F((x_0,...,x_n) * T) with x_i = 1 if `i` is not `None`.

  EXAMPLES:
    sage: K.<t00,t01,t02,t10,t11,t12,t20,t21,t22> = QQ[]
    sage: R.<x0,x1,x2> = K[]
    sage: T = matrix(K, [[t00,t01,t02],[t10,t11,t12],[t20,t21,t22]]); T
    [t00 t01 t02]
    [t10 t11 t12]
    [t20 t21 t22]
    sage: _apply_matrix(T, x0)
    t00*x0 + t10*x1 + t20*x2
    sage: _apply_matrix(T, x1)
    t01*x0 + t11*x1 + t21*x2
    sage: _apply_matrix(T, x2)
    t02*x0 + t12*x1 + t22*x2

    sage: _apply_matrix(T, x0, 1)
    t00*x0 + t20*x2 + t10

    sage: _apply_matrix(T, x0 + x1)
    (t00 + t01)*x0 + (t10 + t11)*x1 + (t20 + t21)*x2

    sage: R.<x0,x1,x2> = QQ[]
    sage: F = x0*(x1 - 2*x0)*(x2 - 3*x0)
    sage: T = matrix(QQ, [[1,2,3],[0,1,0],[0,0,1]]); T
    [1 2 3]
    [0 1 0]
    [0 0 1]
    sage: _apply_matrix(T, F)
    x0*x1*x2

    sage: _apply_matrix(T, F, 0)
    x1*x2
    sage: _apply_matrix(T, F, 1)
    x0*x2
  """

  generators = list(F.parent().gens())
  if i != None:
    generators[i] = F.parent()(1)

  return F(list(vector(generators) * T))


def _ult_line_transformation(base_field, coordinates):
  r"""
  Return a unipotent lower triangular matrix over `base_field`
  transforming the line spanned by some standard basis vector
  to the line spanned by the vector defined by `coordinates`.

  EXAMPLES:
    sage: K.<a,b,c> = QQ[]
    sage: K = K.fraction_field()
    sage: T = _ult_line_transformation(K, [a,b,c]); T
    [  1   0   0]
    [  0   1   0]
    [a/c b/c   1]
    sage: c * vector([0,0,1]) * T
    (a, b, c)
    sage:
    sage: T = _ult_line_transformation(K, [a,b,0]); T
    [  1   0   0]
    [a/b   1   0]
    [  0   0   1]
    sage: b * vector([0,1,0]) * T
    (a, b, 0)
    sage:
    sage: T = _ult_line_transformation(K, [a,0,0]); T
    [1 0 0]
    [0 1 0]
    [0 0 1]
  """

  Vector = [base_field(x) for x in coordinates]
  T = identity_matrix(base_field, len(Vector))
  # Find the maximal index, i_max, with Vector[i_max] != 0
  # and normalize by Vector[i_max]
  i_max, normalized_Vector = _normalize_by_last_nonzero_entry(Vector)
  T[i_max] = normalized_Vector

  return matrix(base_field, T)


def _uut_line_transformation(base_field, coordinates):
  r"""
  Return a unipotent upper triangular matrix over `base_field`
  transforming the line spanned by some standard basis vector
  to the line spanned by the vector defined by `coordinates`.

  EXAMPLES:
    sage: K.<a,b,c> = QQ[]
    sage: K = K.fraction_field()
    sage: T = _uut_line_transformation(K, [a,b,c]); T
    [  1 b/a c/a]
    [  0   1   0]
    [  0   0   1]
    sage: a * vector([1,0,0]) * T
    (a, b, c)
    sage:
    sage: T = _uut_line_transformation(K, [0,b,c]); T
    [  1   0   0]
    [  0   1 c/b]
    [  0   0   1]
    sage: b * vector([0,1,0]) * T
    (0, b, c)
    sage:
    sage: T = _uut_line_transformation(K, [0,0,c]); T
    [1 0 0]
    [0 1 0]
    [0 0 1]
  """

  Vector = [base_field(x) for x in coordinates]
  T = identity_matrix(base_field, len(Vector))
  # Find the minimal index, i_min, with Vector[i_min] != 0
  # and normalize by Vector[i_min]
  i_min, normalized_Vector = _normalize_by_first_nonzero_entry(Vector)
  T[i_min] = normalized_Vector

  return matrix(base_field, T)


def _sorting_permutation_matrix(w):
  r"""
  Return permutation matrix, T, such that vector(w)*T is sorted in
  decreasing order

  INPUT:
  w - vector over Rational Field
  """

  n = len(w)
  if n == 0:
    return identity_matrix(0)

  # Create a list of (value, index) tuples
  indexed_w = [(val, i) for i, val in enumerate(w)]

  # Sort the list in descending order based on the values
  sorted_indexed_w = sorted(indexed_w, key=lambda x: x[0], reverse=True)

  # Extract the sorted indices
  sorted_indices = [index for _, index in sorted_indexed_w]

  # Create the permutation represented as a list.
  permutation_list = [0] * n
  for i, original_index in enumerate(sorted_indices):
    permutation_list[original_index] = i

  # Create the permutation matrix.
  permutation_matrix = zero_matrix(n)
  for i, j in enumerate(permutation_list):
    permutation_matrix[i, j] = 1

  return permutation_matrix


# def _ult_matrix_integralizator(ult_matrix, weight_vector):
#     """
#     Return unipotent matrix, T, which is the conjugation of ult_matrix
#     by a permutation matrix such that for all indices i,j the implication
#         (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
#     holds
# 
#     INPUT:
#         ult_matrix      - unipotent lower triangular matrix
#         weight_vector   - tuple of rational numbers
# 
#     OUTPUT:
#         T   - unipotent matrix such that there exists a permutation P
#               with T = P*ult_matrix*P.inverse() and weight_vector*P is
#               sorted in decreasing order, i.e.
#                 (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
# 
#     MATHEMATICAL INTERPRETATION:
#         First, let
#             A = ult_matrix,
#             w = [w_0,...,w_n] := weight_vector.
#         Let sigma be the permutation with
#             w_{sigma(0)} >= ... >= w_{sigma(n)}.
#         Let P be the matrix with columns given by
#             e_{sigma(0)},...,e_{sigma(n)},
#         where e_i is the i-th standard basis vector. Then we have
#             w * P = [w_{sigma(0)},...,w_{sigma(n)}],
#         i.e.
#             P = _sorting_permutation_matrix(w).
#         Let
#             v := w * P, i.e. v_j = w_{sigma(j)}.
#         Since v is sorted in decreasing order and T is a lower
#         triangular matrix, we have the implications
#             (v_j - v_i < 0) => (j > i) => (A[i][j] = 0).
#         Moreover, the matrix
#             P^{-1} = P.transpose()
#         corresponds to the permutation sigma^{-1}, i.e. the columns
#         of P^{-1} are given by
#             e_{sigma^{-1}(0)},...,e_{sigma^{-1}(n)}
#         and therefore the (i,j)-th entry of
#             A * P^{-1}
#         is equal to
#             A[i][sigma^{-1}(j)].
#         Further, since
#             P = (P^{-1})^{-1} = (P^{-1}).transpose(),
#         the rows of P are given by
#             e_{sigma^{-1}(0)},...,e_{sigma^{-1}(n)}.
#         Thus, the (i,j) entry of
#             T = P * T * P^{-1}
#         is given by
#             T[i][j] = A[sigma^{-1}(i)][sigma^{-1}(j)].
#         Finally, since
#             v_{sigma^{-1}(j)} = w_j
#         we obtain the desired implication
#             (w_j - w_i < 0) => (T[i][j] = 0).
#         The matrix T is unipotent, since unipotent matrices form a
#         normal subgroup of the general linear group and A is unipotent.
#     """
# 
#     # convert all entries to Rationals
#     for i, w in enumerate(weight_vector):
#         weight_vector[i] = QQ(w)
# 
#     permutation_matrix = _sorting_permutation_matrix(weight_vector)
#     T = permutation_matrix * A * permutation_matrix.transpose()
# 
#     return matrix(A.base_ring(), T)


def _integral_line_transformation(base_field, Vector, weight_vector):
  r"""
  Return unipotent matrix, T, over base_field transforming a line spanned
  by a standard basis vector to the line spanned by Vector such that for
  all indices i,j the implication
    (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
  holds.

  INPUT:
  base_field      - field
  Vector          - vector over field
  weight_vector   - vector over Rational Field

  MATHEMATICAL INTERPRETATION:
  First, let
    K := base_field,
    a = [a_0,...,a_n] := Vector,
    w = [w_0,...,w_n] := weight_vector.
  Let sigma be the permutation with
    w_{sigma(0)} >= ... >= w_{sigma(n)}.
  Let P be the matrix with columns given by
    e_{sigma(0)},...,e_{sigma(n)},
  where e_i is the i-th standard basis vector. Then we have
    w * P = [w_{sigma(0)},...,w_{sigma(n)}],
  i.e.
    P = _sorting_permutation_matrix(w).
  Let
    v := w * P, i.e. v_j = w_{sigma(j)}.
  Let
    b := a * P, i.e. b_j = a_{sigma(j)}.
  Let
    T = _ult_line_transformation(K, b).
  Since v is sorted in decreasing order and T is a lower
  triangular matrix, we have the implications
    (v_j - v_i < 0) => (j > i) => (T[i][j] = 0).
  Let e_j be the standard basis vector with
    e_j * T = b.
  Since
    e_{sigma^{-1}} * P = e_j,
    b * P^{-1} = a,
  we have
    e_{sigma^{-1}} * P * T * P^{-1} = a.
  Moreover, the matrix
    P^{-1} = P.transpose()
  corresponds to the permutation sigma^{-1}, i.e. the columns
  of P^{-1} are given by
    e_{sigma^{-1}(0)},...,e_{sigma^{-1}(n)}
  and therefore the (i,j)-th entry of
    T * P^{-1}
  is equal to
    T[i][sigma^{-1}(j)].
  Further, since
    P = (P^{-1})^{-1} = (P^{-1}).transpose(),
  the rows of P are given by
    e_{sigma^{-1}(0)},...,e_{sigma^{-1}(n)}.
  Thus, the (i,j) entry of
    B = P * T * P^{-1}
  is given by
    B[i][j] = T[sigma^{-1}(i)][sigma^{-1}(j)].
  Finally, since
    v_{sigma^{-1}(j)} = w_j
  we obtain the desired implication
    (w_j - w_i < 0) => (B[i][j] = 0).
  The matrix B is unipotent, since unipotent matrices form a
  normal subgroup and T is unipotent.
  """

  # convert all entries to Rationals
  for i, w in enumerate(weight_vector):
    weight_vector[i] = QQ(w)

  permutation_matrix = _sorting_permutation_matrix(weight_vector)
  T = _ult_line_transformation(base_field, vector(Vector) * permutation_matrix)

  return permutation_matrix * T * permutation_matrix.transpose()


def _ult_plane_transformation(linear_form):
  r"""
  Return a unipotent lower triangular matrix with maximal
  number of zeros which transforms `linear_form` to some `x_i`.

  EXAMPLES:
    sage: K.<A,B,C> = QQ[]
    sage: K = K.fraction_field()
    sage: R.<x0,x1,x2> = K[]
    sage: L = A*x0 + B*x1 + C*x2
    sage: T = _ult_plane_transformation(L); T
    [     1      0      0]
    [(-B)/A      1      0]
    [(-C)/A      0      1]
    sage: L(list(vector([x0,x1,x2]) * T))
    A*x0
    sage: 
    sage: L = B*x1 + C*x2
    sage: T = _ult_plane_transformation(L); T
    [     1      0      0]
    [     0      1      0]
    [     0 (-C)/B      1]
    sage: L(list(vector([x0,x1,x2]) * T))
    B*x1
    sage:
    sage: L = C*x2
    sage: T = _ult_plane_transformation(L); T
    [1 0 0]
    [0 1 0]
    [0 0 1]
  """

  base_ring = linear_form.base_ring()
  x0, x1, x2 = linear_form.parent().gens()
  A = linear_form.monomial_coefficient(x0)
  B = linear_form.monomial_coefficient(x1)
  C = linear_form.monomial_coefficient(x2)

  if A != 0:
    T = [[1, 0, 0], [-B/A, 1, 0], [-C/A, 0, 1]]
    return matrix(base_ring, T)
  elif B != 0:
    T = [[1, 0, 0], [0, 1, 0], [0, -C/B, 1]]
    return matrix(base_ring, T)
  else:
    return identity_matrix(base_ring, 3)


def _uut_plane_transformation(linear_form):
  r"""
  Return a unipotent upper triangular matrix with maximal
  number of zeros which transforms `linear_form` to some `x_i`.

  EXAMPLES:
    sage: K.<A,B,C> = QQ[]
    sage: K = K.fraction_field()
    sage: R.<x0,x1,x2> = K[]
    sage: L = A*x0 + B*x1 + C*x2
    sage: T = _uut_plane_transformation(L); T
    [     1      0 (-A)/C]
    [     0      1 (-B)/C]
    [     0      0      1]
    sage: L(list(vector([x0,x1,x2]) * T))
    C*x2
    sage:
    sage: L = A*x0 + B*x1
    sage: T = _uut_plane_transformation(L); T
    [     1 (-A)/B      0]
    [     0      1      0]
    [     0      0      1]
    sage: L(list(vector([x0,x1,x2]) * T))
    B*x1
    sage:
    sage: L = A*x0
    sage: T = _uut_plane_transformation(L); T
    [1 0 0]
    [0 1 0]
    [0 0 1]
  """

  base_ring = linear_form.base_ring()
  x0, x1, x2 = linear_form.parent().gens()
  A = linear_form.monomial_coefficient(x0)
  B = linear_form.monomial_coefficient(x1)
  C = linear_form.monomial_coefficient(x2)

  if C != 0:
    T = [[1, 0, -A/C], [0, 1, -B/C], [0, 0, 1]]
    return matrix(base_ring, T)
  elif B != 0:
    T = [[1, -A/B, 0], [0, 1, 0], [0, 0, 1]]
    return matrix(base_ring, T)
  else:
    return identity_matrix(base_ring, 3)


def _integral_plane_transformation(linear_form, weight_vector):
  r"""
  Return unipotent matrix, T, over base_field transforming a plane given
  by x_i = 0 to the plane defined by linear_form = 0 such that for all
  indices i,j the implication
    (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
  holds.

  INPUT:
  linear_form     - linear form in K[x_0, x_1, x_2]
  weight_vector   - tuple of rational numbers

  MATHEMATICAL INTERPRETATION:
  First, let
    L = linear_form
    w = [w_0,...,w_n] = weight_vector.
  Let sigma be the permutation with
    w_{sigma(0)} >= ... >= w_{sigma(n)}.
  Let P be the matrix with columns given by
    e_{sigma(0)},...,e_{sigma(n)},
  where e_i is the i-th standard basis vector. Then we have
    w * P = [w_{sigma(0)},...,w_{sigma(n)}],
  i.e.
    P = _sorting_permutation_matrix(w).
  Let
    v = w * P, i.e. v_j = w_{sigma(j)}.
  Let
    pL = L((x_0,x_1,x_2)*P^{-1}),
  i.e.
    pL = _apply_matrix(P.inverse(), L).
  Let
    T = _ult_plane_transformation(l)[0]
  and
    TpL = _apply_matrix(T, PL).
  Since v is sorted in decreasing order and T is a lower
  triangular matrix, we have the implications
    (v_j - v_i < 0) => (j > i) => (T[i][j] = 0).
  Let x_i be the generator with
    TpL = a*x_i, a in K.
  Since
    (x_0,x_1,x_2) * P = (x_{sigma(0)},x_{sigma(1)},x_{sigma(2)})
  we have
    TpL((x_0,x_1,x_2)*P) = a*x_{sigma(i)}.
  Let
    PTpL = _apply_matrix(P, TpL).
  All in all,
    PTpL = L((x_0,x_1,x_2)*P*T*P^{-1}),
  i.e.
    PTpL = _apply_matrix(P*T*P.inverse(), L).
  Note, that the matrix
    P^{-1} = P.transpose()
  corresponds to the permutation sigma^{-1}, i.e. the columns
  of P^{-1} are given by
    e_{sigma^{-1}(0)},...,e_{sigma^{-1}(n)}
  and therefore the (i,j)-th entry of
    T * P^{-1}
  is equal to
    T[i][sigma^{-1}(j)].
  Further, since
    P = (P^{-1})^{-1} = (P^{-1}).transpose(),
  the rows of P are given by
    e_{sigma^{-1}(0)},...,e_{sigma^{-1}(n)}.
  Thus, the (i,j) entry of
    B = P * T * P^{-1}
  is given by
    B[i][j] = T[sigma^{-1}(i)][sigma^{-1}(j)].
  Finally, since
    v_{sigma^{-1}(j)} = w_j
  we obtain the desired implication
    (w_j - w_i < 0) => (B[i][j] = 0).
  The matrix B is unipotent, since unipotent matrices form a
  normal subgroup and T is unipotent.
  """

  weight_vector_qq = [QQ(w) for w in weight_vector]
  P = _sorting_permutation_matrix(weight_vector_qq)
  pL = _apply_matrix(P.transpose(), linear_form)
  T = _ult_plane_transformation(pL)

  return P * T * P.transpose()


def _ult_flag_transformation(Vector, linear_form):
  r"""
  Return a unipotent lower triangular matrix transforming a flag given
  by a line spanned by a standard basis vector e_j and a plane x_i = 0
  to the line spanned by Vector and the plane given by linear_form = 0.

  INPUT:
  Vector      - vector with 3 entries
  linear_form - linear form with linear_form(Vector) = 0

  OUTPUT:
  T - unipotent lower triangular matrix with e_j*T = Vector and
  _apply_matrix(T, linear_form) = x_i

  EXAMPLES:
    sage: K.<a,b,c,A,B,C> = QQ[]
    sage: K = K.fraction_field()
    sage: R.<x0,x1,x2> = K[]
    sage: P = [a,b,c]
    sage: L = A*x0 + B*x1 - (a*A/c + b*B/c)*x2
    sage: L(P)
    0
    sage: T = _ult_flag_transformation(P, L); T
    [     1      0      0]
    [(-B)/A      1      0]
    [   a/c    b/c      1]
    sage: c * vector([0,0,1]) * T
    (a, b, c)
    sage: L(list(vector([x0,x1,x2]) * T))
    A*x0
    sage:
    sage: P = [a,b,0]
    sage: L = A*x0 - (a*A/b)*x1 + C*x2
    sage: L(P)
    0
    sage: T = _ult_flag_transformation(P, L); T
    [     1      0      0]
    [   a/b      1      0]
    [(-C)/A      0      1]
    sage: b * vector([0,1,0]) * T
    (a, b, 0)
    sage: L(list(vector([x0,x1,x2]) * T))
    A*x0
    sage:
    sage: P = [a,0,0]
    sage: L = B*x1 + C*x2
    sage: L(P)
    0
    sage: T = _ult_flag_transformation(P, L); T
    [     1      0      0]
    [     0      1      0]
    [     0 (-C)/B      1]
    sage: a * vector([1,0,0]) * T
    (a, 0, 0)
    sage: L(list(vector([x0,x1,x2]) * T))
    B*x1
  """

  Vector = list(Vector)
  if linear_form(Vector) != 0:
    raise ValueError(f"{linear_form} must be zero at {Vector}")

  base_ring = linear_form.base_ring()
  T1 = _ult_line_transformation(base_ring, Vector)
  T2 = _ult_plane_transformation(_apply_matrix(T1, linear_form))

  return T2 * T1


def _uut_flag_transformation(Vector, linear_form):
  r"""
  Return a unipotent upper triangular matrix transforming a flag given
  by a line spanned by a standard basis vector e_j and a plane x_i = 0
  to the line spanned by Vector and the plane given by linear_form = 0.

  INPUT:
  Vector      - vector with 3 entries
  linear_form - linear form with linear_form(Vector) = 0

  OUTPUT:
  T - unipotent upper triangular matrix with e_j*T = Vector and
  _apply_matrix(T, linear_form) = x_i

  EXAMPLES:
    sage: K.<a,b,c,A,B,C> = QQ[]
    sage: K = K.fraction_field()
    sage: R.<x0,x1,x2> = K[]
    sage: P = [a,b,c]
    sage: L = -(b*B/a + c*C/a)*x0 + B*x1 + C*x2
    sage: L(P)
    0
    sage: T = _uut_flag_transformation(P, L); T
    [     1    b/a    c/a]
    [     0      1 (-B)/C]
    [     0      0      1]
    sage: a * vector([1,0,0]) * T
    (a, b, c)
    sage: L(list(vector([x0,x1,x2]) * T))
    C*x2
    sage:
    sage: P = [0,b,c]
    sage: L = A*x0 - (c*C/b)*x1 + C*x2
    sage: L(P)
    0
    sage: T = _uut_flag_transformation(P, L); T
    [     1      0 (-A)/C]
    [     0      1    c/b]
    [     0      0      1]
    sage: b * vector([0,1,0]) * T
    (0, b, c)
    sage: L(list(vector([x0,x1,x2]) * T))
    C*x2
    sage:
    sage: P = [0,0,c]
    sage: L = A*x0 + B*x1
    sage: L(P)
    0
    sage: T = _uut_flag_transformation(P, L); T
    [     1 (-A)/B      0]
    [     0      1      0]
    [     0      0      1]
    sage: c * vector([0,0,1]) * T
    (0, 0, c)
    sage: L(list(vector([x0,x1,x2]) * T))
    B*x1
  """

  Vector = list(Vector)
  if linear_form(Vector) != 0:
    raise ValueError

  base_field = linear_form.base_ring()
  T1 = _uut_line_transformation(base_field, Vector)
  T2 = _uut_plane_transformation(_apply_matrix(T1, linear_form))

  return T2 * T1


def _integral_flag_transformation(Vector, linear_form, weight_vector):
  r"""
  Return unipotent matrix, T, transforming a flag given by a line spanned
  by a standard basis vector e_j and a plane x_i = 0 to the line spanned
  by Vector and the plane given by linear_form = 0 such that for all
  indices i,j the implication
    (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
  holds.

  INPUT:
  Vector        - vector with 3 entries
  linear_form   - linear form with linear_form(Vector) = 0
  weight_vector - tuple of rational numbers

  OUTPUT:
  T - unipotent lower triangular matrix with e_j*T = Vector and
  _apply_matrix(T, linear_form) = x_i

  MATHEMATICAL INTERPRETATION:
  ...
  """

  Vector = list(Vector)
  if linear_form(Vector) != 0:
    raise ValueError

  base_field = linear_form.base_ring()
  T1 = _integral_line_transformation(base_field, Vector, weight_vector)
  T2 = _integral_plane_transformation(_apply_matrix(T1.inverse(), linear_form), weight_vector)[0]

  return T2 * T1
