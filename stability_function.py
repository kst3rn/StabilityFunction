
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.all import *
from linear_valuations import LinearValuation
from transformations import _apply_matrix


class StabilityFunction:
  def __init__( self, homogeneous_form, base_ring_valuation ):
    r"""
    INPUT:
    homogeneous_form - homogeneous polynomial in K[x_0,...,x_n]
    base_ring_valuation - discrete valuation on K
    """

    if homogeneous_form.base_ring() != base_ring_valuation.domain():
      raise ValueError(f"Base rings of {homogeneous_form} and {base_ring_valuation} are not equal")

    self._homogeneous_form    = homogeneous_form
    self._base_ring_valuation = base_ring_valuation

    self._standard_basis = homogeneous_form.parent().gens()
    self._base_ring      = base_ring_valuation.domain()
    self._dimension      = Integer(len(self._standard_basis) - 1)


  def __repr__(self):
    return f"Stability Function of {self.homogeneous_form()} over {self.base_ring()} with {self.base_ring_valuation()}"


  def homogeneous_form(self):
    return self._homogeneous_form


  def base_ring_valuation(self):
    return self._base_ring_valuation


  def dimension(self):
    return self._dimension


  def standard_basis(self):
    return self._standard_basis


  def polynomial_ring(self):
    return self.homogeneous_form().parent()


  def base_ring(self):
    return self._base_ring


  def graded_reduction(self, point_on_BTB):
    linear_valuation = LinearValuation(self.polynomial_ring(),
                                       self.base_ring_valuation(),
                                       point_on_BTB.base_change_matrix(),
                                       point_on_BTB.weight_vector())
    return linear_valuation.graded_reduction(self.homogeneous_form())


  def has_semistable_reduction_at(self, point_on_BTB):
    return self.graded_reduction(point_on_BTB).is_graded_semistable()


  def has_stable_reduction_at(self, point_on_BTB):
    return self.graded_reduction(point_on_BTB).is_graded_stable()


  def initial_form(self, point_on_BTB):
    return point_on_BTB.linear_valuation().initial_form(self.homogeneous_form())


  def active_functions_at(self, base_change_matrix, weight_vector):
    return RestrictedStabilityFunction(self, base_change_matrix).active_functions(weight_vector)


  def descent_direction(self, point_on_BTB, matrix_form = 'ult'):
    r"""
    Return a matrices which describes a base change,
    fixing `point_on_BTB`, to an apartment, where `self`
    can be minimized further.

    INPUT:
    - ``point_on_BTB`` -- point in the domain of `self`.
    - ``matrix_form``  -- one of the strings 'ult', 'uut', 'integral'.

    OUTPUT:
    An invertible matrix.

    MATHEMATICAL INTERPRETATION:
    First, let
      A   = point_on_BTB.base_change_matrix(),
      u   = point_on_BTB.weight_vector(),
      B   = A.inverse(),
      K   = self.base_ring()
      n   = self.dimension()
      E_0 = self.standard_basis() .
    Then the matrix A lies in GL_n(K) and point_on_BTB is represented by the valuation v_{E,u}.
    Now we view E_0 = (x_0,...,x_n) as a vector in Sage and define the basis
      E_1 = (x_0,...,x_n)*B .
    Furthermore, let phi be the stability function, `self`. We want to find a matrix T in GL_n(K)
    such that for the basis E_2 = E_1*T^(-1) we get
      phi(point_on_BTB) < max phi_{E_2},
    where phi_{E_2} is the stability function phi restricted to the apartment given by E_2.
    Note that if n = 2, we can compute E_2 by computing a graded instability (E_2, w) of the graded
    reduction of 'self.homogeneous_form()' with respect to a valuation representing 'point_on_BTB'.
    """

    if self.dimension() != 2:
      raise NotImplementedError
    if self.base_ring() != point_on_BTB.base_ring():
      raise ValueError

    graded_reduction = self.graded_reduction(point_on_BTB)
    rational_graded_instability = graded_reduction.rational_graded_instability(matrix_form)
    if rational_graded_instability is None:
      return None
    return rational_graded_instability.lift_matrix()


  def local_minimum(self, base_change_matrix):
    r"""
    Return the minimum on the apartment defined by
    `base_change_matrix` and the point where `self`
    attains it.

    INPUT:
    - ``base_change_matrix`` -- an invertible matrix.

    EXAMPLES::
      sage: R.<x0,x1,x2> = QQ[]
      sage: F = x2*x1^2 - x0^3 - x0*x2^2
      sage: phi = StabilityFunction(F, QQ.valuation(2))
      sage: T = matrix(QQ, [[1,0,0],[1,1,0],[1,0,1]]); T
      [1 0 0]
      [1 1 0]
      [1 0 1]
      sage: a,b = phi.local_minimum(T)
      sage: a
      -1/6
      sage: b.weight_vector()
      [1/2, 1/3, 0]
      sage: b.base_change_matrix() == T
      True
    """
    return ApartmentStabilityFunction(self, base_change_matrix).minimize()


  def global_minimum(self, matrix_form = 'ult'):
    r"""
    Return the minimum on the Bruhat-Tit building and
    the point where `self` attains it.

    INPUT:
    - ``base_change_matrix`` -- an invertible matrix.
    - ``matrix_form`` -- one of the strings 'ult', 'uut', 'integral'.

    EXAMPLES::
      sage: R.<x0,x1,x2> = QQ[]
      sage: F = x0*x2*(x1^2 + 6*x0*x2)
      sage: phi = StabilityFunction(F, QQ.valuation(2))
      sage: a, b = phi.global_minimum()
      sage: a
      -1/3
      sage: b.weight_vector()
      [0, 1/2, 0]
      sage: b.base_change_matrix()
      [1 0 0]
      [0 1 0]
      [0 0 1]
      sage: b.minimal_simplex_dimension()
      1

      sage: R.<x0,x1,x2> = QQ[]
      sage: F = x2*x1^2 - x0^3 - x0*x2^2
      sage: phi = StabilityFunction(F, QQ.valuation(2))
      sage: a, b = phi.global_minimum()
      sage: a
      -1/6
      sage: b.weight_vector()
      [1/2, 1/3, 0]
      sage: b.base_change_matrix()
      [1 0 0]
      [1 1 0]
      [1 0 1]
      sage: a, b = phi.global_minimum('uut')
      sage: a
      -1/6
      sage: b.weight_vector()
      [0, 1/3, 1/2]
      sage: b.base_change_matrix()
      [1 0 1]
      [0 1 1]
      [0 0 1]
      sage: b.minimal_simplex_dimension()
      2
    """

    global_trafo_matrix = identity_matrix(self.base_ring(), self.dimension() + 1)

    while True:
      # 1. Find minimum on the current apartment defined by global_trafo_matrix.
      min_val, btb_point = self.local_minimum(global_trafo_matrix)

      # 2. If btb_point is in the interior of a chamber,
      #    we are at the global minimum.
      if btb_point.minimal_simplex_dimension() == self.dimension():
        return [min_val, btb_point]

      # 3. If there is no descent direction at btb_point,
      #    we are at the global minimum.
      local_trafo_matrix = self.descent_direction(btb_point, matrix_form)
      if local_trafo_matrix is None:
        return [min_val, btb_point]

      # 4. Update the base change matrix for the next iteration.
      global_trafo_matrix = local_trafo_matrix * global_trafo_matrix


  def evaluate_at(self, point_on_BTB):
    r"""
    Return the value of self at 'point_on_BTB'
    """
    T = point_on_BTB.base_change_matrix()
    w = point_on_BTB.weight_vector()
    return min(affine_function(w) for affine_function in self.affine_functions_on_apartment(T))


  def show(self, base_change_matrix, affine_patch = None):
    return ApartmentStabilityFunction(self, base_change_matrix).show(affine_patch)



class ApartmentStabilityFunction:
  r"""
  Construct the restriction of a stability function to an apartment
  to the following conditions.

  INPUT:
  - ``stability_function`` -- a stability function.
  - ``base_change_matrix`` -- an invertible matrix.
  """

  def __init__(self, stability_function, base_change_matrix):
    r"""
    Construct the restriction of a stability function to an apartment
    to the following conditions.

    INPUT:
    - ``stability_function`` -- a stability function.
    - ``base_change_matrix`` -- an invertible matrix.

    OUTPUT:
    Restriction of `stability_function` to the apartment defined by
    the basis
      (x_0,...,x_n) * T^{-1},
    where
      T = base_change_matrix,
      (x_0,...,x_n) = stability_function.standard_basis().
    """

    if not base_change_matrix.is_invertible():
      raise ValueError

    self._stability_function = stability_function
    self._homogeneous_form = stability_function.homogeneous_form()
    self._base_ring_valuation = stability_function.base_ring_valuation()
    self._base_change_matrix = base_change_matrix


  def __repr__(self):
    return f"Restriction of {self.stability_function()} to the apartment given by self.base_change_matrix()"


  def show(self, affine_patch=None, arg_name='w'):
    r"""
    EXAMPLES::
      sage: R.<x0,x1,x2> = QQ[]
      sage: F = x0*x2*(x1^2 + x0*x2)
      sage: v_2 = QQ.valuation(2)
      sage: phi = StabilityFunction(F, v_2)
      sage: E = identity_matrix(QQ, 3)
      sage: phiE = ApartmentStabilityFunction(phi, E)
      sage: phiE.show()
      (w0, w1, w2) |--> max(1/3*w0 - 2/3*w1 + 1/3*w2, -2/3*w0 + 4/3*w1 - 2/3*w2)
      sage: phiE.show(affine_patch=0)
      (0, w1, w2) |--> max(-2/3*w1 + 1/3*w2, 4/3*w1 - 2/3*w2)
      sage: phiE.show(affine_patch=2)
      (w0, w1, 0) |--> max(1/3*w0 - 2/3*w1, -2/3*w0 + 4/3*w1)
      sage:
      sage: F = x0^2 + x1*x2
      sage: phi = StabilityFunction(F, v_2)
      sage: T = matrix(QQ, [[1,0,0],[0,1,0],[0,2,1]]); T
      [1 0 0]
      [0 1 0]
      [0 2 1]
      sage: phiT = ApartmentStabilityFunction(phi, T)
      sage: phiT.show(affine_patch=0)
      (0, w1, w2) |--> max(2/3*w1 - 4/3*w2 - 1, -1/3*w1 - 1/3*w2, 2/3*w1 + 2/3*w2)
      sage:
      sage: T = matrix(QQ, [[2,0,0],[0,2,0],[0,0,1]]); T
      [2 0 0]
      [0 2 0]
      [0 0 1]
      sage: phiT = ApartmentStabilityFunction(phi, T)
      sage: phiT.show(affine_patch=0)
      (0, w1, w2) |--> max(-1/3*w1 - 1/3*w2 + 1/3, 2/3*w1 + 2/3*w2 - 2/3)
      sage:
      sage: F = x0^2 + 6*x1*x2
      sage: phi = StabilityFunction(F, v_2)
      sage: E = identity_matrix(QQ, 3)
      sage: phiE = ApartmentStabilityFunction(phi, E)
      sage: phiE.show(affine_patch=0)
      (0, w1, w2) |--> max(-1/3*w1 - 1/3*w2 - 1, 2/3*w1 + 2/3*w2)
    """

    N = self.dimension() + 1
    if affine_patch is not None:
      if not isinstance(affine_patch, (int, Integer)):
        raise ValueError(f"{affine_patch} is not an integer")
      elif not 0 <= affine_patch < N:
        raise ValueError(f"{affine_patch} is not between {0} and {N}")
    if not isinstance(arg_name, str):
      raise ValueError(f"{arg_name} is not a string")
    elif not len(arg_name) == 1:
      raise ValueError(f"{arg_name} is not a character")

    R = PolynomialRing(QQ, N, arg_name)
    w = R.gens()
    if affine_patch is not None:
      w = list(w)
      w[affine_patch] = 0
      w = tuple(w)
    aff_forms = []
    for const, lin_form in self.affine_forms(redundancy=False):
      aff_forms.append(const + sum(lin_form[j] * w[j] for j in range(N)))
    print(str(w) + " |--> max" + str(tuple(aff_forms)))


  def stability_function(self):
    return self._stability_function


  def homogeneous_form(self):
    return self._homogeneous_form


  def base_ring_valuation(self):
    return self._base_ring_valuation


  def base_change_matrix(self):
    return self._base_change_matrix


  def dimension(self):
    return self.stability_function().dimension()


  # def active_functions(self, w, flag = True):
  #   r"""
  #   Return the set of active functions a w
  #   """
  # 
  #   d = self.homogeneous_form().degree()
  #   N = self.dimension() + 1
  #   # Compute d/N*v_K( det(A) )
  #   const_A = d/N*self.base_ring_valuation(self._embedding_matrix.det())
  #   affine_functions_values = dict()
  #   F = _apply_matrix(self._embedding_matrix, self.homogeneous_form())
  #   for multi_index, coefficient in F.dict().items():
  #     value_at_w = self.base_ring_valuation()(coefficient) - const_A
  #     for j in range(N):
  #       value_at_w = value_at_w + multi_index[j] * w[j]
  #     affine_functions_values[multi_index] = value_at_w
  #   min_value = min(affine_functions_values.values())
  #   return [key for key, value in affine_functions_values.items() if value == min_value]


  def affine_forms(self, redundancy=True):
    r"""
    Return the affine forms defining `self`.

    OUTPUT:
    {(c_1, (a_01,...,a_n1), (c_2, (a_02,...,a_n2),...}

    EXAMPLES::
      sage: R.<x0,x1,x2> = QQ[]
      sage: v_2 = QQ.valuation(2)
      sage: F = x0*x2*(x1^2 + x0*x2)
      sage: phi = StabilityFunction(F, v_2)
      sage: E = identity_matrix(QQ, 3)
      sage: phiE = ApartmentStabilityFunction(phi, E)
      sage: phiE.affine_forms()
      [(0, (1/3, -2/3, 1/3)), (0, (-2/3, 4/3, -2/3))]
      sage:
      sage: T = matrix(QQ, [[1,0,0],[2,1,0],[3,0,1]]); T
      [1 0 0]
      [2 1 0]
      [3 0 1]
      sage: phiT = ApartmentStabilityFunction(phi, T)
      sage: phiT.affine_forms()
      [(0, (1/3, -2/3, 1/3)),
      (-1, (4/3, -5/3, 1/3)),
      (0, (-2/3, 4/3, -2/3)),
      (-2, (1/3, 1/3, -2/3)),
      (0, (4/3, -2/3, -2/3)),
      (-1, (1/3, 4/3, -5/3)),
      (-2, (4/3, 1/3, -5/3)),
      (0, (4/3, 4/3, -8/3))]
      sage: phiT.affine_forms(redundancy=False)
      [(0, (4/3, 4/3, -8/3)),
      (0, (-2/3, 4/3, -2/3)),
      (0, (4/3, -2/3, -2/3)),
      (-1, (4/3, -5/3, 1/3)),
      (0, (1/3, -2/3, 1/3))]
      sage:
      sage: T = matrix(QQ, [[1,0,0],[2,2,0],[3,0,1]]); T
      [1 0 0]
      [2 2 0]
      [3 0 1]
      sage: phiT = ApartmentStabilityFunction(phi, T)
      sage: phiT.affine_forms()
      [(-2/3, (1/3, -2/3, 1/3)),
      (-5/3, (4/3, -5/3, 1/3)),
      (4/3, (-2/3, 4/3, -2/3)),
      (-2/3, (1/3, 1/3, -2/3)),
      (-8/3, (4/3, -2/3, -2/3)),
      (1/3, (1/3, 4/3, -5/3)),
      (-2/3, (4/3, 1/3, -5/3)),
      (4/3, (4/3, 4/3, -8/3))]
      sage: phiT.affine_forms(redundancy=False)
      [(4/3, (4/3, 4/3, -8/3)),
      (4/3, (-2/3, 4/3, -2/3)),
      (-5/3, (4/3, -5/3, 1/3)),
      (-2/3, (1/3, -2/3, 1/3))]

    .. MATH::
    First, let
      v_K = self.base_ring_valuation(),
      A   = self.base_change_matrix(),
      B   = A.inverse(),
      E_0 = (x_0,...,x_n) = self.standard_basis(),
      F   = self.homogeneous_form().
    Thus, F is a homogeneous polynomial in K[x_0,...,x_n]. Further, we call E_0 the standard
    basis and consider A and B as linear transformations with respect to E_0, i.e.
      A(x_j) = sum_{i=0}^n a_{ij}*x_i  and  B(x_j) = sum_{i=0}^n b_{ij}*x_i.
    Then,
      E := (y_0,...,y_n) := (B(x_0),...,B(x_n))
    is a new basis of K[x_0,...,x_n]. Now if we view E_0 = (x_0,...,x_n) as a vector in Sage,
    we get
      (y_0,...,y_n) = (sum_{i=0}^n b_{i,0}*x_i,...,sum_{i=0}^n b_{i,n}*x_i)
                    = (x_0,...,x_n)*B
    and therefore
      F(x_0,...,x_n) = F((y_0,...,y_n)*B^{-1}) = F((y_0,...,y_n)*A).
    Thus, the homogeneous polynomial
      G(y_0,...,y_n) := F((y_0,...,y_n)*A) in K[y_0,...,y_n]
    describes F with respect to the basis E = (y_0,...,y_n) and A describes the base change.
    Thus, for the valuation v_{E,w} we obtain
      v_{E,w}(F) = min(v_K(a_i) + <i,w> : i in I) with G = sum_{i in I} a_i y^i,
    where i is a multi-index, i.e. I is a subset of NN^{n+1}. Moreover, we have
      omega(v_{E,w}) = 1/(n+1) * (w_0 + ... + w_n - v_K(det(E))).
    Note that per definition det(E) = det(B). Furthermore,
      v_K(det(B)) = v_K(det(A^{-1})) = v_K(det(A)^{-1}) = -v_K(det(A))
    and therefore
      omega(v_{E,w}) = 1/(n+1) * (w_0 + ... + w_n + v_K(det(A))).
    Now let N = n + 1. Then we have
      phi_E(w) = d*omega(v_{E,w}) - v_{E,w}(F)
               = max(d/N*v_K(det(A)) - v_K(a_i) + sum_{j=0}^n (d/N - i_j)*w_j : i in I).
    """

    # Set up variables
    d = Integer(self.homogeneous_form().degree())
    N = self.dimension() + 1   # N = n + 1
    v_K = self.base_ring_valuation()
    val_det_A = d / N * v_K(self.base_change_matrix().det())

    # Compute d*omega(v_{E,w}) - v_{E,w}(F)
    G = _apply_matrix(self.base_change_matrix(), self.homogeneous_form())
    aff_forms = []
    for multi_index, G_coeff in G.dict().items():
      aff_forms.append((val_det_A - v_K(G_coeff),
                        tuple(d / N - i_j for i_j in multi_index)))

    if redundancy:
      return aff_forms

    inequalities_input = []
    for const, lin_form in aff_forms:
      inequalities_input.append([-const] + [-c for c in lin_form] + [1])
    P = Polyhedron(ieqs=inequalities_input, base_ring=QQ)

    minimized_forms = []
    for row in P.inequalities_list():
      const_term = row[0]
      w_coeffs = row[1:-1]
      z_coeff = row[-1]
      scale = -1 / z_coeff
      new_const = const_term * scale
      new_lin_form = tuple(c * scale for c in w_coeffs)
      minimized_forms.append((new_const, new_lin_form))

    return minimized_forms


  def minimize(self):
    r"""
    Return the minimum of `self` and the
    point where it is attained.

    OUTPUT:
    A pair `(a, b)` consisting of a rational number `a`,
    which equals the minimum of `self` and a point `b` on
    the Bruhat-Tits building where `self` attains `a`.

    EXAMPLES::
      sage: R.<x0,x1,x2> = QQ[]
      sage: F = x0*x2*(x1^2 + 6*x0*x2)
      sage: phi = StabilityFunction(F, QQ.valuation(2))
      sage: T = matrix(QQ, [[1,0,0],[2,2,0],[3,0,1]]); T
      [1 0 0]
      [2 2 0]
      [3 0 1]
      sage: phiT = ApartmentStabilityFunction(phi, T)
      sage: a, b = phiT.minimize()
      sage: a
      -1/3
      sage: b
      Point on the Bruhat-Tits Building of SL(3) over Rational Field with 2-adic valuation
      sage: b.weight_vector()
      [1/2, 0, 1/2]
    """

    MILP = MixedIntegerLinearProgram(solver='PPL', maximization=False)
    v = MILP.new_variable()
    t = v['maximum_of_affine_forms']    
    MILP.set_objective(t)
    MILP.add_constraint(v[0] == 0)

    N = self.dimension() + 1
    for const_coeff, lin_form in self.affine_forms():
      MILP_term = const_coeff + sum(lin_form[j] * v[j] for j in range(N))
      MILP.add_constraint(MILP_term <= t)
    MILP.solve()
    solution_dict = MILP.get_values(v)
    minimum = solution_dict['maximum_of_affine_forms']
    weight_vector = [solution_dict[i] for i in range(N)]
    b = BTB_Point(self.base_ring_valuation(),
                  self.base_change_matrix(),
                  weight_vector)

    return (minimum, b)


  def optimal_polyhedron(self, affine_patch=0):
    r"""
    Return the polyhedron consistion of optimal
    solution of `self`.

    .. MATH::
    By the special structure of `self` as a piecewise affine
    affine
      RR^{n+1} ---> RR
    it factors over RR^{n+1}/RR(1,...,1). In particular, the
    maximum of `self` is the same on any affine patch.
    """

    N = self.dimension() + Integer(1)
    if affine_patch is not None:
      if not isinstance(affine_patch, (int, Integer)):
        raise ValueError(f"{affine_patch} is not an integer")
      elif not 0 <= affine_patch < N:
        raise ValueError(f"{affine_patch} is not between {0} and {N}")

    # 1. Find the maximum value.
    maximum, _ = self.maximize()

    # 2. Construct the inequalities for the optimal face.
    # The optimal face is the set of points w where phi(w) = maximum.
    # Since phi(w) = min( c_i + L_i(w) ), this is equivalent to
    # c_i + L_i(w) >= maximum   for all i
    #
    # Polyhedron representation needs inequalities in the form
    # b + a_0*w_0 + a_1*w_1 + ... >= 0
    #
    # Our inequalities are:
    # (c_i - maximum) + L_i[0]*w_0 + L_i[1]*w_1 + ... >= 0
    ieqs = []
    for const, lin_form in self.affine_forms():
      # b = const - maximum
      # [a_0, a_1, ...] = lin_form
      ieqs.append([const - maximum] + list(lin_form))

    # 3. Add the equality constraint.
    #
    # Polyhedron representation needs equalities in the form
    # b + a_0*w_0 + a_1*w_1 + ... == 0
    #
    # Our equality is:
    # 0 + 0*w_0 + ... + 1*w_{affine_patch} + ... + 0*w_{N-1} == 0
    eqns = [[Integer(0)] * (N + 1)]
    if affine_patch is not None:
      eqns[0][affine_patch + 1] = Integer(1)

    # 4. Return the polyhedron
    return Polyhedron(ieqs=ieqs, eqns=eqns, base_ring=QQ)


  def bounded_suplevel_sets(self):
    r"""
    Return `True` if the suplevel sets are bounded
    and `False` otherwise.

    .. MATH::
    The suplevel sets are bounded if and only if the
    set of optimal solutions is bounded.
    """

    return self.optimal_polyhedron().is_compact()


  def integral_points(self):
    r"""
    Return the list of all integral points
    in the level set of `self.maximize()`.
    """

    if not self.bounded_suplevel_sets():
      raise ValueError(f"Suplevel sets are unbounded")

    return self.optimal_polyhedron().integral_points()


  def semistable_models(self):
    r"""
    Return the generator yielding all semistable models
    on the apartment given by `self.base_change_matrix()`.
    """

    a, b = self.maximize()
    if not self.stability_function().has_semistable_reduction_at(b):
      return None

    v_K = self.base_ring_valuation()
    T = self.base_change_matrix()
    F = self.homogeneous_form()
    for weight_vector in self.integral_points():
      yield BTB_Point(v_K, T, weight_vector).hypersurface_model(F)



class BTB_Point:
  def __init__(self, base_ring_valuation, base_change_matrix, weight_vector):
    self._base_ring_valuation = base_ring_valuation
    self._base_change_matrix = base_change_matrix

    # convert all entries to rationals
    weight_vector_qq = [QQ(w) for w in weight_vector]
    normalized_weight_vector = []
    w_min = min(weight_vector_qq)
    for w in weight_vector_qq:
      normalized_weight_vector.append(w - w_min)
    self._weight_vector = normalized_weight_vector


  def __repr__(self):
    return f"Point on the Bruhat-Tits Building of SL({len(self._weight_vector)}) over {self._base_ring_valuation.domain()} with {self._base_ring_valuation}"


  def base_ring_valuation(self):
    return self._base_ring_valuation


  def base_change_matrix(self):
    return self._base_change_matrix


  def weight_vector(self):
    return self._weight_vector


  def affine_patch(self):
    return self._weight_vector.index(0)


  def base_ring(self):
    return self._base_change_matrix.base_ring()


  def linear_valuation(self):
    r"""
    Return a linear valuation representing self
    """
    K = self.base_ring_valuation().domain()
    N = len(self.weight_vector())
    R = PolynomialRing(K, N, 'x')
    return LinearValuation(R, self.base_ring_valuation(), self.base_change_matrix(), self.weight_vector())


  def minimal_simplex_dimension(self, ramification_index=None):
    r"""
    Return the dimension of the minimal simplex containing `self`.

    EXAMPLES::
      sage: w = [1/2, 1/3, 5/3]
      sage: E = identity_matrix(QQ, 3)
      sage: P = BTB_Point(QQ.valuation(2), E, w); P
      Point on the Bruhat-Tits Building of SL(3) over Rational Field with 2-adic valuation
      sage: P.minimal_simplex_dimension()
      2
      sage: P.minimal_simplex_dimension(ramification_index=2)
      2
      sage: P.minimal_simplex_dimension(ramification_index=3)
      1
      sage: P.minimal_simplex_dimension(ramification_index=6)
      0
      sage: w = [0, 1/2, 3/2]
      sage: P = BTB_Point(QQ.valuation(2), E, w)
      sage: P.minimal_simplex_dimension()
      1
      sage: P.minimal_simplex_dimension(ramification_index=2)
      0
      sage: w = [0, 1, 3]
      sage: P = BTB_Point(QQ.valuation(2), E, w)
      sage: P.minimal_simplex_dimension()
      0

    MATHEMATICAL INTERPRETATION (ToDo. But at this point we just give an example.): 
    If we start with the point (0, 1/2, 3/2, 1/3), we first move it to
    (0, 1/2, 1/2, 1/3). Then we see that {(0,0,0,0), (0, 1, 1, 0), (0, 1, 1, 1)}
    is the minimal simplex containing (0, 1/2, 1/2, 1/3). This is because we have
    to consider the lines <(0, 1, 0, 0)>, <(0, 0, 1, 0)>, <(0, 0, 0, 1)>,
    <(0, 1, 1, 0)>, <(0, 1, 0, 1)>, <(0, 0, 1, 1)>, <(0, 1, 1, 1)> and to find the
    minimal direct sum, containing (0, 1/2, 1/2, 1/3). In our case this is the sum
    <(0, 1, 1, 0)> oplus <(0, 0, 0, 1> = <(0,1,1,0), (0, 0, 0, 1)>.
    Since the simplices need to have a cyclic basis, we must choose
    (0, 1, 1, 0), (0, 1, 1, 1) as a basis of the direct sum above.
    In general, we first have to normalize the valuation to have the value group ZZ,
    i.e. to divide the vector 'self.weight_vector()' by v_K(pi_K). Let w_normalized be
    this normalized vector. Now we have to move w_normalized inside the cube [0,1]^N.
    Let w be this translated vector. Now remove all zeros from w. Let w_without_zeros be
    this modified vector. Now compute the cardinality of the set set(w_without_zeros).
    This cardinality is exactly the dimension of the minimal simplex containing self.
    Note that this is the same as the cardinality of the set set(w) minus 1, i.e.
    len(set(w)) - 1.
    """

    # normalize the value group to be ZZ
    if ramification_index is None:
      value_groug_generator = self.base_ring_valuation().value_group().gen()
    else:
      value_groug_generator = ZZ(1) / ramification_index
    norm_weight_vector = []
    for c in self._weight_vector:
      norm_weight_vector.append(c / value_groug_generator)

    # translate inside the unit cube
    trans_norm_weight_vector = []
    for c in norm_weight_vector:
      trans_norm_weight_vector.append(c - floor(c))

    return len(set(trans_norm_weight_vector)) - 1


  def is_vertex(self):
    return self.minimal_simplex_dimension() == 0


  def ramification_index(self):
    r"""
    Return minimal positive integer `r` such that all coordinates
    of `self.weight_vector()` are contained in 1/r * ZZ.

    EXAMPLES::
      sage: E = identity_matrix(QQ, 3)
      sage: w = [1/2, 1/3, 5/3]
      sage: P = BTB_Point(QQ.valuation(2), E, w)
      sage: P.ramification_index()
      6
      sage:
      sage: w = [1/2, 1/3, 5/11]
      sage: P = BTB_Point(QQ.valuation(2), E, w)
      sage: P.ramification_index()
      66
    """
    return lcm(b.denominator() for b in self.weight_vector())    


  def move_to_origin(self):
    r"""
    Return ...

    EXAMPLES::
      sage: w = [0, 1, 2]
      sage: E = identity_matrix(QQ, 3)
      sage: b = BTB_Point(QQ.valuation(2), E, w)
      sage: b_new = b.move_to_origin()
      sage: b_new.weight_vector()
      [0, 0, 0]
      sage: b_new.base_change_matrix()
      [1 0 0]
      [0 2 0]
      [0 0 4]
      sage:
      sage: w = [0, 1, 3]
      sage: T = matrix(QQ, [[1,0,0],[2,1,0],[5,0,1]]); T
      [1 0 0]
      [2 1 0]
      [5 0 1]
      sage: b = BTB_Point(QQ.valuation(2), T, w)
      sage: b_new = b.move_to_origin()
      sage: b_new.weight_vector()
      [0, 0, 0]
      sage: b_new.base_change_matrix()
      [ 1  0  0]
      [ 4  2  0]
      [40  0  8]
    ..MATH::
    """
    if not self.is_vertex():
      raise ValueError(f"self is not a vertex")

    g = self.base_ring_valuation().value_group().gen()
    normalized_weight_vector = [a / g for a in self.weight_vector()]
    pi_K = self.base_ring_valuation().uniformizer()
    pi_K_powers = [pi_K**a for a in normalized_weight_vector]

    D = diagonal_matrix(self.base_ring(), pi_K_powers)
    return BTB_Point(self.base_ring_valuation(),
                     D * self.base_change_matrix(),
                     [QQ(0)] * len(self.weight_vector()))


  def hypersurface_model(self, F):
    r"""
    Retrun the hypersurface model of `F` corresponding to `self`.
    """
    if self.minimal_simplex_dimension() != 0:
      raise TypeError(f"{self} is not a vertex")

    b = self.move_to_origin()
    T = b.base_change_matrix()
    v = b.linear_valuation()
    pi_K = b.base_ring_valuation().uniformizer()
    G = _apply_matrix(T, F)
    normalized_exponent = v(F) / v(pi_K)
    return G / pi_K**normalized_exponent

