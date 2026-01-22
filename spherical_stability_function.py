from sage.all import * 
import math
import numpy as np
import matplotlib.pyplot as plt

class SphericalStabilityFunction:
  r"""
  Construct a concave piecewise affine function to the following conditions

  INPUT:
  homogeneous_form - homogeneous polynomial in K[x_0,...,x_n]
  """

  def __init__(self, homogeneous_form):
    r"""
    Construct a concave piecewise affine function to the following conditions

    INPUT:
    homogeneous_form - homogeneous polynomial in K[x_0,...,x_n]
    """
    if not homogeneous_form.is_homogeneous():
      raise ValueError

    self._homogeneous_form = homogeneous_form
    self._dimension = homogeneous_form.parent().ngens()
    self._base_ring = homogeneous_form.base_ring()


  def __call__(self, transformation_matrix, position=None):
    r"""
    Return `self.evaluate(transformation_matrix, position)` if
    position is not None. Otherwise return the symbolic presentation
    of `self` on the apartment given by transformation_matrix.

    EXAMPLES:
    sage: R.<x0,x1,x2> = GF(3)[]
    sage: f = x0^2 + x1^2 + x2^2
    sage: mu = SphericalStabilityFunction(f)
    sage: T = matrix(GF(3), [[1,0,2],[1,1,0],[0,2,1]]); T
    [1 0 2]
    [1 1 0]
    [0 2 1]
    sage: mu(T)
    min(2*w0, w0 + w1, 2*w1, w0 + w2, w1 + w2, 2*w2)/||w||
    sage: mu(T, [2,-1,-1])
    -1/3*sqrt(6)
    sage: mu(T, 2)
    -2/3*sqrt(6)*sin(2)
    sage: mu(T, 3/2*pi)
    -1/3*sqrt(6)
    """
    mu = ApartmentSphericalStabilityFunction(self.homogeneous_form(),
                                             transformation_matrix)
    if position is None:
      return mu
    return mu(position)


  def dimension(self):
    return self._dimension


  def base_ring(self):
    return self._base_ring

  def homogeneous_form(self):
    return self._homogeneous_form


  def evaluate(self, transformation_matrix, position):
    r"""
    Evaluate `self` on the apartment given by transformation_matrix
    at the point given by position.

    INPUT:
    - ``transformation_matrix`` -- invertible matrix
    - ``position`` -- either a list of rational numbers or an angle in radians

    EXAMPLES::
    sage: R.<x0,x1,x2> = GF(3)[]
    sage: f = x0^2 + x1^2 + x2^2
    sage: mu = SphericalStabilityFunction(f)
    sage: T = matrix(GF(3), [[1,0,2],[1,1,0],[0,2,1]]); T
    [1 0 2]
    [1 1 0]
    [0 2 1]
    sage: mu.evaluate(T, [1,0,-1])
    -sqrt(2)
    sage: mu.evaluate(T, [2,-1,-1])
    -1/3*sqrt(6)
    """
    if not transformation_matrix.is_invertible():
      raise ValueError

    mu = ApartmentSphericalStabilityFunction(self.homogeneous_form(),
                                             transformation_matrix)
    return mu.evaluate(position)


  def active_functions(self, transformation_matrix, position):
    r"""
    Return the set of multi indices corresponding to active functions
    of `self` on the apartment given by `transformation_matrix` at the
    point given by `position`.

    EXAMPLES::
    sage: R.<x0,x1,x2> = GF(3)[]
    sage: F = x0^2 + x1^2 + x2^2
    sage: T = matrix(GF(3), [[1,2,0], [1,1,1], [0,1,1]]); T
    [1 2 0]
    [1 1 1]
    [0 1 1]
    sage: mu = SphericalStabilityFunction(F)
    sage: mu.active_functions(T, pi/6)
    {(0, 0, 2), (0, 1, 1)}
    sage: mu.active_functions(T, [2,-1,-1])
    {(0, 0, 2), (0, 1, 1)}
    sage: mu.active_functions(T, pi/2)
    {(0, 0, 2)}
    sage: mu.active_functions(T, [1,1,-2])
    {(0, 0, 2)}
    sage: mu.active_functions(T, 5*pi/6)
    {(0, 0, 2), (1, 0, 1), (2, 0, 0)}
    sage: mu.active_functions(T, [-1,2,-1])
    {(0, 0, 2), (1, 0, 1), (2, 0, 0)}
    sage: mu.active_functions(T, 7*pi/6)
    {(2, 0, 0)}
    sage: mu.active_functions(T, [-2,1,1])
    {(2, 0, 0)}
    sage: mu.active_functions(T, 3*pi/2)
    {(2, 0, 0)}
    sage: mu.active_functions(T, [-1,-1,2])
    {(2, 0, 0)}
    sage: mu.active_functions(T, 11*pi/6)
    {(0, 1, 1)}
    sage: mu.active_functions(T, [1,-2,1])
    {(0, 1, 1)}
    """
    if not transformation_matrix.is_invertible():
      raise ValueError

    mu = ApartmentSphericalStabilityFunction(self.homogeneous_form(),
                                             transformation_matrix)
    return mu.active_functions(position)


  def stability_status(self):
    r"""
    Return 'unstable', 'strictly semistable' or 'stable', depending on
    whether the projective hypersurface defined by `self.homogeneous_form()`
    is unstable, strictly semistable or stable, respectively.

    EXAMPLES:
    sage: R.<x0,x1,x2> = GF(17^3)[]
    sage: f = x0^3*(x0 + x1 + x2)
    sage: mu = SphericalStabilityFunction(f)
    sage: mu.stability_status()
    'unstable'

    sage: R.<x0,x1,x2> = GF(2)[]
    sage: f = (x0^2 + x1*x2)^2
    sage: mu = SphericalStabilityFunction(f)
    sage: mu.stability_status()
    'strictly semistable'

    sage: R.<x0,x1,x2> = GF(3)[]
    sage: f = x0*x1*x2*(x0 + x1 + x2)
    sage: mu = SphericalStabilityFunction(f)
    sage: mu.stability_status()
    'stable'
    """

    R = self.base_ring()
    F = self.homogeneous_form()
    n = self.dimension() - 1
    ult_matrices = _unipotent_lower_triangular_matrices(R, n)

    signs = set()
    for T in ult_matrices:
      mu = ApartmentSphericalStabilityFunction(F, T)
      apartment_sign = mu.sign_of_maximum()
      if apartment_sign == 1:
        return 'unstable'
      signs.add(apartment_sign)
    if max(signs) == 0:
      return 'strictly semistable'
    return 'stable'


  def plot_cartesian(self, transformation_matrix,
                     plot_individual_Li=True,
                     base_radius=1.0):
    r"""
    Generates a 3D cylindrical plot of f(theta) = min_{i in J} L_i(theta)
    on the unit circle S_H, using the direct orthonormal parameterization:
    w(theta) = cos(theta)v1 + sin(theta)v2.

    L_i(theta) = A1 * cos(theta) + A_2 * sin(theta) where
    v1 = (1/sqrt(2))*(1,-1,0), v2 = (1/sqrt(6))*(1,1,-2)
    A_1 = <i, v1> = (i0 - i1) / sqrt(2)
    A_2 = <i, v2> = (i0 + i1 - 2*i2) / sqrt(6)

    INPUT:
    transformation_matrix - invertible matrix over the base ring of
                            self.homogeneous_form().
    """
    if self.dimension() != 3:
      raise NotImplementedError(
          "This plotting method is implemented only for dimension 3 "
          "(homogeneous polynomials in 3 variables x0, x1, x2)."
      )

    if not transformation_matrix.is_invertible():
      raise ValueError

    mu = ApartmentSphericalStabilityFunction(self.homogeneous_form(),
                                             transformation_matrix)
    mu.plot_cartesian(plot_individual_Li, base_radius)


  def plot_polar(self, transformation_matrix, plot_individual_Li=True):
    r"""
    Generates a 2D polar plot of r = f(theta) = min_{i in J} L_i(theta).
    The angle theta corresponds to the parameterization of the unit circle S_H
    w(theta) = cos(theta)v1 + sin(theta)v2, and r is the value of the function.

    L_i(theta) = A_1 * cos(theta) + A_2 * sin(theta) where
    v1 = (1/sqrt(2))*(1,-1,0), v2 = (1/sqrt(6))*(1,1,-2)
    A_1 = <i, v1> = (i0 - i1) / sqrt(2)
    A_2 = <i, v2> = (i0 + i1 - 2*i2) / sqrt(6)

    INPUT:
    transformation_matrix - invertible matrix over the base ring of
                            self.homogeneous_form().
    plot_individual_Li    - boolean, whether to plot individual L_i(theta) curves.
    """

    if self.dimension() != 3:
      raise NotImplementedError(
          "This plotting method is implemented only for dimension 3 "
          "(homogeneous polynomials in 3 variables x0, x1, x2)."
      )

    if not transformation_matrix.is_invertible():
      raise ValueError

    mu = ApartmentSphericalStabilityFunction(self.homogeneous_form(),
                                             transformation_matrix)
    mu.plot_polar(plot_individual_Li)



class ApartmentSphericalStabilityFunction:
  r"""
  Construct spherical stability function restricted to an apartment
  to the following conditions.

  INPUT:
    - ``homogeneous_form`` -- homogeneous polynomial in K[x_0,...,x_n].
    - ``transformation_matrix`` -- invertible matrix over K.

  OUTPUT:
    The restriction of the spherical stability function
    of `F` to the apartment given by `(x_0,...,x_n) * T`, where
    `F` is the homogeneous form `homogeneous_form` and `T` the
    transformation matrix `transformation_matrix`.

  EXAMPLES::
    sage: R.<x0,x1,x2> = GF(3)[]
    sage: F = x0^2 + x1^2 + x2^2
    sage: T = matrix(GF(3), [[1,0,2],[1,1,0],[0,2,1]])
    sage: mu = ApartmentSphericalStabilityFunction(F, T)
    sage: mu
    max(-2*w0, -w0 - w1, -2*w1, -w0 - w2, -w1 - w2, -2*w2)/||w||
    sage:
    sage: T = matrix(GF(3), [[1,2,0], [1,1,1], [0,1,1]])
    sage: mu = ApartmentSphericalStabilityFunction(F, T)
    sage: mu
    max(-2*w0, -w0 - w2, -w1 - w2, -2*w2)/||w||
    sage:
    sage: T = identity_matrix(GF(3), 3)
    sage: F = x0^4 + x1^4 + x1^2*x2^2
    sage: mu = ApartmentSphericalStabilityFunction(F, T)
    sage: mu
    max(-4*w0, -4*w1, -2*w1 - 2*w2)/||w||

  ..MATH::
  Let F \in K[x_0,...,x_n] be a homogeneous form.
  Let T \in GL_{n+1}(K) be an invertible matrix.
  Let (y_0,...,y_n) = (x_0,...,x_n) * T^{-1}.
  Note that to write F with respect to the basis (y_0,...,y_n) means
  to consider the homogeneous form
  G(y_0,...,y_n) = F((y_0,...,y_n) * T).
  Write
  G = \sum_{i \in J} a_i x_0^{i_0} ... x_n^{i_n}
  with
  J = {i : a_i \neq 0} \subset \NN^{n+1}.
  Let H = {w \in \RR^{n+1} : w_0 + ... + w_n = 0} \setminus \{0\}.
  The spherical stability function of `F` restricted to the apartment
  given by `(y_0,...,y_n)` is the function
  H ---> \RR, w \mapsto max_{i \in J}(-<i,w>) / ||w||,
  where ||w|| is the euclidean norm of w and
  <i,w> = i_0*w_0 + ... + i_n*w_n.
  """

  def __init__(self, homogeneous_form, transformation_matrix):
    r"""
    Construct spherical stability function restricted to an apartment
    to the following conditions.

    INPUT:
      - ``homogeneous_form`` -- homogeneous polynomial in K[x_0,...,x_n].
      - ``transformation_matrix`` -- invertible matrix over K.

    OUTPUT:
      The restriction of the spherical stability function
      of `F` to the apartment given by `(x_0,...,x_n) * T`, where
      `F` is the homogeneous form `homogeneous_form` and `T` the
      transformation matrix `transformation_matrix`.

    EXAMPLES::
      sage: R.<x0,x1,x2> = GF(3)[]
      sage: F = x0^2 + x1^2 + x2^2
      sage: T = matrix(GF(3), [[1,0,2],[1,1,0],[0,2,1]])
      sage: mu = ApartmentSphericalStabilityFunction(F, T)
      sage: mu
      max(-2*w0, -w0 - w1, -2*w1, -w0 - w2, -w1 - w2, -2*w2)/||w||
      sage:
      sage: T = matrix(GF(3), [[1,2,0], [1,1,1], [0,1,1]])
      sage: mu = ApartmentSphericalStabilityFunction(F, T)
      sage: mu
      max(-2*w0, -w0 - w2, -w1 - w2, -2*w2)/||w||
      sage:
      sage: T = identity_matrix(GF(3), 3)
      sage: F = x0^4 + x1^4 + x1^2*x2^2
      sage: mu = ApartmentSphericalStabilityFunction(F, T)
      sage: mu
      max(-4*w0, -4*w1, -2*w1 - 2*w2)/||w||

    ..MATH::
    Let F \in K[x_0,...,x_n] be a homogeneous form.
    Let T \in GL_{n+1}(K) be an invertible matrix.
    Let (y_0,...,y_n) = (x_0,...,x_n) * T^{-1}.
    Note that to write F with respect to the basis (y_0,...,y_n) means
    to consider the homogeneous form
    G(y_0,...,y_n) = F((y_0,...,y_n) * T).
    Write
    G = \sum_{i \in J} a_i x_0^{i_0} ... x_n^{i_n}
    with
    J = {i : a_i \neq 0} \subset \NN^{n+1}.
    Let H = {w \in \RR^{n+1} : w_0 + ... + w_n = 0} \setminus \{0\}.
    The spherical stability function of `F` restricted to the apartment
    given by `(y_0,...,y_n)` is the function
    H ---> \RR, w \mapsto max_{i \in J}(-<i,w>) / ||w||,
    where ||w|| is the euclidean norm of w and
    <i,w> = i_0*w_0 + ... + i_n*w_n.
    """
    if not transformation_matrix.is_invertible():
      raise ValueError("The transformation matrix is not invertible")

    self._homogeneous_form = homogeneous_form
    self._transformation_matrix = transformation_matrix
    G = _apply_matrix(transformation_matrix, homogeneous_form)
    self._multi_indices = list(G.dict().keys())


  def __repr__(self):
    linear_functions = []
    for multi_index in self.multi_indices():
      linear_functions.append(-sum(multi_index[i] * var('w' + str(i))
                                   for i in range(self.dimension())))
    return "max" + str(tuple(linear_functions)) + "/||w||"


  def __call__(self, position):
    r"""
    Return self.evaluate(position).

    EXAMPLES:
    sage: R.<x0,x1,x2> = GF(3)[]
    sage: f = x0^2 + x1^2 + x2^2
    sage: T = matrix(GF(3), [[1,0,2],[1,1,0],[0,2,1]]); T
    [1 0 2]
    [1 1 0]
    [0 2 1]
    sage: mu = ApartmentSphericalStabilityFunction(f, T)
    sage: mu
    min(2*w0, w0 + w1, 2*w1, w0 + w2, w1 + w2, 2*w2)/||w||
    sage: mu([2,-1,-1])
    -1/3*sqrt(6)
    sage: mu(2)
    -2/3*sqrt(6)*sin(2)
    sage: mu(3/2*pi)
    -1/3*sqrt(6)
    """
    return self.evaluate(position)


  def homogeneous_form(self):
    return self._homogeneous_form


  def transformation_matrix(self):
    return self._transformation_matrix


  def dimension(self):
    return self.homogeneous_form().parent().ngens()


  def multi_indices(self):
    return self._multi_indices


  def evaluate(self, position):
    r"""
    Evaluate `self` at the point given by `position`.

    INPUT:
    - ``position`` -- either a list of rational numbers or an angle in radians.

    EXAMPLES::
      sage: R.<x0,x1,x2> = GF(3)[]
      sage: F = x0^2 + x1^2 + x2^2
      sage: T = matrix(GF(3), [[1,2,0], [1,1,1], [0,1,1]])
      sage: mu = ApartmentSphericalStabilityFunction(F, T)
      sage: mu([1,1,-2]); RR(mu([1,1,-2])); RR(mu(pi/6 + 1*pi/3))
      2/3*sqrt(6)
      1.63299316185545
      1.63299316185545
      sage: mu([1,-2,1]); RR(mu([1,-2,1])); RR(mu(pi/6 + 5*pi/3))
      1/6*sqrt(6)
      0.408248290463863
      0.408248290463863
    """
    if isinstance(position, list):
      weight_vector = [QQ(i) for i in position]
      weight_vector_norm = vector(QQ, weight_vector).norm()
      if weight_vector_norm == 0:
        raise ValueError("Weight vector cannot be the zero vector.")
      L = []
      for multi_index in self.multi_indices():
        f = sum(multi_index[j] * weight_vector[j]
                for j in range(self.dimension()))
        L.append(-f / weight_vector_norm)
      return max(L)

    if self.dimension() != 3:
      raise NotImplementedError

    return max(A1 * cos(position) + A2 * sin(position)
               for A1, A2, _ in self._L_theta_coefficients())


  def active_functions(self, position):
    r"""
    Return the set of multi indices corresponding to active functions
    of self on the apartment given by transformation_matrix at the
    point given by position.

    EXAMPLES:
    sage: R.<x0,x1,x2> = GF(3)[]
    sage: F = x0^2 + x1^2 + x2^2
    sage: T = matrix(GF(3), [[1,2,0], [1,1,1], [0,1,1]])
    sage: mu = SphericalStabilityFunction(F)
    sage: mu.active_functions(pi/6)
    {(0, 0, 2), (0, 1, 1)}
    sage: mu.active_functions([2,-1,-1])
    {(0, 0, 2), (0, 1, 1)}
    sage: mu.active_functions(pi/2)
    {(0, 0, 2)}
    sage: mu.active_functions([1,1,-2])
    {(0, 0, 2)}
    sage: mu.active_functions(5*pi/6)
    {(0, 0, 2), (1, 0, 1), (2, 0, 0)}
    sage: mu.active_functions([-1,2,-1])
    {(0, 0, 2), (1, 0, 1), (2, 0, 0)}
    sage: mu.active_functions(7*pi/6)
    {(2, 0, 0)}
    sage: mu.active_functions([-2,1,1])
    {(2, 0, 0)}
    sage: mu.active_functions(3*pi/2)
    {(2, 0, 0)}
    sage: mu.active_functions([-1,-1,2])
    {(2, 0, 0)}
    sage: mu.active_functions(11*pi/6)
    {(0, 1, 1)}
    sage: mu.active_functions([1,-2,1])
    {(0, 1, 1)}
    """
    self_value = self(position)
    active_funcs = set()

    if isinstance(position, list):
      weight_vector = [QQ(i) for i in position]
      weight_vector_norm = vector(QQ, weight_vector).norm()
      for multi_index in self.multi_indices():
        L_value = sum(multi_index[j] * weight_vector[j]
                      for j in range(self.dimension())) / weight_vector_norm
        if self_value == -L_value:
          active_funcs.add(multi_index)
      return active_funcs

    if self.dimension() != 3:
      raise NotImplementedError

    for A1, A2, mult_idx in self._L_theta_coefficients():
      L_value = A1 * cos(position) + A2 * sin(position)
      if self_value == L_value:
        active_funcs.add(mult_idx)
    return active_funcs


  def sign_of_maximum(self):
    r"""
    Return the sign of the maximal value of self.

    ToDo: simplify the code and the number of iterations as well as the
    general complexity. Use for that the L1-morm condition ||w||_1 = 1,
    which can be realized by introducing positive and negative variables,
    w_i = w_i_pos - w_i_neg.

    EXAMPLES:
    sage: K = GF(17^3)
    sage: a = K.gen()
    sage: R.<x0,x1,x2> = K[]
    sage: f = x0^3*(x0 + x1 + x2)
    sage: T = matrix(K, [[1,5*a^2 + 1,1],[0,1,a^2 + 2*a + 7],[0,0,1]]); T
    [              1      5*z3^2 + 1               1]
    [              0               1 z3^2 + 2*z3 + 7]
    [              0               0               1]
    sage: mu = ApartmentSphericalStabilityFunction(f, T)
    sage: mu.sign_of_maximum()
    1

    sage: R.<x0,x1,x2> = GF(2)[]
    sage: f = (x0^2 + x1*x2)^2
    sage: T = matrix(GF(2), [[1,0,1],[0,1,0],[0,0,1]]); T
    [1 0 1]
    [0 1 0]
    [0 0 1]
    sage: mu = ApartmentSphericalStabilityFunction(f, T)
    sage: mu.sign_of_maximum()
    0

    sage: R.<x0,x1,x2> = GF(3)[]
    sage: f = x0*x1*x2*(x0 + x1 + x2)
    sage: T = identity_matrix(GF(3), 3)
    sage: mu = ApartmentSphericalStabilityFunction(f, T)
    sage: mu.sign_of_maximum()
    -1
    """

    all_sings = set()
    # positive faces
    for p_position in range(self.dimension()):
      MILP = MixedIntegerLinearProgram(solver='PPL')
      v = MILP.new_variable()
      t = v['minimum']
      MILP.set_objective(t)

      # Conditions to be on [-1,1]x...x[-1,1]x{1}x[-1,1]x...x[-1,1]
      for i in range(self.dimension()):
        if i == p_position:
          MILP.add_constraint(v[i] == 1)
        MILP.add_constraint(-1 <= v[i] <= 1)

      # Condition to be in H
      MILP.add_constraint(sum(v[i] for i in range(self.dimension())) == 0)

      # All linear functions are bounded by minimum,
      for multi_index in self.multi_indices():
        lin_func = sum(Integer(index) * v[j] for j, index in enumerate(multi_index))
        MILP.add_constraint(t <= lin_func)

      MILP.solve()
      values = MILP.get_values(v)
      local_sign = sign(values['minimum'])
      if local_sign == 1:
        return 1
      all_sings.add(local_sign)

    # negative faces
    for n_position in range(self.dimension()):
      MILP = MixedIntegerLinearProgram(solver='PPL')
      v = MILP.new_variable()
      t = v['minimum']
      MILP.set_objective(t)

      # Conditions to be on [-1,1]x...x[-1,1]x{-1}x[-1,1]x...x[-1,1]
      for i in range(self.dimension()):
        if i == n_position:
          MILP.add_constraint(v[i] == -1)
        MILP.add_constraint(-1 <= v[i] <= 1)

      # Condition to be in H
      MILP.add_constraint(sum(v[i] for i in range(self.dimension())) == 0)

      # All linear functions are bounded by minimum,
      for multi_index in self.multi_indices():
        lin_func = sum(Integer(index) * v[j] for j, index in enumerate(multi_index))
        MILP.add_constraint(t <= lin_func)

      MILP.solve()
      values = MILP.get_values(v)
      local_sign = sign(values['minimum'])
      if local_sign == 1:
        return 1
      all_sings.add(local_sign)

    return max(all_sings)


  def plot_cartesian(self, plot_individual_Li=True, base_radius=1.0):
    r"""
    Generates a 3D cylindrical plot.
    """
    if self.dimension() != 3:
      raise NotImplementedError(
          "This plotting method is implemented only for dimension 3."
      )

    d_poly = self.homogeneous_form().degree()

    coeffs_L_theta = []
    sqrt2 = math.sqrt(2)
    sqrt6 = math.sqrt(6)

    for i_vec_sage in self.multi_indices():
      i_vec = [float(c) for c in i_vec_sage]
      i0, i1, i2 = i_vec
      A_1 = (i0 - i1) / sqrt2
      A_2 = (i0 + i1 - 2*i2) / sqrt6
      coeffs_L_theta.append((A_1, A_2))

    def L_theta_numerical_inner(theta_val, A_coeff, B_coeff):
      return A_coeff * math.cos(theta_val) + B_coeff * math.sin(theta_val)

    def f_theta_numerical_inner(theta_val):
      min_val = float('inf')
      if not coeffs_L_theta:
          return 0.0
      for A1, A2 in coeffs_L_theta:
        val = L_theta_numerical_inner(theta_val, A1, A2)
        if val < min_val:
          min_val = val
      return min_val

    theta_angles = np.linspace(0, 2 * math.pi, 400)
    f_values = np.array([f_theta_numerical_inner(t_ang) for t_ang in theta_angles])

    x_coords = base_radius * np.cos(theta_angles)
    y_coords = base_radius * np.sin(theta_angles)
    z_coords_f = f_values

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    num_to_label_Li = 5
    if plot_individual_Li and len(coeffs_L_theta) > 0:
      try:
        colormap = plt.get_cmap('viridis')
        colors = [colormap(i) for i in np.linspace(0, 1, len(coeffs_L_theta))]
      except:
        colors = ['darkgrey'] * len(coeffs_L_theta)

      for k, (A1, A2) in enumerate(coeffs_L_theta):
        z_coords_Lk = np.array([L_theta_numerical_inner(t_ang, A1, A2)
                                for t_ang in theta_angles])
        line_Lk, = ax.plot(x_coords, y_coords, z_coords_Lk,
                           linestyle='--', linewidth=1.0,
                           color=colors[k % len(colors)],
                           alpha=0.6, zorder=1)

    main_line_width = 2.5
    line_f, = ax.plot(x_coords, y_coords, z_coords_f,
                      color='blue', linewidth=main_line_width,
                      alpha=1.0, zorder=2)

    # --- Plotting special apartment vertices ---
    special_thetas = [
        math.pi/6, math.pi/2, 5*math.pi/6,
        7*math.pi/6, 3*math.pi/2, 11*math.pi/6
    ]
    first_special_point = True
    # Make marker size comparable to main line width
    marker_s = main_line_width**2
    
    for theta_star in special_thetas:
      f_theta_star = f_theta_numerical_inner(theta_star)
      x_star = base_radius * math.cos(theta_star)
      y_star = base_radius * math.sin(theta_star)
      z_star = f_theta_star

      if first_special_point:
          sc_point = ax.scatter([x_star], [y_star], [z_star], color='red', s=marker_s,
                                zorder=3, depthshade=False)
          first_special_point = False
      else:
          ax.scatter([x_star], [y_star], [z_star], color='red', s=marker_s, 
                     zorder=3, depthshade=False)
    # --- End of plotting special apartment vertices ---

    ax.view_init(elev=25., azim=-135)
    plt.tight_layout() 
    plt.show()


  def plot_polar(self, plot_individual_Li=True, sectors=None):
    r"""
    Generates a 2D polar plot of `self` restricted to the unit sphere.

    ..MATH::
    Let H = {w \in \RR^{n+1} : w_0 + ... + w_n = 0}.
    Let \mu : H \to \RR be `self`.
    Now consider the orthonormal basis of H given by
    v_1 = 1/sqrt(2)*(1,-1,0),
    v_2 = 1/sqrt(6)*(1,1,-2)
    The postcomposition of \mu with the isomorphism
    A : \RR^2 ---> H, a \mapsto a_1*v_1 + a_2*v_2
    can be described as follows. We have
    \mu(w) = \max_i(-<i, w>) / ||w||
    and therefore for a \in \RR^2 the function value \mu(A(a)) is equal
    to
    \max_i(a1*(i_1 - i_0)/sqrt(2), a2*(2*i_2 - i_0 - i_1)/sqrt(6)).
    The graph of \mu \circ A restricted to the unit sphere S^1 can be
    described via the parameterization
    [0, 2pi) \to \RR, theta \mapsto \mu \circ A(cos(theta), sin(theta)).

    The vertices of the apartment correspond to the following angles:
    (1)   \pi / 6 ~ ( 2, -1, -1)
    (2)   \pi / 2 ~ ( 1,  1, -2)
    (3)  5\pi / 6 ~ (-1,  2, -1)
    (4)  7\pi / 6 ~ (-2,  1,  1)
    (5)  3\pi / 2 ~ (-1, -1,  2)
    (6) 11\pi / 6 ~ ( 1, -2,  1)
    """
    if self.dimension() != 3:
      raise NotImplementedError(
          "This plotting method is implemented only for dimension 3."
      )

    d_poly = self.homogeneous_form().degree()
    coeffs_L_theta = []
    sqrt2 = math.sqrt(2)
    sqrt6 = math.sqrt(6)

    for i_vec_sage in self.multi_indices():
      i_vec = [float(c) for c in i_vec_sage]
      i0, i1, i2 = i_vec
      A_1 = (i1 - i0) / sqrt2
      A_2 = (2*i2 - i0 - i1) / sqrt6
      coeffs_L_theta.append((A_1, A_2))

    def L_theta_numerical_inner(theta_val, A_coeff, B_coeff):
      return A_coeff * math.cos(theta_val) + B_coeff * math.sin(theta_val)

    def f_theta_numerical_inner(theta_val):
      min_val = float('inf')
      if not coeffs_L_theta: 
          return 0.0
      for A1, A2 in coeffs_L_theta:
        val = L_theta_numerical_inner(theta_val, A1, A2)
        if val < min_val:
          min_val = val
      return min_val

    theta_angles = np.linspace(0, 2 * math.pi, 400)
    f_values = np.array([f_theta_numerical_inner(t_ang) for t_ang in theta_angles])

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})

    num_to_label_Li = 5
    if plot_individual_Li and len(coeffs_L_theta) > 0:
      try:
        colormap = plt.get_cmap('viridis')
        colors = [colormap(i) for i in np.linspace(0, 1, len(coeffs_L_theta))]
      except:
        colors = ['darkgrey'] * len(coeffs_L_theta)

      for k, (A1, A2) in enumerate(coeffs_L_theta):
        r_coords_Lk = np.array([L_theta_numerical_inner(t_ang, A1, A2) for t_ang in theta_angles])
        ax.plot(theta_angles, r_coords_Lk, linestyle='--',
                linewidth=1.0, color=colors[k % len(colors)])

    main_line_width = 2.0
    ax.plot(theta_angles, f_values, color='blue',
            linewidth=main_line_width, alpha=1.0)

    # --- Sector Shading Logic ---
    if sectors is not None:
      try:
        sec_a, sec_b = sectors
        # Map indices 1..6 to their respective angles
        # Formula based on vertices: angle_k = pi/6 + (k-1)*pi/3
        # 1 -> pi/6, 2 -> pi/2, etc.
        angle_a = math.pi/6 + (sec_a - 1) * (math.pi/3)
        angle_b = math.pi/6 + (sec_b - 1) * (math.pi/3)

        # Generate theta range for the sector
        sector_thetas = np.linspace(angle_a, angle_b, 200)

        # Get current plot limits to fill the background completely
        r_min, r_max = ax.get_ylim()

        # Shade the sector. zorder=0 places it behind the main plot lines.
        ax.fill_between(sector_thetas, r_min, r_max, color='lightgrey', alpha=0.5, zorder=0)

      except (ValueError, TypeError):
        print("Warning: 'sectors' must be a tuple of two integers.")

    # --- Plotting special apartment vertices ---
    special_thetas = [
        math.pi/6, math.pi/2, 5*math.pi/6,
        7*math.pi/6, 3*math.pi/2, 11*math.pi/6
    ]
    marker_s = main_line_width**3

    for theta_star in special_thetas:
      r_star = f_theta_numerical_inner(theta_star)
      ax.scatter([theta_star], [r_star], color='red', s=marker_s, zorder=3)
    # --- End of plotting special apartment vertices ---

    plt.tight_layout()
    plt.show()


  def _L_theta_coefficients(self):
    r"""
    Return the coefficients of the parameterization of `self`
    on the unit sphere.

    ..MATH::
    Assume self.dimension() == 3.
    Let F = self.homogeneous_form() \in K[x_0, x_1, x_2].
    Let d = F.degree().
    Let I = {i \in NN^3 : i_0 + i_1 + i_2 = d}.
    Let G = _apply_matrix(transformation_matrix, F).
    Now we can write
    G = \sum_{i \in J} a_i x_0^{i_0} x_1^{i_1} x_2^{i_2}
    with
    J = {i : a_i \neq 0} \subset I.
    Let H = {w \in \RR^3 : w_0 + w_1 + w_2 = 0}.
    Let S_H = {w \in H : ||w|| = 1}.
    Then `self` is given by the function
    \mu(w) = max_{i\in J}(-<i, w> / ||w||),
    where <i,w> = i_0*w_0 + i_1*w_1 + i_2*w_2, for w \in H.
    Now we consider the orthonormal basis
    v_1 = 1/sqrt(2)*(1,-1,0),
    v_2 = 1/sqrt(6)*(1,1,-2)
    of H. So, we have the isomorphism
    A : \RR^2 ---> H, a \mapsto a_1*v_1 + a_2*v_2.
    For any i \in J we obtain the linear function
    L_i : \RR^2 ---> \RR, a \mapsto -<i, A(a)>.
    The graph of the restriction of L_i to the unit sphere S^1 can be
    described vie the parameterization
    [0, 2pi) ---> \RR, theta \mapsto L_i(cos(theta), sin(theta)).
    Note that this graph is the intersection of the plane given by the
    graph of L_i with the unit cylinder, i.e. is an ellipse. We have
    L_i(cos(theta), sin(theta)) = A_1*cos(theta) + A_2*sin(theta)
    for
    A_1 = (i_1 - i_0)/sqrt(2) and A_2 = (2*i_2 - i_0 - i_1)/sqrt(6).
    """
    if self.dimension() != 3:
      raise NotImplementedError

    coeffs_L_theta = []
    sqrt2 = sqrt(QQ(2))
    sqrt6 = sqrt(QQ(6))

    for _mult_idx in self.multi_indices():
      mult_idx = [QQ(c) for c in _mult_idx]
      i0, i1, i2 = mult_idx
      A_1 = (i1 - i0) / sqrt2
      A_2 = (2*i2 - i0 - i1) / sqrt6
      coeffs_L_theta.append((A_1, A_2, mult_idx))

    return coeffs_L_theta



# ================== helper functions ==================

def _apply_matrix(T, F, affine_patch = None):
    """
    Return F((x_0,...,x_n) * T) or its dehomogenization
    at affine_patch, i.e. x_{affine_patch} = 1 if
    affine_patch != None

    INPUT:
        T            - matrix over K
        F            - polynomial in K[x_0,...,x_n]
        affine_patch - integer between 0 and n

    OUTPUT:
        F((x_0,...,x_n) * T) with x_{affine_patch} = 1 if
        affine_patch != None

    MATHEMATICAL INTERPRETATION:
        ToDo...
    """

    generators = list(F.parent().gens())
    if affine_patch != None:
        generators[affine_patch] = F.parent()(1)
    return F(list( vector(generators) * T ))


def _unipotent_lower_triangular_matrices(R, n):
  r"""
  Creates an iterator for all (n+1)x(n+1) unipotent lower triangular
  matrices over the finite ring R.

  INPUT:
    K - finite ring
    n - positive integer defining the matrix size (n+1)

  YIELDS:
    The next (n+1)x(n+1) unipotent lower triangular matrix
  """
  dim = n + 1
  # The positions of the entries below the main diagonal are fixed.
  # In 0-indexed coordinates, these are the positions (i, j) where i > j.
  lower_triangular_indices = []
  for i in range(dim):
    for j in range(i):
      lower_triangular_indices.append((i, j))
  num_free_entries = len(lower_triangular_indices) # This is n*(n+1)/2

  # Create an iterator over the Cartesian product R x R x ... x R.
  combinations_iterator = R**num_free_entries

  # Iterate through each unique combination of values
  for combo in combinations_iterator:
    # Start with a fresh identity matrix for each combination
    M = identity_matrix(R, dim)
    # Fill the lower triangular part with the values from the current combination
    for i, pos in enumerate(lower_triangular_indices):
      M[pos] = combo[i]
    yield M
