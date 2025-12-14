
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from transformations import _apply_matrix
from sage.all import *



class LinearValuation:
  r"""
  Construct the linear valuation on `polynomial_ring` diagonalized
  by the weight system given by the basis
  `(x_0,...,x_n) * base_change_matrix.inverse()`
  and the weight vector `weight_vector`.

  INPUT:
  - ``polynomial_ring`` -- a polynomial ring R[x_0,...,x_n].
  - ``base_valuation`` -- a valuation on the base ring of `R`.
  - ``base_change_matrix`` -- an invertible matrix over `R`.
  - ``weight_vector`` -- a tuple of rational numbers (w_0,...,w_n).

  OUTPUT:
  - ``v_{E,w}`` -- the linear valuation on `polynomial_ring` with
                    w = weight_vector,
                    E = (x_0,...,x_n)*base_change_matrix.inverse().

  EXAMPLES::
    sage: K = QQ
    sage: R.<x0,x1,x2> = K[]
    sage: v_K = K.valuation(3)
    sage: T = matrix(K, [[1,0,0],[1,1,0],[0,2,1]]); T
    [1 0 0]
    [1 1 0]
    [0 2 1]
    sage: w = [2,0,5]
    sage: v = LinearValuation(R, v_K, T, w); v
    Linear valuation on Multivariate Polynomial Ring in x0, x1, x2 over Rational Field

  The valuation is diagonalized by the basis self.adapted_basis().
    sage: v.adapted_basis()
    (x0 - x1 + 2*x2, x1 - 2*x2, x2)
    sage: y0, y1, y2 = v.adapted_basis()
    sage: [v(y0), v(y1), v(y2)] == w
    True

  The valuation of the monomials in y0, y1, y2 is a consequence
  of the axiom v(a*b) = v(a) + v(b).
    sage: v(3*y0*y2^2) == v_K(3) + w[0] + 2*w[2]
    True
    sage: v(9*y0^7*y1) == v_K(9) + 7*w[0] + w[1]
    True

  The valuation of the sum of monomials in y0, y1, y2 is equal to
  the minimum of valuations of these monomials.
    sage: v(3*y0*y2^2 + 9*y0^7*y1) == min(v(3*y0*y2^2), v(9*y0^7*y1))
    True

  To compute the valuation of a polynomial f in x0, x1, x2 one has to
  compute its coefficients with respect to y0, y1, y2. This is done by
  acting with `T` on the variables from the right, i.e. f((x0,x1,x2)*T)).
    sage: v.standard_basis()
    (x0, x1, x2)
    sage: f = x1^2 + 2*x1*x2
    sage: g = f(list(vector([x0,x1,x2]) * T)); g
    x1^2 + 6*x1*x2 + 8*x2^2
    sage: g(y0,y1,y2) == f
    True
    sage: v(f) == min(v(y1^2), v(6*y1*y2), v(8*y2^2))
    True

  ..MATH::
  First, let
    K   = polynomial_ring.base_ring(),
    E_0 = (x_0,...,x_n) = polynomial_ring.gens(),
    v_K = base_valuation,
    T   = base_change_matrix,
    B   = base_change_matrix.inverse(),
    w   = weight_vector.
  Thus, we can write
    polynomial_ring = K[x_0,...,x_n].
  We call E_0 the standard basis and consider T and B as linear
  transformations, with respect to E_0, i.e.
    T(x_j) = sum_{i=0}^n a_{ij}*x_i  and  B(x_j) = sum_{i=0}^n t_{ij}*x_i .
  Then E = (y_0,...,y_n) = ( B(x_0),...,B(x_n) ) is a new basis of K[x_0,...,x_n].
  and we set
    LinearValuation( K[x_0,...,x_n], v_K, T, w ) = v_{E,u} .
  Now let F be a polynomial in K[x_0,...,x_n]. To compute v_{E,u}(F) we have to
  write F with respect to the basis E. This will be done as follows. If we view
    E_0 = (x_0,...,x_n)
  as a vector in Sage, we get
    (x_0,...,x_n)*B = (sum_{i=0}^n t_{i,0}*x_i,...,sum_{i=0}^n t_{i,n}*x_i)
                    = (y_0,...,y_n)
  and therefore
    F(x_0,...,x_n) = F( (y_0,...,y_n)*B^{-1} ) = F( (y_0,...,y_n)*T ).
  Thus, the polynomial
    G(y_0,...,y_n) = F( (y_0,...,y_n)*T ) in K[y_0,...,y_n]
  describes F with respect to the basis (y_0,...,y_n) and T describes the base change.
  Thus, with the notation
    <i,u> = i_0*u_0 + ... + i_n*u_n,
  we obtain
    v_{E,u}(F) = min( v_K(a_i) + <i,u> : i in I ) with G = sum_{i in I} a_i y^i,
  where i is a multi-index, i.e. I is a subset of NN^{n+1}.
  """

  def __init__(self,
               polynomial_ring,
               base_valuation,
               base_change_matrix,
               weight_vector):
    r"""
    Construct the linear valuation on `polynomial_ring` diagonalized
    by the weight system given by the basis
    `(x_0,...,x_n) * base_change_matrix.inverse()`
    and the weight vector `weight_vector`.

    INPUT:
    - ``polynomial_ring`` -- a polynomial ring R[x_0,...,x_n].
    - ``base_valuation`` -- a valuation on the base ring of `R`.
    - ``base_change_matrix`` -- an invertible matrix over `R`.
    - ``weight_vector`` -- a tuple of rational numbers (w_0,...,w_n).

    OUTPUT:
    - ``v_{E,w}`` -- the linear valuation on `polynomial_ring` with
                     w = weight_vector,
                     E = (x_0,...,x_n)*base_change_matrix.inverse().

    EXAMPLES::
      sage: K = QQ
      sage: R.<x0,x1,x2> = K[]
      sage: v_K = K.valuation(3)
      sage: T = matrix(QQ, [[1,0,0],[1,1,0],[0,2,1]]); T
      [1 0 0]
      [1 1 0]
      [0 2 1]
      sage: w = [2,0,5]
      sage: v = LinearValuation(R, v_K, T, w); v
      Linear valuation on Multivariate Polynomial Ring in x0, x1, x2 over Rational Field

    The valuation is diagonalized by the basis self.adapted_basis().
      sage: v.adapted_basis()
      (x0 - x1 + 2*x2, x1 - 2*x2, x2)
      sage: y0, y1, y2 = v.adapted_basis()
      sage: [v(y0), v(y1), v(y2)] == w
      True

    The valuation of the monomials in y0, y1, y2 is a consequence
    of the axion v(a*b) = v(a) + v(b).
      sage: v(3*y0*y2^2) == v_K(3) + w[0] + 2*w[2]
      True
      sage: v(9*y0^7*y1) == v_K(9) + 7*w[0] + w[1]
      True

    The valuation of the sum of monomials in y0, y1, y2 is equal to
    the minimum of valuations of these monomials.
      sage: v(3*y0*y2^2 + 9*y0^7*y1) == min(v(3*y0*y2^2), v(9*y0^7*y1))
      True

    To compute the valuation of a polynomial f in x0, x1, x2 one has to
    compute its coefficients with respect to y0, y1, y2. This is done by
    acting with `T` on the variables from the right, i.e. f((x0,x1,x2)*T)).
      sage: v.standard_basis()
      (x0, x1, x2)
      sage: f = x1^2 + 2*x1*x2
      sage: g = f(list(vector([x0,x1,x2]) * T)); g
      x1^2 + 6*x1*x2 + 8*x2^2
      sage: g(y0,y1,y2) == f
      True
      sage: v(f) == min(v(y1^2), v(6*y1*y2), v(8*y2^2))
      True

    ..MATH::
    First, let
      K   = polynomial_ring.base_ring(),
      E_0 = (x_0,...,x_n) = polynomial_ring.gens(),
      v_K = base_valuation,
      T   = base_change_matrix,
      B   = base_change_matrix.inverse(),
      w   = weight_vector.
    Thus, we can write
      polynomial_ring = K[x_0,...,x_n].
    We call E_0 the standard basis and consider T and B as linear
    transformations, with respect to E_0, i.e.
      T(x_j) = sum_{i=0}^n a_{ij}*x_i  and  B(x_j) = sum_{i=0}^n t_{ij}*x_i .
    Then E = (y_0,...,y_n) = ( B(x_0),...,B(x_n) ) is a new basis of K[x_0,...,x_n].
    and we set
      LinearValuation( K[x_0,...,x_n], v_K, T, w ) = v_{E,u} .
    Now let F be a polynomial in K[x_0,...,x_n]. To compute v_{E,u}(F) we have to
    write F with respect to the basis E. This will be done as follows. If we view
      E_0 = (x_0,...,x_n)
    as a vector in Sage, we get
      (x_0,...,x_n)*B = (sum_{i=0}^n t_{i,0}*x_i,...,sum_{i=0}^n t_{i,n}*x_i)
                      = (y_0,...,y_n)
    and therefore
      F(x_0,...,x_n) = F( (y_0,...,y_n)*B^{-1} ) = F( (y_0,...,y_n)*T ).
    Thus, the polynomial
      G(y_0,...,y_n) = F( (y_0,...,y_n)*T ) in K[y_0,...,y_n]
    describes F with respect to the basis (y_0,...,y_n) and T describes the base change.
    Thus, with the notation
      <i,u> = i_0*u_0 + ... + i_n*u_n,
    we obtain
      v_{E,u}(F) = min( v_K(a_i) + <i,u> : i in I ) with G = sum_{i in I} a_i y^i,
    where i is a multi-index, i.e. I is a subset of NN^{n+1}.
    """

    self._domain = polynomial_ring
    self._base_valuation = base_valuation
    self._base_change_matrix = base_change_matrix
    self._weight_vector = vector(QQ, weight_vector)


  def __repr__(self):
    return f"Linear valuation on {self.domain()}"


  def __call__(self, polynomial):
    return self.evaluate(polynomial)


  def domain(self):
    return self._domain


  def base_valuation(self):
    return self._base_valuation


  def base_change_matrix(self):
    return self._base_change_matrix


  def weight_vector(self):
    return self._weight_vector


  def base_ring(self):
    return self.domain().base_ring()


  def domain_ngens(self):
    return len(self.weight_vector())


  def standard_basis(self):
    return self.domain().gens()


  def adapted_basis(self):
    T = self.base_change_matrix().inverse()
    return tuple(_apply_matrix(T, x) for x in self.standard_basis())


  def base_value_group(self):
    return self.base_valuation().value_group()


  def base_residue_ring(self):
    return self.base_valuation().residue_ring()


  def evaluate(self, f):
    r"""
    Return the value of `self` at `f`.

    INPUT:
    - ``f`` -- a polynomial in the domain of `self`.

    OUTPUT:
    - a rational number

    ..MATH::
    First, let
      K   = self.polynomial_ring.base_ring(),
      E_0 = (x_0,...,x_n) = self.polynomial_ring.gens(),
      v_K = self.base_valuation,
      T   = self.base_change_matrix,
      B   = self.base_change_matrix.inverse(),
      u   = self.weight_vector.
    Thus, we can write
      self.domain() = K[x_0,...,x_n].
    We call E_0 the standard basis and consider T and B as linear
    transformations, with respect to E_0, i.e.
      T(x_j) = sum_{i=0}^n a_{ij}*x_i 
    and
      B(x_j) = sum_{i=0}^n t_{ij}*x_i.
    Then,
      E = (y_0,...,y_n) = (B(x_0),...,B(x_n))
    is a new basis of K[x_0,...,x_n] and
      LinearValuation( K[x_0,...,x_n], v_K, T, u ) = v_{E,u}.
    Let F be a polynomial in K[x_0,...,x_n]. To evaluate v_{E,u} at
    F we have to write F with respect to the basis E. This will be
    done as follows. We have
      (x_0,...,x_n)*B = (sum_{i=0}^n B_{i,0}*x_i,...,sum_{i=0}^n B_{i,n}*x_i)
                      = (y_0,...,y_n)
    and therefore
      F(x_0,...,x_n) = F((y_0,...,y_n)*B^{-1}) = F((y_0,...,y_n)*T).
    Thus, the polynomial
      G(y_0,...,y_n) = F((y_0,...,y_n)*T) in K[y_0,...,y_n]
    describes F with respect to the basis (y_0,...,y_n) and T describes the base change.
    For a multi index set I subset NN^(n+1) we can write
      G = sum_{i in I} a_i y_0^i_0 * ... * y_n^i_n
    Thus, with the notation
      <i,u> = i_0*u_0 + ... + i_n*u_n,
    we obtain
      v_{E,u}(F) = min( v_K(a_i) + <i,u> : i in I ).
    """

    if f == 0:
      return +Infinity
    else:
      F = self.domain()(f)

    N = self.domain_ngens()
    G = _apply_matrix(self.base_change_matrix(), F)
    values = set()
    for multi_index, G_coefficient in G.dict().items():
      value = self.base_valuation()(G_coefficient)
      value += vector(QQ, multi_index) * self.weight_vector()
      values.add(value)
    return min(values)


  def initial_form(self, f):
    r"""
    Return the initial form of `f` with respect to `self`.

    INPUT:
    - ``f`` -- element in the domain of `self`.

    OUTPUT:
    A string.
    """

    import re
    f_valuation = self(f)
    G = _apply_matrix(self.base_change_matrix(), f)

    if G.is_zero():
      return "0"

    gens = self.standard_basis()
    G_initial = self.domain()(0)
    N = self.domain_ngens()
    for multi_index, coefficient in G.dict().items():
      value = self.base_valuation()(coefficient)
      value += sum(multi_index[j] * self._weight_vector[j] for j in range(N))
      if value == f_valuation:
        G_initial += coefficient * prod(gens[j]**multi_index[j] for j in range(N))

    initial_form_str = str(G_initial)
    var_names = [str(g) for g in gens]
    sorted_var_names = sorted(var_names, key=len, reverse=True)
    for v_name in sorted_var_names:
      pattern = r'\b' + re.escape(v_name) + r'\b'
      replacement = f"T({v_name})"
      initial_form_str = re.sub(pattern, replacement, initial_form_str)
    return initial_form_str


  def graded_reduction(self, f):
    r"""
    Return the graded reduction of `f` with respect to `self`.

    INPUT:
    - ``f`` -- element of the domain of `self`.

    EXAMPLES:

    MATHEMATICAL INTERPRETATION:
    First, let
      K   = self.base_ring(),
      E_0 = (x_0,...,x_n) = self.standard_basis(),
      v_K = self.base_valuation,
      k   = v_K.residue_field(),
      p_K = v_K.uniformizer()
      T   = self.base_change_matrix,
      B   = self.base_change_matrix.inverse(),
      u   = self.weight_vector.
    Thus, we can write
      self.domain() = K[x_0,...,x_n].
    We call E_0 the standard basis and consider T and B as linear
    transformations, with respect to E_0, i.e.
      T(x_j) = sum_{i=0}^n a_{ij}*x_i  and  B(x_j) = sum_{i=0}^n t_{ij}*x_i.
    Then E = (y_0,...,y_n) = (B(x_0),...,B(x_n)) is a new basis of K[x_0,...,x_n]
    and
      LinearValuation( K[x_0,...,x_n], v_K, T, u ) = v_{E,u}.
    The graded reduction ring of v_{E,u} is the graded ring
      R = k[t,t^(-1)][Y_0,...,Y_n]
    where
      t = [p_K], Y_0 = [y_0],..., Y_n = [y_n]
    are the graded reduction of
      p_K, y_0,..., y_n,
    respectively. Moreover, the grading of R is given by the degrees of
    the generators
      |t| = v_K(p_K), |Y_0| = v_{E,u}(y_0),..., |Y_n| = v_{E,u}(y_n).
    To compute the graded reduction  [F] of F we first have to write F with
    respect to the basis E. This will be done as follows. We have
      (x_0,...,x_n)*B = (sum_{i=0}^n B_{i,0}*x_i,...,sum_{i=0}^n B_{i,n}*x_i)
                      = (y_0,...,y_n)
    and therefore
      F(x_0,...,x_n) = F( (y_0,...,y_n)*B^(-1) ) = F( (y_0,...,y_n)*T ).
    Since T = B^(-1), the polynomial
      G(y_0,...,y_n) = F( (y_0,...,y_n)*T ) in K[y_0,...,y_n]
    describes F with respect to the basis (y_0,...,y_n) and T describes the
    base change. Note that mathematically the rings
      K[x_0,...,x_n] and K[y_0,...,y_n]
    are the same as well as the polynomials F and G are the same. We only
    perform a change of coordinates. So, for a multi index set I subset
    NN^{n+1} let
      F = G = sum_{i in I} g_i y_0^i_0 * ... * y_n^i_n.
    Now, with the notation
      <i,u> = i_0*u_0 + ... + i_n*u_n
    and
      J = {i in I : v_{E,u}(F) = v_{E,u}(g_i y_0^i_0 * ... * y_n^i_n)},
    we have
      v_{E,u}(F) = min(v_K(g_i) + <i,u> : i in I)
    and therefore the graded reduction [F] is given by
      sum_{i in J} [g_i / p_K^(v_{E,u}(F)-<i,u>)]*t^(v_{E,u}(F)-<i,u>) * Y_0^i_1*...*Y_n^i_n,
    where
      [g_i / p_K^(v_{E,u}(F)-<i,u>)]
    is the reduction of
      g_i / p_K^(v_{E,u}(F)-<i,u>)
    with respect to v_K.
    """

    return GradedReduction(f, self)



class GradedReduction:
  r"""
  Construct ...
  """

  def __init__(self, polynom, linear_valuation):
    r"""
    Construct ...
    """

    self._polynom = linear_valuation.domain()(polynom)
    self._linear_valuation = linear_valuation
    self._weight_vector = linear_valuation.weight_vector()
    self._base_residue_ring = linear_valuation.base_residue_ring()
    self._base_ring_grading = linear_valuation.base_value_group().gen()

    n_w_v = [u / self._base_ring_grading for u in self._weight_vector]
    self._normalized_weight_vector  = tuple(n_w_v)

    N = linear_valuation.domain_ngens()
    G = _apply_matrix(linear_valuation.base_change_matrix(), self._polynom)
    K = LaurentPolynomialRing(self._base_residue_ring, 't')
    grR = PolynomialRing(K, N, 'y')
    nrR = PolynomialRing(self._base_residue_ring, N, 'z')
    graded_reduction = grR(0)
    normalized_reduction = nrR(0)
    polynom_value = linear_valuation(polynom)
    v_K = linear_valuation.base_valuation()
    pi_K = v_K.uniformizer()
    for mult_index, G_coeff in G.dict().items():
      coeff_value = v_K(G_coeff)
      term_value = coeff_value + vector(QQ, mult_index) * self._weight_vector
      if term_value == polynom_value:
        coeff_nr_value = coeff_value / self._base_ring_grading
        pi_K_power = pi_K**coeff_nr_value
        coeff_nr_red = v_K.reduce(G_coeff / pi_K_power)
        nr_monom = prod(nrR.gen(j)**mult_index[j] for j in range(N))
        normalized_reduction += coeff_nr_red * nr_monom
        coeff_gr_red = coeff_nr_red * K.gen()**coeff_nr_value
        gr_monom = prod(grR.gen(j)**mult_index[j] for j in range(N))
        graded_reduction += coeff_gr_red * gr_monom
    self._graded_reduction = graded_reduction
    self._normalized_reduction = normalized_reduction


  def __repr__(self):
    return str(self.graded_reduction_polynomial())


  def lift(self):
    return self._polynom


  def weight_vector(self):
    return self._weight_vector


  def normalized_weight_vector(self):
    return self._normalized_weight_vector


  def base_residue_ring(self):
    return self._base_residue_ring


  def standard_basis(self):
    return self.graded_reduction_ring.gens()


  def base_ring_grading(self):
    return self._base_ring_grading


  def normalized_reduction_polynomial(self):
    return self.r_polynom


  def valuation(self):
    return self._linear_valuation


  def graded_reduction_polynomial(self):
    return self._graded_reduction


  def normalized_reduction_polynomial(self):
    return self._normalized_reduction


  def parent(self):
    return self.graded_reduction_polynomial().parent()


  def graded_base_ring(self):
    return self.graded_reduction_polynomial().base_ring()


  def is_graded_semistable(self):
    r"""
    Return `False` if `self` has a graded instability
    and `True` otherwise.
    """
    from plane_curves import ProjectivePlaneCurve
    reduced_curve = ProjectivePlaneCurve(self.normalized_reduction_polynomial())
    return reduced_curve.is_semistable()


  def is_graded_stable(self):
    r"""
    Return `False` of `self` has a graded semiinstability
    and `True` otherwise.
    """
    from plane_curves import ProjectivePlaneCurve
    reduced_curve = ProjectivePlaneCurve(self.normalized_reduction_polynomial())
    return reduced_curve.is_stable()


  def graded_instability(self, matrix_form='ult'):
    r"""
    Return a graded instability of `self` if it exists
    and `None` otherwise.

    MATHEMATICAL INTERPRETATION:
    First, let
      (y_0,...,y_n) = self.standard_basis(),
      (z_0,...,z_n) = self.RR_standard_basis(),
      f = self.graded_reduction_polynomial(),
      g = self.normalized_reduction_polynomial().
    Let (E,w) with E = (e_0,...,e_n) be an instability of g over k, i.e.
    there exists an invertible matrix
      M in GL_n(k)
    which describes the base change from E to (z_0,...,z_n). So, if we
    view E and (z_0,...,z_n) as a vectors in SageMath, we can write
      E * M = (z_0,...,z_n).
    Note that M, as a matrix over k, desribes a k-linear transforamtion
    of degree zero, i.e. it does not change the grading.
    Now let D be the invertible diagonal matrix
      D = diag(t^(u_0/s), ..., t^(u_n/s)) in GL_3(k[T,T^(-1)]).
    We define the basis
      E_new = (e_0_new,...,e_n_new)
    by
    E_new * D^(-1) = E.
    Note that D describes a k-linear transforamtion which changes the degree
    of homogeneous elements. Furthermore, we obtain
      E_new * D^(-1) * M * D = E * M * D = (z_0,...,z_n) * D
                             = (y_0,...,y_n).
    In particular, (E_new,w) is an instability of f over k[T,T^(-1)]. Note also
    that D^(-1) * M * D describes a k-linear isomorphism of degree zero and
    hence induces an isomorphism of graded rings,
    k[T,T^(-1)][y_0,...,y_n] ---> k[T,T^(-1)][y_0,...,y_n].
    REMARK. Of course we have
      E * M * D = (z_0,...,z_n) * D = (y_0,...,y_n)
    and therefore (E,w) is already an instability of f over k[T,T^(-1)]. But the
    matrix M*D describes a k-linear transforamtion, which shifts degrees. For
    reasons of computational clarity, we want to avoid transformations that
    lead to degree shifts. Thus, we only consider transformations of the form
      D^(-1) * M * D
    as explained above. The (i,j)-th entry of such a matrix is equal to
    m_{ij} * t^{(u_j - u_i) / s}.
    """

    if len(self.normalized_weight_vector()) != 3:
      raise NotImplementedError

    from plane_curves import ProjectivePlaneCurve
    reduced_curve = ProjectivePlaneCurve(self.normalized_reduction_polynomial())
    instability = reduced_curve.instability()
    if instability is None:
      return None
    T = instability.base_change_matrix(matrix_form)
    return GradedInstability(self, T)


  def rational_graded_instability(self, matrix_form = 'ult'):
    r"""
    Return a rational graded instability of `self` if it exists
    and `None` otherwise.
    """

    if len(self.normalized_weight_vector()) != 3:
      raise NotImplementedError
    if matrix_form not in {'ult', 'uut', 'integral'}:
      raise ValueError

    w = self.normalized_weight_vector()
    differences = [w[2] - w[0], w[1] - w[0], w[2] - w[1]]
    if all(x not in ZZ for x in differences):
      return None

    from plane_curves import ProjectivePlaneCurve
    reduced_curve = ProjectivePlaneCurve(self.normalized_reduction_polynomial())

    if all(x in ZZ for x in differences):
      instability = reduced_curve.instability()
      if instability is None:
        return None
      T = instability.base_change_matrix(matrix_form)
      return GradedInstability(self, T)

    for i, j in [(1,0), (2,0), (2,1)]:
      diff = w[j] - w[i]
      if diff in ZZ:
        a = reduced_curve.elementary_instability_direction((i,j))
        if a is not None:
          T = [[1,0,0],[0,1,0],[0,0,1]]
          if matrix_form == 'ult' or matrix_form == 'integral' and diff >= 0:
            T[i][j] = -a
          elif matrix_form == 'uut' or matrix_form == 'integral' and diff < 0:
            T[j][i] = Integer(1) / (-a)
          T = matrix(self.base_residue_ring(), T)
          return GradedInstability(self, T)

    return None


class GradedInstability:

  def __init__(self, graded_reduction, instability_matrix):
    r"""
    INPUT:
    - ``graded_reduction`` - object in the class `GradedReduction`.
    - ``instability_matrix``  - invertible matrix over the degree zero
                                subring of the graded reduction ring
                                representiing a graded instability.
    """

    self.graded_reduction = graded_reduction
    self.graded_base_ring = graded_reduction.graded_base_ring()
    self.graded_reduction_ring = graded_reduction.parent()
    self.base_ring_grading = graded_reduction.base_ring_grading()
    self.grading = graded_reduction.weight_vector()
    self.norm_grading = graded_reduction.normalized_weight_vector()
    self.instability_matrix = instability_matrix


  def __repr__(self):
    return f"Graded Instability of {self.graded_reduction}"


  def is_rational(self):
    r"""
    Return `True` if the instability matrix is defined over
    the graded reduction ring and `False` otherwise.
    """

    for i, row in enumerate(self.instability_matrix):
      for j, entry in enumerate(row):
        if entry != 0:
          t_exponent = self.norm_grading[j] - self.norm_grading[i]
          if t_exponent not in ZZ:
            return False
    return True


  def base_change_matrix(self):
    t = self.graded_base_ring.gen()
    N = len(self.norm_grading)
    graded_trafo_matrix = [[0 for i in range(N)] for j in range(N)]

    for i, row_of_instability_matrix in enumerate(self.instability_matrix):
      for j, entry_in_row in enumerate(row_of_instability_matrix):
        if entry_in_row != 0:
          t_exponent = self.norm_grading[j] - self.norm_grading[i]
          if t_exponent.is_integral():
            graded_trafo_matrix[i][j] = entry_in_row * t**t_exponent
          else:
            return None
    return matrix(self.graded_base_ring, graded_trafo_matrix)


  def lift_matrix(self):
    if not self.is_rational():
      return None

    base_ring = self.graded_reduction.lift().base_ring()
    _valuation_ = self.graded_reduction.valuation()
    base_valuation = _valuation_.base_valuation()
    prime_element = base_valuation.uniformizer()
    N = len(self.grading)
    lifted_graded_trafo_matrix = [[0 for i in range(N)] for j in range(N)]

    for i, row_of_instability_matrix in enumerate(self.instability_matrix):
      for j, entry_in_row in enumerate(row_of_instability_matrix):
        if entry_in_row != 0:
          exponent = self.norm_grading[j] - self.norm_grading[i]
          if exponent.is_integral():
            lifted_graded_trafo_matrix[i][j] = base_valuation.lift(entry_in_row) * prime_element**exponent

    return matrix(base_ring, lifted_graded_trafo_matrix)


  def print_matrix(self):
    r"""
    Just to see the matrix in case if not self.is_rational(). 
    """

    if self.is_rational():
      print(self.base_change_matrix())
    else:
      N = len(self.grading)
      t = var('t')
      array_2d = [[0 for i in range(N)] for j in range(N)]
      for i in range(N):
        for j in range(N):
          t_exponent = (self.grading[j] - self.grading[i]) / self.base_ring_grading
          if self.instability_matrix[i][j] == 0:
            array_2d[i][j] = 0
          elif t_exponent == 0:
            array_2d[i][j] = self.instability_matrix[i][j]
          else:
            array_2d[i][j] = "(" + str(self.instability_matrix[i][j]) + ")*" + str(t) + "^(" + str(t_exponent) + ")"

      v_K = self.graded_reduction.valuation().base_valuation()
      print("Let t be the reduction of the uniformizer " + str(v_K.uniformizer()) + ". Then the base change to instability is given by the following matrix:")

      # Calculate the maximum width of each column
      col_widths = [max( len(str(item)) for item in col ) for col in zip(*array_2d)]
      # Print the array with aligned columns
      for row in array_2d:
        for i, item in enumerate(row):
          print(str(item).ljust(col_widths[i] + 2), end="")  # Add spacing
        print()

