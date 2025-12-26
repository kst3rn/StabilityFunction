from copy import copy
import itertools
from sage.matrix.constructor import matrix, zero_matrix, identity_matrix
from sage.modules.free_module_element import vector
from sage.rings.ring import Ring


def _apply_matrix(T, F, i=None):
  r"""
  Return F((x_0,...,x_n) * T) or its dehomogenization
  at `i`, i.e. x_{i} = 1 if `i` is not `None`.

  INPUT:
  - ``T`` -- matrix over K.
  - ``F`` -- polynomial in K[x_0,...,x_n].
  - ``i`` -- an integer between 0 and n or `None`.

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
  num_gens = F.parent().ngens()
  if not (T.nrows() == num_gens and T.ncols() == num_gens):
    raise ValueError(f"Matrix T must be {num_gens}x{num_gens}")
  generators = list(F.parent().gens())
  if i is not None:
    if not (0 <= i < num_gens):
      raise ValueError(f"i = {i} is not between 0 and {num_gens-1}")
    generators[i] = F.parent()(1)
  return F( list(vector(generators) * T) )


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
  T2 = _integral_plane_transformation(_apply_matrix(T1, linear_form), weight_vector)[0]

  return T2 * T1


def _move_point_and_line_to_001_and_x0(base_ring, P, L):
  r"""
  Return an invertible matrix `T` over `base_ring` such that
  (0,0,1)*T = P and the linear form L((x0,x1,x2)*T) is equal
  to C*x_0.

  INPUT:
  - ``base_ring`` -- a ring.
  - ``P`` -- a list/tuple of 3 coordinates representing a point.
  - ``L`` -- a linear form A0*x0 + A1*x1 + A2*x2.

  OUTPUT:
  An invertible 3x3 matrix.

  EXAMPLES::
    sage: R.<x0,x1,x2> = QQ[]
    sage: _move_point_and_line_to_001_and_x0(QQ, [0,0,1], x0)
    [1 0 0]
    [0 1 0]
    [0 0 1]
    sage: P = [0,0,1]
    sage: L = x1
    sage: T = _move_point_and_line_to_001_and_x0(QQ, P, L); T
    [0 1 0]
    [1 0 0]
    [0 0 1]
    sage: vector([0,0,1])*T
    (0, 0, 1)
    sage: L(list(vector([x0,x1,x2]) * T))
    x0
    sage: P = [3,2,1]
    sage: L = 2*x0 - 3*x1
    sage: T = _move_point_and_line_to_001_and_x0(QQ, P, L); T
    [1 0 0]
    [3 2 0]
    [3 2 1]
    sage: vector([0,0,1])*T
    (3, 2, 1)
    sage: L(list(vector([x0,x1,x2]) * T))
    2*x0
  """

  if not isinstance(base_ring, Ring):
    raise ValueError(f"{base_ring} is not a ring")

  R = L.parent()
  if R.base_ring() is not base_ring:
    raise ValueError(f"The base ring of {L} is not equal to {base_ring}")

  P = [base_ring(x) for x in P]
  A = [L.coefficient(x) for x in R.gens()]

  if not any(P):
    raise ValueError(f"The point must be nonzero. Provided: {P}")
  if not any(A):
    raise ValueError(f"The linear form must be nonzero. Provided: {L}")

  val_at_P = sum(a*p for a, p in zip(A,P))
  if val_at_P != 0:
    raise ValueError(f"The homogeneous form {L} does not vanish at {P}.")

  k = -1
  for idx in range(3):
    if A[idx] != 0:
      k = idx
      break

  row0 = [base_ring(0)] * 3
  row0[k] = base_ring(1)
  others = [idx for idx in range(3) if idx != k]
  i, j = others[0], others[1]
  row1 = [base_ring(0)] * 3

  if P[j] != 0:
    row1[k] = -A[i]
    row1[i] =  A[k]
    row1[j] =  0
  else:
    row1[k] = -A[j]
    row1[j] =  A[k]
    row1[i] =  0

  return matrix(base_ring, [row0, row1, P])


def _unipotent_lower_triangular_matrices(R, n):
  r"""
  Creates an iterator for all n by n unipotent lower triangular
  matrices over the finite ring `R`.

  INPUT:
  - ``R`` -- finite ring.
  - ``n`` -- positive integer defining the matrix size `n`.

  YIELDS:
  The next n by n unipotent lower triangular matrix.
  """
  indices = [(r, c) for r in range(n) for c in range(r)]
  template = identity_matrix(R, n)
  iterator = itertools.product(R, repeat=len(indices))
  for values in iterator:
    if not any(values):
      continue
    M = copy(template)
    for (r, c), val in zip(indices, values):
      M[r, c] = val
    yield M


def _unipotent_integral_matrices(R, n, weight_vector):
  r"""
  Creates an iterator for certain n by n unipotent matrices T
  over the finite ring `R` such that the implication
    (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
  holds.
  """
  P = _sorting_permutation_matrix(weight_vector)
  P_inverse = P.transpose()
  for M in _unipotent_lower_triangular_matrices(R, n):
    yield P * M * P_inverse


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
