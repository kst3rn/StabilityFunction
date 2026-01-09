from warnings import warn
from itertools import combinations, product
from sage.all import gcd, PolynomialRing, GF, QQ, ZZ, ceil, matrix, GaussValuation, vector, Infinity
from semistable_model.stability import StabilityFunction
from semistable_model.curves import ProjectivePlaneCurve
from semistable_model.stability import minimum_as_valuative_function


def semistable_reduction_field(homogeneous_form,
                               base_ring_valuation,
                               ramification_index=None):
  r"""
  Try to find a minimal extension of the base field
  of `homogeneous_form` such that the projective curve
  defined by `homogeneous_form` over this extension
  has a model with semistable reduction.
  If `ramification_index` is not `None` try to find
  an extension of provided ramification index.
  """
  if not homogeneous_form.is_homogeneous():
    raise ValueError(f"{homogeneous_form} is not homogeneous.")
  if not base_ring_valuation.domain() == QQ:
    raise NotImplementedError(f"The base ring must be {QQ}")
  if not base_ring_valuation.residue_field() == GF(2):
    raise NotImplementedError(f"The residue field is not {GF(2)}")
  if homogeneous_form.degree() != 4:
    warn(
      f"Provided homogeneous form has degree {homogeneous_form.degree()}, but "
      "currently this function is designed for quartics (degree 4). "
      "For other degrees, the algorithm may enter an infinite loop "
      "or take an excessive amount of time.",
      UserWarning,
      stacklevel=2
      )

  if ramification_index is not None:
    return extension_search(homogeneous_form,
                            base_ring_valuation,
                            ramification_index)
  i = 1
  while True:
    L = extension_search(homogeneous_form, base_ring_valuation, 2*i)
    if L is not None:
      return L
    i += 1


def extension_search(homogeneous_form,
                     base_ring_valuation,
                     ramification_index):
  r"""
  Try to find an extension of the base ring of `homogeneous_form`.

  EXAMPLES::
    sage: R.<x,y,z> = QQ[]
    sage: F = y^4 + 2*x^3*z + x*y^2*z + 2*x*z^3
    sage: extension_search(F, QQ.valuation(2), 2)
    Number Field in piK with defining polynomial x^2 + 2
    sage:
    sage: F = 16*x^4 + y^4 + 8*y^3*z + 16*x*y*z^2 + 4*x*z^3
    sage: extension_search(F, QQ.valuation(2), 2)
    None
    sage: extension_search(F, QQ.valuation(2), 4)
    Number Field in piL with defining polynomial x^12 + 2*x^6 + 2
    sage:
    sage: F = 4*x^4 + 4*x*y^3 + y^4 + 2*x*z^3 + 4*y*z^3 + z^4
    sage: extension_search(F, QQ.valuation(2), 4)
    Number Field in piK with defining polynomial x^4 + 2*x^3 + 2*x^2 + 2
    sage:
    sage: F = -2*x^3*y - 12*y^4 - 4*x^3*z - 3*x^2*y*z - 12*y^3*z - 4*x^2*z^2 - 12*x*y*z^2 + 16*y^2*z^2 + 5*y*z^3
    sage: extension_search(F, QQ.valuation(2), 4)
    Number Field in piL with defining polynomial x^16 + 2
  """

  if not homogeneous_form.is_homogeneous():
    raise ValueError(f"{homogeneous_form} is not homogeneous")
  if homogeneous_form.base_ring() != base_ring_valuation.domain():
    raise ValueError(f"The base ring of {homogeneous_form} is not {base_ring_valuation.domain()}")
  if homogeneous_form.base_ring() is not QQ:
    raise NotImplementedError(f"The base ring must be {QQ}")
  if base_ring_valuation.residue_field() is not GF(2):
    raise NotImplementedError(f"The residue field of {base_ring_valuation} is not {GF(2)}")

  K = homogeneous_form.base_ring()
  phi = StabilityFunction(homogeneous_form, base_ring_valuation)
  minimum, btb_point = phi.global_minimum('uut')
  if phi.has_semistable_reduction_at(btb_point):
    if btb_point.is_vertex():
      return K
    piK = base_ring_valuation.uniformizer()
    r_K = base_ring_valuation(piK).denominator()
    r_L = btb_point.ramification_index()
    r = r_L / gcd(r_K, r_L)
    S = PolynomialRing(K, 'x')
    s = S.gen()
    L = K.extension(s**r - piK, 'piL')
    return L.absolute_field('piL')

  R = homogeneous_form.parent()
  S = PolynomialRing(K, 'x')
  v0 = GaussValuation(S, base_ring_valuation)
  R_S = R.change_ring(S)
  F_S = R_S(homogeneous_form)
  s = S.gen()
  step = ZZ(1)/ZZ(ramification_index)
  fixed_valuation = v0.augmentation(s, step)

  w = btb_point.weight_vector()
  w = [QQ(w[i] / step) for i in range(len(w))]
  M = phi.normalized_descent_direction(btb_point, 'integral')
  local_trafo_matrix = [[0,0,0],[0,0,0],[0,0,0]]
  for i, j in product(range(3), range(3)):
    if not M[i][j].is_zero():
      if w[j] - w[i] < 0:
        return None
      local_trafo_matrix[i][j] = M[i][j] * s**(w[j] - w[i])
  local_trafo_matrix = matrix(S, local_trafo_matrix)

  T = btb_point.base_change_matrix()
  global_trafo_matrix = local_trafo_matrix * T
  return _search_tree(F_S, fixed_valuation, step, minimum, global_trafo_matrix, 0, depth_limit=+Infinity)


def _search_tree(F, fixed_valuation, step, minimum, global_trafo_matrix, depth, depth_limit):
  r"""
  Heuristic search.
  """

  depth = depth + 1
  if depth > depth_limit:
    return None

  x0, x1, x2 = F.parent().gens()
  h, e = minimum_as_valuative_function(
    F(list(vector([x0, x1, x2]) * global_trafo_matrix)),
    fixed_valuation)

  local_max_val = h.local_maxima()
  max_local_max = max(a for a, b in local_max_val)
  point_with_biggest_local_max = [b for a, b in local_max_val if a == max_local_max]
  discoids = [b.discoid() for b in point_with_biggest_local_max]
  min_degree_discoid = min(discoids, key=lambda pair: pair[0].degree())
  center, radius = min_degree_discoid
  adjusted_radius = _ceil_step(radius, step)
  center = F.base_ring()(center)
  j = 0
  new_radius = adjusted_radius - j * step
  if new_radius <= fixed_valuation.value_group().gen():
    return None

  while True:
    K = QQ.extension(center, 'piK')
    piK = K.gen()
    R_K = F.parent().change_ring(K)
    F_K = R_K(F)
    phi_typeI = StabilityFunction(F_K, K.valuation(2))
    aI, bI = phi_typeI.local_minimum(_evaluate_matrix(global_trafo_matrix, piK))
    if phi_typeI.has_semistable_reduction_at(bI):
      if bI.minimal_simplex_dimension(ZZ(1) / step) != 0:
        v_K = phi_typeI.base_ring_valuation()
        piK = v_K.uniformizer()
        r_K = v_K(piK).denominator()
        r_L = bI.ramification_index()
        r = r_L / gcd(r_K, r_L)
        S = PolynomialRing(K, 'x')
        s = S.gen()
        L = K.extension(s**r - piK, 'piL')
        return L.absolute_field('piL')
      return K

    new_typeII_valuation = fixed_valuation.augmentation(center, new_radius)
    phi_typeII = StabilityFunction(F, new_typeII_valuation)
    new_minimum, new_btb_point = phi_typeII.local_minimum(global_trafo_matrix)
    if new_minimum >= minimum:
      break
    elif new_btb_point.minimal_simplex_dimension(ZZ(1) / step) == 2:
      continue

    j = j + 1
    new_radius = radius - j * step

    w_normalized = [QQ(x / step) for x in new_btb_point.weight_vector()]
    for i, j in combinations(range(3), 2):
      w_difference = w_normalized[j] - w_normalized[i]
      if w_difference in ZZ:
        local_trafo_matrix = [[1,0,0],[0,1,0],[0,0,1]]
        if w_difference >= 0:
          local_trafo_matrix[i][j] = F.base_ring().gen()**w_difference
        else:
          local_trafo_matrix[j][i] = F.base_ring().gen()**(-w_difference)
        local_trafo_matrix = matrix(F.base_ring(), local_trafo_matrix)
        new_global_trafo_matrix = local_trafo_matrix * global_trafo_matrix
        result = _search_tree(F, fixed_valuation, step, new_minimum, new_global_trafo_matrix, depth, depth_limit)
        if result is not None:
          return result

#    result = _search_tree(F, fixed_valuation, step, new_minimum, global_trafo_matrix, depth, depth_limit)
#    if result is not None:
#      return result



def _ceil_step(x, r):
  r"""
  Return the ceiling of x on the grid r * ZZ.
  Input:
  x: The number to round.
  r: The step size.
  Output:
  A Rational number in ZZ[r]
  """
  if r < 0:
    raise ValueError(f"{r} is not positive")

  x = QQ(x)
  r = QQ(r)
  return ceil(x / r) * r


def _evaluate_matrix(T, a):
  r"""
  a
  """
  M = [[0,0,0],[0,0,0],[0,0,0]]
  for i in range(3):
    for j in range(3):
      if T[i][j].degree() > 0:
        M[i][j] = T[i][j](a)
      else:
        M[i][j] = T[i][j]
  return matrix(a.parent(), M)

