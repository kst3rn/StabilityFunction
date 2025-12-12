from itertools import combinations
from sage.all import PolynomialRing, GF, QQ, ZZ, ceil, matrix, GaussValuation, vector, Infinity
from stability_function import StabilityFunction
from plane_curves import ProjectivePlaneCurve
from parametric_optimization import minimum_as_valuative_function

def _tree_search():
  return None

def find_base_ring_extension(homogeneous_form, base_ring_valuation, ramification_index):
  r"""
  Try to find an extension of the base ring of `homogeneous_form`.

  EXAMPLES::
    sage: R.<x0,x1,x2> = QQ[]
    sage: F = x1^4 + 2*x0^3*x2 + x0*x1^2*x2 + 2*x0*x2^3
    sage: find_base_ring_extension(F, QQ.valuation(2), 2)
    Number Field in piK with defining polynomial s^2 + 2

    sage: R.<x,y,z> = QQ[]
    sage: R.<x,y,z> = QQ[]
    sage: F = 16*x**4 + y**4 + 8*y**3*z + 16*x*y*z**2 + 4*x*z**3
    sage: find_base_ring_extension(F, QQ.valuation(2), 4)
    'Any extension of Number Field in piK with defining polynomial s^4 + 2*s^2 + 2 making the point [0, 3/2, 4/3] integral'
    sage: F = 4*x**4 + 4*x*y**3 + y**4 + 2*x*z**3 + 4*y*z**3 + z**4
    sage: find_base_ring_extension(F, QQ.valuation(2), 4)
    Number Field in piK with defining polynomial s^4 + 2*s^3 + 2*s^2 + 2

    sage: R.<x0,x1,x2> = QQ[]
    sage: G = -2*x0^3*x1 - 12*x1^4 - 4*x0^3*x2 - 3*x0^2*x1*x2 - 12*x1^3*x2 - 4*x0^2*x2^2 - 12*x0*x1*x2^2 + 16*x1^2*x2^2 + 5*x1*x2^3
    sage: find_base_ring_extension(G, QQ.valuation(2), 4)
    'Any extension of Number Field in piK with defining polynomial s^4 + 2 making the point [0, 3/8, 25/16] integral'
  """

  if homogeneous_form.base_ring() != base_ring_valuation.domain():
    raise ValueError(f"The base ring of {homogeneous_form} is not {base_ring_valuation.domain()}")
  if homogeneous_form.base_ring() is not QQ:
    raise ValueError(f"The base ring must be {QQ}")
  if base_ring_valuation.residue_field() is not GF(2):
    raise ValueError(f"The residue field of {base_ring_valuation} is not {GF(2)}")

  K = homogeneous_form.base_ring()
  phi = StabilityFunction(homogeneous_form, base_ring_valuation)
  minimum, btb_point = phi.global_minimum('uut')
  if phi.has_semistable_reduction_at(btb_point):
    if btb_point.minimal_simplex_dimension() != 0:
      r = btb_point.ramification_index()
      L = K
      # make `L` to an extension of K such that b becomes integral
      return f"Any extension of {K} making the point {btb_point.weight_vector()} integral"
    return K

  R = homogeneous_form.parent()
  S = PolynomialRing(K, 's')
  s = S.gen()
  R_S = R.change_ring(S)
  F_S = R_S(homogeneous_form)
  T = btb_point.base_change_matrix()

  v0 = GaussValuation(S, base_ring_valuation)
  step = ZZ(1)/ZZ(ramification_index)
  fixed_valuation = v0.augmentation(s, step)

  # create here initialization matrix:
  w_normalized = [QQ(x / step) for x in btb_point.weight_vector()]
  F_b = phi.graded_reduction(btb_point)
  f = F_b.normalized_reduction_polynomial()
  X_b = ProjectivePlaneCurve(f)
  local_trafo_matrix = [[1,0,0],[0,1,0],[0,0,1]]
  # combinations(range(3), 2) yields (0,1), (0,2), (1,2)
  for i, j in combinations(range(3), 2):
    w_difference = w_normalized[j] - w_normalized[i]
    if w_difference in ZZ and X_b.elementary_instability_direction((i,j)) is not None:
      if w_difference >= 0:
        local_trafo_matrix[i][j] = s**w_difference
      else:
        local_trafo_matrix[j][i] = s**(-w_difference)
      break        
  local_trafo_matrix = matrix(S, local_trafo_matrix)
  global_trafo_matrix = local_trafo_matrix * T

  return _search_tree(F_S, fixed_valuation, step, minimum, global_trafo_matrix, 0, depth_limit=+Infinity)


def _search_tree(F, fixed_valuation, step, minimum, global_trafo_matrix, depth, depth_limit): # test with F = 16*x**4 + y**4 + 8*y**3*z + 16*x*y*z**2 + 4*x*z**3
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
  while True:
    K = QQ.extension(center, 'piK')
    piK = K.gen()
    R_K = F.parent().change_ring(K)
    F_K = R_K(F)
    phi_typeI = StabilityFunction(F_K, K.valuation(2))
    aI, bI = phi_typeI.local_minimum(_evaluate_matrix(global_trafo_matrix, piK))
    if phi_typeI.has_semistable_reduction_at(bI):
      if bI.minimal_simplex_dimension(ZZ(1) / step) != 0:
        return f"Any extension of {K} making the point {bI.weight_vector()} integral"
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

