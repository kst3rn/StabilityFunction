from sage.all import PolynomialRing, GF, QQ, ZZ, ceil, matrix, GaussValuation
from StabilityFunction import StabilityFunction
from parametric_optimization import minimum_as_valuative_function

def _tree_search():
  return None

def find_base_ring_extension_of(homogeneous_form, base_ring_valuation, ramification_index):
  r"""
  Try to find an extension of the base ring of `homogeneous_form`.
  """

  if homogeneous_form.base_ring() != base_ring_valuation.domain():
    raise ValueError(f"The base ring of {homogeneous_form} is not {base_ring_valuation.domain()}")
  if homogeneous_form.base_ring() is not QQ:
    raise ValueError(f"The base ring must be {QQ}")
  if base_ring_valuation.residue_field() is not GF(2):
    raise ValueError(f"The residue field of {base_ring_valuation} is not {GF(2)}")

  K = homogeneous_form.base_ring()
  phi = StabilityFunction(homogeneous_form, base_ring_valuation)
  minimum, btb_point = phi.global_minimum()
  if phi.has_semistable_reduction_at(btb_point):
    if btb_point.minimal_simplex_dimension() != 0:
      L = K
      # make `L` to an extension of K such that b becomes integral
      return L
    return K

  R = homogeneous_form.parent()
  S = PolynomialRing(K, 's')
  s = S.gen()
  R_S = R.change_ring(S)
  F_S = R_S(homogeneous_form)
  T = btb_point.base_change_matrix()

  # Check totally ramified extensions of degree 2.
  v0 = GaussValuation(S, base_ring_valuation)
  step = ZZ(1)/ZZ(ramification_index)
  fixed_valuation = v0.augmentation(s, step)
  # create here initialization:
  # (1) Determine what shape of base change matrix is needed (get normalized reduction polynomial and check for all elementary instabilities, if btb_point.minimal_simplex_dimension == 1 and otherwise for any instability (maybe no 'diagonal').)


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
  discoids = [b.discoid() for b in biggest_local_maxima]
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
    phi_typeI = StabilityFunction(F_new, new_extension.valuation(2))
    aI, bI = phi_typeI.local_minimum(_evaluate_matrix(global_trafo_matrix, piK))
    if phi_typeI.has_semistable_reduction_at(bI):
      return new_extension

    new_typeII_valuation = fixed_valuation.augmentation(center, new_radius)
    phi_typeII = StabilityFunction(F, new_typeII_valuation)
    new_minimum, new_btb_point = phi.local_minimum(global_trafo_matrix)
    if new_minimum >= minimum:
      break

    j = j + 1
    new_radius = radius - j * step

    w_normalized = [QQ(x / step) for x in new_btb_point.weight_vector()]
    local_trafo_matrix = [[1,0,0],[0,1,0],[0,0,1]]
    for i in range(3):
      for j in range(i + 1, 3):
        w_difference = w_normalized[j] - w_normalized[i]
        if w_difference in ZZ:
          if w_difference >= 0:
            local_trafo_matrix[i][j] = F.parent().gen()**(w_difference)
          else:
            local_trafo_matrix[j][i] = F.parent().gen()**(-w_difference)
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

