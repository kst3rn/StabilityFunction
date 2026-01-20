from itertools import permutations
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class BruhatTitsBuilding:
  r"""
  Construct a Bruhat-Tits building.

  INPUT:
  - ``base_ring_valuation`` -- discrete valuation on a field.
  - ``dimension`` -- positive integer.
  """

  def __init__(self, base_ring_valuation, dimension):
    r"""
    Construct a Bruhat-Tits building.

    INPUT:
    - ``base_ring_valuation`` -- discrete valuation on a field.
    - ``dimension`` -- positive integer.
    """
    self._base_ring_valuation = base_ring_valuation
    self._dimension = dimension


  def base_ring_valuation(self):
    r"""
    Return the base ring valuation of `self`.
    """
    return self._base_ring_valuation


  def dimension(self, ):
    r"""
    Return the dimension of `self`.
    """
    return self._dimension


  def dim_plus1(self):
    r"""
    Return the dimension plus one.
    """
    return self.dimension() + 1



class Apartment(BruhatTitsBuilding):
  r"""
  Construct an apartment.

  INPUT:
  - ``base_ring_valuation`` -- discrete valuation on a field.
  - ``base_change_matrix`` -- invertible matrix.
  """
  def __init__(self, base_ring_valuation, base_change_matrix):
    r"""
    Construct an apartment.

    INPUT:
    - ``base_ring_valuation`` -- discrete valuation on a field.
    - ``base_change_matrix`` -- invertible matrix.
    """
    if not base_change_matrix.is_invertible():
      raise ValueError("The base change matrix must be invertible.")
    dimension = base_change_matrix.nrows() - 1
    super().__init__(base_ring_valuation, dimension)
    self._base_change_matrix = base_change_matrix


  def base_change_matrix(self):
    r"""
    Return the base change matrix of `self`.
    """
    return self._base_change_matrix


  def intersection(self, apartment, complement=False, affine_patch=None):
    r"""
    Print the intersection of `self` with `apartment`
    on `affine_patch`. If `complement` is `True`, print
    the complement of the intersection.

    EXAMPLES::
      sage: L.<a> = QQ.extension(x^8 + 8*x^6 + 120*x^4 + 32*x^2 + 16)
      sage: v_L = L.valuation(2)
      sage: pi1 = 1/16*a^7 + 1/2*a^5 + 15/2*a^3 + 5/2*a
      sage: pi2 = -1/16*a^7 - 1/2*a^5 - 15/2*a^3 - 3/2*a
      sage: T1 = matrix(L, [[1,0,0],[-1/pi2,1,0],[-pi2^2,0,1]])
      sage: M1 = matrix(L, [[1,0,0],[-1/pi1,1,0],[-pi1^2,0,1]])
      sage: A1 = Apartment(v_L, T1)
      sage: A2 = Apartment(v_L, M1)
      sage: A1.intersection(A2)
      Intersection of the following conditions:
      w0 <= w1
      w0 <= w2 + 1
      sage: A1.intersection(A2, complement=True, affine_patch=2)
      Union of the following conditions:
      w0 > w1
      w0 > 1
    """
    if not isinstance(complement, bool):
      raise ValueError(f"{complement} is not a boolean.")
    if affine_patch is not None:
      if not isinstance(affine_patch, (int, Integer)):
        raise ValueError(f"{affine_patch} is not an integer.")
      elif not 0 <= affine_patch <= self.dimension():
        raise ValueError(
          f"{affine_patch} is not between 0 and {self.dimension()}"
          )

    v_K = self.base_ring_valuation()
    T1 = self.base_change_matrix()
    T2 = apartment.base_change_matrix()
    T = T1.inverse() * T2

    if any(v_K(T[i][i]) for i in range(self.dim_plus1())):
      raise NotImplementedError(
        "Provided base change matrix combination is not allowed. "
        "There is a diagonal element of T1^{-1} * T2 with nonzero valuation."
        )

    R = PolynomialRing(QQ, self.dim_plus1(), 'w')
    w = list(R.gens())
    if affine_patch is not None:
      w[affine_patch] = R(0)

    if complement == False:
      print("Intersection of the following conditions:")
    else:
      print("Union of the following conditions:")
    for i, j in permutations(range(self.dim_plus1()), 2):
      val = v_K(T[i][j])
      if val == +Infinity:
        continue
      if complement == False:
        print(f"{w[j]} <= {v_K(T[i][j]) + w[i]}")
      else:
        print(f"{w[j]} > {v_K(T[i][j]) + w[i]}")
