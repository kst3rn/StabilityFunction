

class BruhatTitsBuilding:
  r"""
  Construct a Bruhat-Tits building.
  """
  def __init__(self, base_ring_valuation, dimension):
    r"""
    Construct a Bruhat-Tits building.
    """
    self._base_ring_valuation = base_ring_valuation
    self._dimension = dimension



class Apartment:
  r"""
  Construct an apartment.
  """
  def __init__(self, base_ring_valuation, base_change_matrix):
    r"""
    Construct an apartment.
    """
    self._base_ring_valuation = base_ring_valuation
    self._base_change_matrix = base_change_matrix
    self._dimension = base_change_matrix.nrows()


  def base_ring_valuation(self):
    r"""
    Return the base ring valuation of `self`.
    """
    return self._base_ring_valuation


  def base_change_matrix(self):
    r"""
    Return the base change matrix of `self`.
    """
    return self._base_change_matrix


  def dimension(self):
    r"""
    Return the dimension of `self`.
    """
    return self._dimension


  def intersection(self, apartment, affine_patch=None):
    r"""
    Return the intersection of `self` with `apartment`
    on `affine_patch`.
    """
    v_K = self.base_ring_valuation()
    T1 = self.base_change_matrix()
    T2 = apartment.base_change_matrix()
    T = T1.inverse() * T2
    if any(v_K(T[i][i]) for i in range(self.dimension())):
      raise NotImplementedError("Provided base change matrix combination is not allowed.")
