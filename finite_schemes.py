from sage.all import GF, lcm


def _rationalizator_over_GF(J):
  r"""
  Return the minimal extension of the base field of `J` making
  all closed points of V(J) rational.

  INPUT:
  - ``J`` -- an ideal of dimension zero in K[x,y].

  OUTPUT:
  A finite field extension of K making all closed points
  of the variety defined by J rational.

  EXAMPLES::
    sage: R.<x,y> = GF(2)[]
    sage: J = R.ideal(x^2 + x + 1, y^2 + y + 1)
    sage: L = _rationalizator_over_GF(J); L
    Finite Field in z2 of size 2^2
    sage: RL = R.change_ring(L)
    sage: JL = J.change_ring(RL)
    sage: J.minimal_associated_primes()
    [Ideal (y^2 + y + 1, x + y) of Multivariate Polynomial Ring in x, y over Finite Field of size 2,
    Ideal (y^2 + y + 1, x + y + 1) of Multivariate Polynomial Ring in x, y over Finite Field of size 2]
    sage: JL.minimal_associated_primes()
    [Ideal (y + (z2 + 1), x + z2) of Multivariate Polynomial Ring in x, y over Finite Field in z2 of size 2^2,
    Ideal (y + (z2 + 1), x + (z2 + 1)) of Multivariate Polynomial Ring in x, y over Finite Field in z2 of size 2^2,
    Ideal (y + z2, x + z2) of Multivariate Polynomial Ring in x, y over Finite Field in z2 of size 2^2,
    Ideal (y + z2, x + (z2 + 1)) of Multivariate Polynomial Ring in x, y over Finite Field in z2 of size 2^2]
    sage: 
    sage: J = R.ideal(x^3 + x + 1, y + x)
    sage: L = _rationalizator_over_GF(J); L
    Finite Field in z3 of size 2^3
    sage: RL = R.change_ring(L)
    sage: JL = J.change_ring(RL)
    sage: J.minimal_associated_primes()
    [Ideal (y^3 + y + 1, x + y) of Multivariate Polynomial Ring in x, y over Finite Field of size 2]
    sage: JL.minimal_associated_primes()
    [Ideal (x + y, y + (z3^2 + z3)) of Multivariate Polynomial Ring in x, y over Finite Field in z3 of size 2^3,
    Ideal (x + y, y + z3) of Multivariate Polynomial Ring in x, y over Finite Field in z3 of size 2^3,
    Ideal (x + y, y + (z3^2)) of Multivariate Polynomial Ring in x, y over Finite Field in z3 of size 2^3]
  """
  base_ring = J.base_ring()
  if not (base_ring.is_field() and base_ring.is_finite()):
    raise ValueError(f"{base_ring} is not a finite field.")

  primes = J.minimal_associated_primes()
  degrees = [P.vector_space_dimension() for P in primes]
  if not degrees:
    return base_ring
  L_degree = lcm(degrees)
  q = base_ring.order()
  return GF(q**L_degree)
