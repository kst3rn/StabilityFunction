from sage.all import GF, lcm


class FiniteProjectiveScheme:
  r"""
  Construct a finite scheme to the following conditions.

  INPUT:
  - ``defining_ideal`` -- homogeneous ideal in K[x_0, x_1, x_2].
  """

  def __init__(self, defining_ideal):
    r"""
    Construct a finite scheme to the following conditions.

    INPUT:
    - ``defining_ideal`` -- homogeneous ideal in K[x_0, x_1, x_2].
    """
    if not defining_ideal.dimension() - 1 == 0:
      raise ValueError(f"{defining_ideal} does not define a projective scheme of dimension 0.")
    self._defining_ideal = defining_ideal


  def __repr__(self):
    return f"Finite Projective Scheme with defining ideal {self.defining_ideal()}"


  def defining_ideal(self):
    return self._defining_ideal


  def base_ring(self):
    return self.defining_ideal().base_ring()


  def base_change(self, L):
    r"""
    Return the base change of `self` to `L`.

    EXAMPLES::
      sage: R.<x,y,z> = GF(2)[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x*z^2 + y*z^2])
      sage: X = FiniteProjectiveScheme(J)
      sage: 
      sage: R.<x,y,z> = GF(2)[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x*z^2 + y*z^2])
      sage: X = FiniteProjectiveScheme(J); X
      Finite Projective Scheme with defining ideal Ideal (x^3 + x*z^2 + z^3, x*z^2 + y*z^2) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 2
      sage: X.base_change(GF(2^3))
      Finite Projective Scheme with defining ideal Ideal (x^3 + x*z^2 + z^3, x*z^2 + y*z^2) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3
    """
    R = self.defining_ideal().ring()
    R_L = R.change_ring(L)
    J_L = self.defining_ideal().change_ring(R_L)
    return FiniteProjectiveScheme(J_L)


  def closed_points(self, defining_ideals=True):
    r"""
    Return the closed points of `self`.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: J = R.ideal([x^2 + y*z + z^2, y^2 + y*z + z^2])
      sage: X = FiniteProjectiveScheme(J)
      sage: X.closed_points()
      [Ideal (y^2 + y*z + z^2, x - y) of Multivariate Polynomial Ring in x, y, z over Rational Field,
      Ideal (y^2 + y*z + z^2, x + y) of Multivariate Polynomial Ring in x, y, z over Rational Field]
      sage: X.closed_points(defining_ideals=False)
      [Finite Projective Scheme with defining ideal Ideal (y^2 + y*z + z^2, x - y) of Multivariate Polynomial Ring in x, y, z over Rational Field,
      Finite Projective Scheme with defining ideal Ideal (y^2 + y*z + z^2, x + y) of Multivariate Polynomial Ring in x, y, z over Rational Field]

      sage: R.<x,y,z> = GF(2)[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x*z^2 + y*z^2])
      sage: X = FiniteProjectiveScheme(J)
      sage: X.closed_points()
      [Ideal (z, x) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 2,
      Ideal (y^3 + y*z^2 + z^3, x + y) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 2]
    """
    ideals = self.defining_ideal().minimal_associated_primes()
    if defining_ideals:
      return ideals
    return [FiniteProjectiveScheme(I) for I in ideals]


  def splitting_field(self):
    r"""
    Return the minimal extension of the base field
    of `self` making all closed points rational.

    EXAMPLES::
      sage: R.<x,y,z> = GF(2)[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x*z^2 + y*z^2])
      sage: X = FiniteProjectiveScheme(J)
      sage: L = X.splitting_field(); L
      Finite Field in z3 of size 2^3
      sage: X_L = X.base_change(L)
      sage: X_L.closed_points()
      [Ideal (z, x) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3,
      Ideal (y + (z3^2 + z3)*z, x + (z3^2 + z3)*z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3,
      Ideal (y + (z3^2)*z, x + (z3^2)*z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3,
      Ideal (y + z3*z, x + z3*z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3]
      sage: X.closed_points()
      [Ideal (z, x) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 2,
      Ideal (y^3 + y*z^2 + z^3, x + y) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 2]
      sage:
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x + y])
      sage: X = FiniteProjectiveScheme(J)
      sage: L = X.splitting_field(); L
      Finite Field in z3 of size 2^3
      sage: X_L = X.base_change(L)
      sage: X_L.closed_points()
      [Ideal (x + y, (z3 + 1)*y + z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3,
      Ideal (x + y, (z3^2 + 1)*y + z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3,
      Ideal (x + y, (z3^2 + z3 + 1)*y + z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z3 of size 2^3]
      sage: X.closed_points()
      [Ideal (y^3 + y*z^2 + z^3, x + y) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 2]
      sage:
      sage: J = R.ideal([x^2 + y*z + z^2, y^2 + y*z + z^2])
      sage: X = FiniteProjectiveScheme(J)
      sage: L = X.splitting_field(); L
      Finite Field in z2 of size 2^2
      sage: X_L = X.base_change(L)
      sage: X_L.closed_points()
      [Ideal (y + (z2 + 1)*z, x + (z2 + 1)*z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z2 of size 2^2,
      Ideal (y + z2*z, x + z2*z) of Multivariate Polynomial Ring in x, y, z over Finite Field in z2 of size 2^2]
      sage: X.closed_points()
      [Ideal (y^2 + y*z + z^2, x + y) of Multivariate Polynomial Ring in x, y, z over Finite Field of size 2]
    """
    if not (self.base_ring().is_field() and self.base_ring().is_finite()):
      raise NotImplementedError(f"{self.base_ring()} is not a finite field.")

    residue_field_degrees = []
    for P in self.closed_points(defining_ideals=True):
      f = P.hilbert_polynomial()
      deg = f.leading_coefficient() * f.degree().factorial()
      residue_field_degrees.append(deg)
    if not residue_field_degrees:
      return self.base_ring()
    splitting_field_degree = lcm(residue_field_degrees)
    q = self.base_ring().order()
    return GF(q**splitting_field_degree)
