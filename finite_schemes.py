from sage.all import GF, lcm


class FiniteScheme:
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

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x + y])
      sage: X = FiniteScheme(J); X
      Finite Scheme V₊(x^3 + x*z^2 + z^3, x + y) over Rational Field
    """
    if defining_ideal.dimension() - 1 > 0:
      raise ValueError(f"{defining_ideal} does not define a projective scheme of dimension 0.")
    self._defining_ideal = defining_ideal


  def __repr__(self):
    return f"Finite Scheme V\u208A{tuple(self.defining_ideal().gens())} over {self.base_ring()}"


  def defining_ideal(self):
    return self._defining_ideal


  def base_ring(self):
    return self.defining_ideal().base_ring()


  def base_change(self, L):
    r"""
    Return the base change of `self` to `L`.

    EXAMPLES::
      sage: R.<x,y,z> = GF(2)[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x + y])
      sage: X = FiniteScheme(J); X
      Finite Scheme V₊(x^3 + x*z^2 + z^3, x + y) over Finite Field of size 2
      sage: X.base_change(GF(2^3))
      Finite Scheme V₊(x^3 + x*z^2 + z^3, x + y) over Finite Field in z3 of size 2^3
    """
    R = self.defining_ideal().ring()
    R_L = R.change_ring(L)
    J_L = self.defining_ideal().change_ring(R_L)
    return FiniteScheme(J_L)


  def closed_points(self, defining_ideals=True):
    r"""
    Return the closed points of `self`.

    EXAMPLES::
      sage: R.<x,y,z> = QQ[]
      sage: J = R.ideal([x^2 + y*z + z^2, y^2 + y*z + z^2])
      sage: X = FiniteScheme(J)
      sage: X.closed_points()
      [Ideal (y^2 + y*z + z^2, x + y) of Multivariate Polynomial Ring in x, y, z over Rational Field,
      Ideal (y^2 + y*z + z^2, x - y) of Multivariate Polynomial Ring in x, y, z over Rational Field]
      sage: X.closed_points(defining_ideals=False)
      [Finite Scheme V₊(y^2 + y*z + z^2, x - y) over Rational Field,
      Finite Scheme V₊(y^2 + y*z + z^2, x + y) over Rational Field]

      sage: R.<x,y,z> = GF(2)[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x*z^2 + y*z^2])
      sage: X = FiniteScheme(J)
      sage: X.closed_points(defining_ideals=False)
      [Finite Scheme V₊(z, x) over Finite Field of size 2,
      Finite Scheme V₊(y^3 + y*z^2 + z^3, x + y) over Finite Field of size 2]
    """
    ideals = self.defining_ideal().minimal_associated_primes()
    irr_id = self.defining_ideal().ring().irrelevant_ideal()
    ideals = [J for J in ideals if J != irr_id]
    if defining_ideals:
      return ideals
    return [FiniteScheme(I) for I in ideals]


  def splitting_field(self):
    r"""
    Return the minimal extension of the base field
    of `self` making all closed points rational.

    EXAMPLES::
      sage: R.<x,y,z> = GF(2)[]
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x*z^2 + y*z^2])
      sage: X = FiniteScheme(J)
      sage: L = X.splitting_field(); L
      Finite Field in z3 of size 2^3
      sage: X_L = X.base_change(L)
      sage: X_L.closed_points(defining_ideals=False)
      [Finite Scheme V₊(z, x) over Finite Field in z3 of size 2^3,
      Finite Scheme V₊(y + (z3^2 + z3)*z, x + (z3^2 + z3)*z) over Finite Field in z3 of size 2^3,
      Finite Scheme V₊(y + (z3^2)*z, x + (z3^2)*z) over Finite Field in z3 of size 2^3,
      Finite Scheme V₊(y + z3*z, x + z3*z) over Finite Field in z3 of size 2^3]
      sage: X.closed_points(defining_ideals=False)
      [Finite Scheme V₊(z, x) over Finite Field of size 2,
      Finite Scheme V₊(y^3 + y*z^2 + z^3, x + y) over Finite Field of size 2]
      sage:
      sage: J = R.ideal([x^3 + x*z^2 + z^3, x + y])
      sage: X = FiniteScheme(J)
      sage: L = X.splitting_field(); L
      Finite Field in z3 of size 2^3
      sage: X_L = X.base_change(L)
      sage: X_L.closed_points(defining_ideals=False)
      [Finite Scheme V₊(x + y, (z3 + 1)*y + z) over Finite Field in z3 of size 2^3,
      Finite Scheme V₊(x + y, (z3^2 + 1)*y + z) over Finite Field in z3 of size 2^3,
      Finite Scheme V₊(x + y, (z3^2 + z3 + 1)*y + z) over Finite Field in z3 of size 2^3]
      sage: X.closed_points(defining_ideals=False)
      [Finite Scheme V₊(y^3 + y*z^2 + z^3, x + y) over Finite Field of size 2]
      sage:
      sage: J = R.ideal([x^2 + y*z + z^2, y^2 + y*z + z^2])
      sage: X = FiniteScheme(J)
      sage: L = X.splitting_field(); L
      Finite Field in z2 of size 2^2
      sage: X_L = X.base_change(L)
      sage: X_L.closed_points(defining_ideals=False)
      [Finite Scheme V₊(y + (z2 + 1)*z, x + (z2 + 1)*z) over Finite Field in z2 of size 2^2,
      Finite Scheme V₊(y + z2*z, x + z2*z) over Finite Field in z2 of size 2^2]
      sage: X.closed_points(defining_ideals=False)
      [Finite Scheme V₊(y^2 + y*z + z^2, x + y) over Finite Field of size 2]
    """
    if not (self.base_ring().is_field() and self.base_ring().is_finite()):
      raise NotImplementedError(f"{self.base_ring()} is not a finite field.")

    residue_field_extension_degrees = []
    for P in self.closed_points(defining_ideals=True):
      f = P.hilbert_polynomial()
      deg = f.leading_coefficient() * f.degree().factorial()
      residue_field_extension_degrees.append(deg)
    if not residue_field_extension_degrees:
      return self.base_ring()
    splitting_field_ext_deg = lcm(residue_field_extension_degrees)
    q = self.base_ring().order()
    return GF(q**splitting_field_ext_deg)
