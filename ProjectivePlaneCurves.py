
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from functools import cached_property
from sage.all import *


class ProjectivePlaneCurve:

    def __init__(self, polynomial):
        """
        INPUT:
            ...
        """

        if not polynomial.is_homogeneous():
            raise TypeError

        if not len(polynomial.parent().gens()) == 3:
            raise ValueError

        if not polynomial.base_ring().is_field():
            raise ValueError

        self.polynomial = polynomial
        self.degree = self.polynomial.degree()
        self.polynomial_ring = self.polynomial.parent()
        self.base_ring = self.polynomial.base_ring()
        self.projective_plane = ProjectiveSpace(self.polynomial_ring)
        self.plane_curve = self.projective_plane.curve(self.polynomial)
        self.standard_basis = self.polynomial_ring.gens()


    def __repr__(self):
        return f"Projective Plane Curve with defining polynomial {self.polynomial}"


    def get_base_ring(self):
        return self.base_ring


    def get_polynomial(self):
        return self.polynomial


    def get_standard_basis(self):
        return self.standard_basis


    def tangent_cone_at(self, P):
        """
        Return the tangent cone of self at P
        """
        return PPC_TangentCone(self, P)


    def is_smooth(self):
        return self.plane_curve.is_smooth()


    def is_reduced(self):
        return not any(multiplicity > 1 for factor, multiplicity in self._decompose)


    def is_irreducible(self):
        if len(self._decompose) > 1:
            return False
        return self.is_reduced()


    def is_semistable(self):
        if len(self.instabilities()) == 0:
            return True
        return False


    def is_stable(self):
        r"""
        ToDo
        """


    def reduced_subscheme(self):
        r"""
        Return the reduced subscheme of self as a projective plane curve
        """
        f = prod(factor for factor, multiplicity in self._decompose)
        return ProjectivePlaneCurve(f)


    def irreducible_components(self):
        return [ProjectivePlaneCurve(factor)
                for factor, multiplicity in self._decompose]


    def nonreduced_components(self):
        return [(multiplicity, ProjectivePlaneCurve(factor))
                for factor, multiplicity in self._decompose
                if multiplicity > 1]


    def points_with_high_multiplicity(self):
        """
        Return a list of points on self with multiplicity > self.degree/2

        OUTPUT:
            [(point_1, m_1), (point_2, m_2), ...] - where m_i is the
                                                    multiplicity of point_i
        """

        L = []
        for P in self.plane_curve.singular_points():
            m = self.plane_curve.multiplicity(P)
            if m > self.degree / Integer(2):
                L.append((P, m))
        return L


    def lines_with_high_multiplicity(self):
        """
        Return a list of lines in self of multiplicity > 0

        OUTPUT:
            [(line_1, m_1, G_1), ...] - where m_i is the multiplicity of line_i,
                                        and either line_i^m_i * G_i = self.polynomial
                                        if 0 < m_i <= self.degree/3 or G_i = None else.
        """

        L = []
        polynomial_factors = list(self.polynomial.factor())
        for i, (factor, factor_multiplicity) in enumerate(polynomial_factors):
            if factor.degree() > 1:
                continue

            if factor_multiplicity > self.degree / Integer(3):
                L.append((factor, factor_multiplicity, None))

            else:
                G = Integer(1)
                for j, (Gfactor, Gfactor_multiplicity) in enumerate(polynomial_factors):
                    if j != i:
                        G = G * Gfactor**Gfactor_multiplicity
                L.append((factor, factor_multiplicity, G))

        return L


    def pseudo_instabilities(self):
        """
        Return the list of all pseudo-instabilities.

        MATHEMATICAL INTERPRETATION:
            We will explain how to compute all such pseudo-instabilities and why there are
            only finitely many of them.
        """

        list_of_pseu_inst = []

        # CASES (b) and (d)
        for P, m in self.points_with_high_multiplicity():
            if m > 2 * self.degree / 3:  # CASE (b)
                list_of_pseu_inst.append(PseudoInstability(self, Point = P))
            elif m > self.degree / 2:  # CASE (d) 1/2
                for L, L_multiplicity in PPC_TangentCone(self, P).embedded_lines():
                    if L_multiplicity > m / 2:
                        list_of_pseu_inst.append(PseudoInstability(self, P, L))

        # CASES (a) and (c)
        for L, L_multiplicity, G in self.lines_with_high_multiplicity():
            if L_multiplicity > self.degree / Integer(3):  # CASE (a)
                list_of_pseu_inst.append(PseudoInstability(self, Line = L))
            else: # CASE (c) 1/2
                L_curve = self.projective_plane.curve(L)
                G_curve = self.projective_plane.curve(G)
                for P in L_curve.intersection_points(G_curve):
                    if L_curve.intersection_multiplicity(G_curve, P) > (self.degree - L_multiplicity) / Integer(2):
                        list_of_pseu_inst.append(PseudoInstability(self, P, L))

        return list_of_pseu_inst


    def instabilities(self):
        """
        Return the list of all pseudo-instabilities which yield instabilities
        """

        return [I for I in self.pseudo_instabilities() if I.is_instability()]


    @cached_property
    def _decompose(self):
        r"""
        Return the factored form of self.polynomial.
        """
        return list(self.polynomial.factor())



class PPC_TangentCone:

    def __init__(self, projective_plane_curve, Point):
        """
        At the moment only rational points are allowed
        """

        # Convert to list
        Point = list(Point)

        if projective_plane_curve.get_polynomial()(Point) != 0:
            raise ValueError

        self.projective_plane_curve = projective_plane_curve
        self.base_ring = projective_plane_curve.get_base_ring()
        self.polynomial_ring = PolynomialRing(self.base_ring, ['x', 'y'])
        self.gen1, self.gen2 = self.polynomial_ring.gens()

        # Coerce coordinates to self.base_ring
        for i in range(3):
            Point[i] = self.base_ring(Point[i])

        # Find the maximal index, i_max, with Point[i_max] != 0 and normalize by Point[i_max]
        self.affine_patch, self.normalized_point = _normalize_by_last_nonzero_entry(Point)


    def __repr__(self):
        return f"Tangent cone of {self.projective_plane_curve} at {self.normalized_point}"


    def get_polynomial(self):

        PPC_equation = self.projective_plane_curve.get_polynomial()
        dehomogenization = [self.gen1, self.gen2]
        dehomogenization.insert(self.affine_patch, self.polynomial_ring(0))
        dehomogenization_translated = []

        for i in range(3):
            dehomogenization_translated.append(dehomogenization[i] + self.normalized_point[i])

        f = PPC_equation(dehomogenization_translated)

        f_homo_comp_dict = f.homogeneous_components()
        minimal_degree = min(f_homo_comp_dict.keys())
        tangent_cone_polynomial = f_homo_comp_dict[minimal_degree]

        return tangent_cone_polynomial


    def embedded_polynomial(self):
        """
        MATHEMATICAL INTERPRETATION:
            First, let
                F = self.projective_plane_curve.get_polynomial(),
                K = self.projective_plane_curve.get_base_ring(),
                P = self.normalized_point,
                j = self.affine_patch,
                x_0, x_1, x_2 = self.projective_plane_curve.get_standard_basis().
            Then P[j] = 1 and for all j < i <= 2 we have
                P[i] = 0.
            Further, we set y_0 = x_0, y_1 = x_1, y_2 = x_2, y_j = 1 and 
                f = F(y_0 + P[0]*y_j, y_1 + P[1]*y_j, y_2 + P[2]*y_j).
            This is the dehomogenization of F to the affine patch j and the
            subsequent translation of the point P to the origin in this affine
            patch. Thus,
                f = _apply_matrix(T, F, affine_patch = j)
            with
                T = _ult_line_transformation(K, P).
            Now let TanCon be the homogeneous component of f of the lowest degree,
            i.e. TanCon is the homogeneous polynomial defining the tangent cone of
            self.projective_plane_curve at the point P. We view TanCon as a polynomial
            in K[x_0, x_1, x_2], although it does not depend on x_j. Let e_j be the
            j-th standard basis vector. Then we have
                TanCon(e_j) = 0 and e_j*T = P.
            Now we define
                TanCon_P = TanCon(x_0 - P[0]*x_j, x_1 - P[1]*x_j, x_2 - P[2]*x_j),
            i.e.
                TanCon_P = _apply_matrix(T.inverse(), TanCan).
            In particular,
                TanCon_P(P) = TanCon(e_j) = 0.
            Moreover, the dehomogenization of TanCon_P to the affine patch j is the
            translation of TanCon from the origin to the point P.            
        """

        T = _ult_line_transformation(self.base_ring, self.normalized_point)
        T_inverse = T.inverse()
        F = self.projective_plane_curve.get_polynomial()
        f = _apply_matrix(T, F, self.affine_patch)

        f_homo_comp_dict = f.homogeneous_components()
        minimal_degree = min(f_homo_comp_dict.keys())
        tangent_cone_polynomial = f_homo_comp_dict[minimal_degree]

        return _apply_matrix(T.inverse(), tangent_cone_polynomial)


    def get_lines(self):
        """
        Return linear factors of the polynomial self.get_polynomial()
        """

        L = []
        tangent_cone_polynomial = self.get_polynomial()
        for factor, factor_multiplicity in list(tangent_cone_polynomial.factor()):
            if factor.degree() == 1:
                L.append((factor, factor_multiplicity))
        return L


    def embedded_lines(self):
        """
        Return linear factors of the polynomial self.embedded_polynomial()
        """

        L = []
        tangent_cone_polynomial = self.embedded_polynomial()
        for factor, factor_multiplicity in list(tangent_cone_polynomial.factor()):
            if factor.degree() == 1:
                L.append((factor, factor_multiplicity))
        return L



class PseudoInstability:

    def __init__(self, projective_plane_curve, Point = None, Line = None):
        if point == None and line == None:
            raise ValueError

        if Point != None:
            Point = list(Point)

        self.proj_plane_curve = projective_plane_curve
        self.point = Point
        self.line = Line
        self.base_ring = self.proj_plane_curve.get_base_ring()

    def __repr__(self):
        if self.point == None:
            return f"Pseudo-instability of {self.proj_plane_curve} given by {self.line}"
        elif self.line == None:
            return f"Pseudo-instability of {self.proj_plane_curve} given by {self.point}"
        else:
            return f"Pseudo-instability of {self.proj_plane_curve} given by {self.point} and {self.line}"


    def flag(self):
        """
        Return the flag defining self
        """
        if self.point == None:
            return self.line
        elif self.line == None:
            return self.point
        return (self.point, self.line)


    def point_transformation(self, matrix_form = 'uut'):
        """
        Return unipotent matrix transforming a standard basis vector to self.point
        """

        if self.point == None:
            return identity_matrix(self.base_ring, 3)

        if matrix_form == 'uut':
            return _uut_line_transformation(self.base_ring, self.point)
        elif matrix_form == 'ult':
            return _ult_line_transformation(self.base_ring, self.point)
        elif isinstance(matrix_form, list):
            return _integral_line_transformation(self.base_ring, self.point, matrix_form)
        else:
            raise ValueError


    def base_change_matrices(self, matrix_form = 'uut'):
        """
        Return unipotent matrix transforming a standard basis vector and the line
        spanned by it to self.flag()
        """

        if self.point == None:
            if matrix_form == 'uut':
                return _uut_plane_transformation(self.line)
            elif matrix_form == 'ult':
                return _ult_plane_transformation(self.line)
            elif isinstance(matrix_form, list):
                return _integral_plane_transformation(self.line, matrix_form)
            else:
                raise ValueError
        elif self.line == None:
            return [self.point_transformation(matrix_form)]
        else:
            if matrix_form == 'uut':
                return [_uut_flag_transformation(self.point, self.line)]
            elif matrix_form == 'ult':
                return [_ult_flag_transformation(self.point, self.line)]
            elif isinstance(matrix_form, list):
                return [_integral_flag_transformation(self.point, self.line, matrix_form)]
            else:
                raise ValueError


    def get_base_change_matrix(self, matrix_form = 'uut'):
        return self.base_change_matrices(matrix_form)[0]


    def is_instability(self):
        """
        Return True or False depending on whether self corresponds to
        an instability of self.proj_plane_curve or not

        MATHEMATICAL INTERPRETATION:
            First, let
                K = self.base_ring,
                T = self.flag_transformation(),
                F = self.proj_plane_curve.get_polynomial().
            Furthermore, let
                (x0, x1, x2) = self.proj_plane_curve.get_standard_basis()
            and
                G = F((x0,x1,x2)*T),
            i.e.
                G = _apply_matrix(T, F).
            For a multi index set I subset NN^3 we can write
                G = sum_{i in I} a_i x0^i0 * x1^i1 * x2^i2
            with a_i != 0 for all i in I.
            Note that with respect to the new basis (x0,x1,x2)*T the flag
            (self.point, self.line) is given by (e_j, x_i). Thus, it yields
            an instability, if there exists a balanced weight vector
                (w0, w1, w2) in QQ^3, i.e. w0 + w1 + w2 = 0,
            such that
                i0*w0 + i1*w1 + i2*w2 > 0
            for all i in I. 
            REMARK. Any nonzero multiple of a balanced weight vector is again
            a balanced weight vector. Thus, it suffices to consider
                (w0, w1, w2) in QQ^3
            wtih
                -1 <= w0, w1, w2 <= 1.
            Thus, we only have to maximize the function
                min(i0*w0 + i1*w1 + i2*w2 : i in I)
            under the constraints -1 <= w0, w1, w2 <= 1 and to check
            whether the maximum is > 0 or not.
            REMARK. If self.point or self.line is None, then self is an
            instability, see [Proposition 2.6, SternWewers].
        """

        if self.point == None:
            return True

        if self.line == None:
            return True

        T = self.get_base_change_matrix()
        F = self.proj_plane_curve.get_polynomial()
        G = _apply_matrix(T, F)

        MILP = MixedIntegerLinearProgram(solver='PPL')
        v = MILP.new_variable()

        t = v['maximum']
        w0 = v['w0']
        w1 = v['w1']
        w2 = v['w2']

        MILP.set_objective(t)

        MILP.add_constraint(-1 <= w0 <= 1)
        MILP.add_constraint(-1 <= w1 <= 1)
        MILP.add_constraint(-1 <= w1 <= 1)
        MILP.add_constraint(w0 + w1 + w2 == 0)

        for i in G.dict():
            MILP.add_constraint(t <= i[0] * w0 + i[1] * w1 + i[2] * w2)

        MILP.solve()
        values = MILP.get_values(v)

        return values['maximum'] > 0



class PPC_Instability:

    def __init__(self, base_change_matrix, weight_vector, geometric_type = 'not specified'):

        self.base_change_matrix = base_change_matrix
        self.weight_vector = weight_vector
        self.geometric_type = geometric_type


    def get_base_change_matrix(self):
        return self.base_change_matrix



    def get_weight_vector(self):
        return self.weight_vector



    def get_geometric_type(self):
        return self.geometric_type



# ================== helper functions ==================

def _min_index_of_nonzero_entry(L):
    """
    Return minimal integer i between 0 and len(L) with L[i] != 0

    INPUT:
        L - list of rational numbers

    OUTPUT:
        i - minimal integer between 0 and len(L) with L[i] != 0
    """

    return next((i for i in range(len(L)) if L[i] != 0), None)


def _max_index_of_nonzero_entry(L):
    """
    Return maximal integer i between 0 and len(L) with L[i] != 0

    INPUT:
        L - list of rational numbers

    OUTPUT:
        i - maximal integer between 0 and len(L) with L[i] != 0
    """

    return next((i for i in reversed(range(len(L))) if L[i] != 0), None)


def _normalize_by_first_nonzero_entry(L):
    """
    Return ...
    """

    i_min = _min_index_of_nonzero_entry(L)
    first_nonzero_entry = L[i_min]

    return (i_min, [x / first_nonzero_entry for x in L])


def _normalize_by_last_nonzero_entry(L):
    """
    Return ...
    """

    i_max = _max_index_of_nonzero_entry(L)
    last_nonzero_entry = L[i_max]

    return (i_max, [x / last_nonzero_entry for x in L])


def _apply_matrix(T, F, affine_patch = None):
    """
    Return F((x_0,...,x_n) * T) or its dehomogenization
    at affine_patch, i.e. x_{affine_patch} = 1 if
    affine_patch != None

    INPUT:
        T            - matrix over K
        F            - polynomial in K[x_0,...,x_n]
        affine_patch - integer between 0 and n

    OUTPUT:
        F((x_0,...,x_n) * T) with x_{affine_patch} = 1 if
        affine_patch != None

    MATHEMATICAL INTERPRETATION:
        ToDo...
    """

    generators = list(F.parent().gens())
    if affine_patch != None:
        generators[affine_patch] = F.parent()(1)
    return F(list( vector(generators) * T ))


def _ult_line_transformation(base_field, Vector):
    """
    Return unipotent lower triangular matrix over base_field transforming
    a line spanned by a standard basis vector to the line spanned by Vector
    """

    Vector = list(Vector)
    T = identity_matrix(base_field, len(Vector))

    # Find the maximal index, i_max, with Vector[i_max] != 0 and normalize by Vector[i_max]
    i_max, normalized_Vector = _normalize_by_last_nonzero_entry(Vector)

    T[i_max] = normalized_Vector

    return matrix(base_field, T)


def _uut_line_transformation(base_field, Vector):
    """
    Return unipotent upper triangular matrix over base_field transforming
    a line spanned by a standard basis vector to the line spanned by Vector
    """

    Vector = list(Vector)
    T = identity_matrix(base_field, len(Vector))

    # Find the minimal index, i_min, with Vector[i_min] != 0 and normalize by Vector[i_min]
    i_min, normalized_Vector = _normalize_by_first_nonzero_entry(Vector)

    T[i_min] = normalized_Vector

    return matrix(base_field, T)


def _sorting_permutation_matrix(w):
    """
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
    """
    Return unipotent matrix, T, over base_field transforming a line spanned
    by a standard basis vector to the line spanned by Vector such that for
    all indices i,j the implication
        (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
    holds

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
    """
    Return a list of unipotent lower triangular matrices with maximal
    number of zeros which transform the plane defined by linear_form
    to a plane defined by x_i = 0

    MATHEMATICAL INTERPRETATION:
        First, let
            L = linear_form .
        Then we can write
            L = A*x + B*y + C*z .
        Depending on whether
            A != 0 or A != 0 and B != 0 or A = 0 and B != 0 or A = 0 and B = 0
        the unipotent lower triangular matrices with maximal number of zeros, T,
        moving L to a coordinate axis, via L((x, y, z) * T), are given by the four
        matrices:

            [[  1   0, 0],
             [-B/A, 1, 0],
             [-C/A, 0, 1 ]]

            [[  1,    0,  0],
             [-B/A,   1,  0],
             [0,    -C/B, 1]]

            [[1,   0,  0],
             [0,   1,  0],
             [0, -C/B, 1]]

            [[1, 0, 0],
             [0, 1, 0],
             [0, 0, 1]]

        Note that in the case A != 0 and B != 0 both matrices

            [[  1   0, 0],
             [-B/A, 1, 0],
             [-C/A, 0, 1 ]]

        and

            [[  1,    0,  0],
             [-B/A,   1,  0],
             [0,    -C/B, 1]]

        are equal if and only if C = 0.
    """

    base_ring = linear_form.base_ring()

    x0, x1, x2 = linear_form.parent().gens()
    A = linear_form.monomial_coefficient(x0)
    B = linear_form.monomial_coefficient(x1)
    C = linear_form.monomial_coefficient(x2)
    L = []

    if A != 0:
        T = [[1, 0, 0], [-B/A, 1, 0], [-C/A, 0, 1]]
        T = matrix(base_ring, T)
        L.append(T)
        if B != 0 and C != 0:
            T = [[1, 0, 0], [-B/A, 1, 0], [0, -C/B, 1]]
            T = matrix(base_ring, T)
            L.append(T)

    elif B != 0:
        T = [[1, 0, 0], [0, 1, 0], [0, -C/B, 1]]
        T = matrix(base_ring, T)
        L.append(T)

    else:
        T = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        T = matrix(base_ring, T)
        L.append(T)

    return L


def _uut_plane_transformation(linear_form):
    """
    Return a list of unipotent upper triangular matrices with maximal
    number of zeros which transform the plane defined by linear_form
    to a plane defined by x_i = 0

    MATHEMATICAL INTERPRETATION:
        First, let
            L = linear_form .
        Then we can write
            L = A*x + B*y + C*z .
        Depending on whether
            C != 0 or C != 0 and B != 0 or C = 0 and B != 0 or C = 0 and B = 0
        the unipotent upper triangular matrices with maximal number of zeros, T,
        moving L to a coordinate axis, via L((x, y, z) * T), are given by the four
        matrices:

            [[1, 0, -A/C],
             [0, 1, -B/C],
             [0, 0,   1 ]]

            [[1, -A/B,   0],
             [0,   1,  -B/C],
             [0,   0,    1]]

            [[1, -A/B, 0],
             [0,   1,  0],
             [0,   0,  1]]

            [[1, 0, 0],
             [0, 1, 0],
             [0, 0, 1]]

        Note that in the case C != 0 and B != 0 both matrices

            [[1, 0, -A/C],
             [0, 1, -B/C],
             [0, 0,   1 ]]

        and

            [[1, -A/B,  0],
             [0,   1, -B/C],
             [0,   0,   1]]

        are equal if and only if A = 0.
    """

    base_ring = linear_form.base_ring()

    x0, x1, x2 = linear_form.parent().gens()
    A = linear_form.monomial_coefficient(x0)
    B = linear_form.monomial_coefficient(x1)
    C = linear_form.monomial_coefficient(x2)
    L = []

    if C != 0:
        T = [[1, 0, -A/C], [0, 1, -B/C], [0, 0, 1]]
        T = matrix(base_ring, T)
        L.append(T)
        if B != 0 and A != 0:
            T = [[1, -A/B, 0], [0, 1, -B/C], [0, 0, 1]]
            T = matrix(base_ring, T)
            L.append(T)

    elif B != 0:
        coordinate_axis_index = 1
        T = [[1, -A/B, 0], [0, 1, 0], [0, 0, 1]]
        T = matrix(base_ring, T)
        L.append(T)

    else:
        coordinate_axis_index = 0
        T = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        T = matrix(base_ring, T)
        L.append(T)

    return L


def _integral_plane_transformation(linear_form, weight_vector):
    """
    Return unipotent matrix, T, over base_field transforming a plane given
    by x_i = 0 to the plane defined by linear_form = 0 such that for all
    indices i,j the implication
        (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
    holds

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

    # convert all entries to Rationals
    for i, w in enumerate(weight_vector):
        weight_vector[i] = QQ(w)

    P = _sorting_permutation_matrix(weight_vector)
    pL = _apply_matrix(P.transpose(), linear_form)
    list_of_matrices = [P*T*P.transpose() for T in _ult_plane_transformation(pL)]

    return list_of_matrices


def _ult_flag_transformation(Vector, linear_form):
    """
    Return unipotent lower triangular matrix transforming a flag given
    by a line spanned by a standard basis vector e_j and a plane x_i = 0
    to the line spanned by Vector and the plane given by linear_form = 0

    INPUT:
        Vector      - vector with 3 entries
        linear_form - linear form with linear_form(Vector) = 0

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
    T1 = _ult_line_transformation(base_field, Vector)
    T2 = _ult_plane_transformation(_apply_matrix(T1.inverse(), linear_form))[0]

    return T2 * T1


def _uut_flag_transformation(Vector, linear_form):
    """
    Return unipotent upper triangular matrix transforming a flag given
    by a line spanned by a standard basis vector e_j and a plane x_i = 0
    to the line spanned by Vector and the plane given by linear_form = 0

    INPUT:
        Vector      - vector with 3 entries
        linear_form - linear form with linear_form(Vector) = 0

    OUTPUT:
        T - unipotent upper triangular matrix with e_j*T = Vector and
            _apply_matrix(T, linear_form) = x_i

    MATHEMATICAL INTERPRETATION:
        ...
    """

    Vector = list(Vector)
    if linear_form(Vector) != 0:
        raise ValueError

    base_field = linear_form.base_ring()
    T1 = _uut_line_transformation(base_field, Vector)
    T2 = _uut_plane_transformation(_apply_matrix(T1.inverse(), linear_form))[0]

    return T2 * T1


def _integral_flag_transformation(Vector, linear_form, weight_vector):
    """
    Return unipotent matrix, T, transforming a flag given by a line spanned
    by a standard basis vector e_j and a plane x_i = 0 to the line spanned
    by Vector and the plane given by linear_form = 0 such that for all
    indices i,j the implication
        (weight_vector[j] - weight_vector[i] < 0) => (T[i][j] = 0)
    holds

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
    T2 = _integral_plane_transformation(_apply_matrix(T1.inverse(), linear_form), weight_vector)[0]

    return T2 * T1
