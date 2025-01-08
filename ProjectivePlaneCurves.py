from sage.all import *



class ProjectivePlaneCurve:

    def __init__(self, polynomial):
        """
        INPUT:
            ...
        """

        if not polynomial.is_homogeneous():
            return TypeError

        if not len(polynomial.parent().gens()) == 3:
            return ValueError

        self.polynomial = polynomial
        self.degree = self.polynomial.degree()
        self.polynomial_ring = self.polynomial.parent()
        self.base_ring = self.polynomial.base_ring()
        self.projective_plane = ProjectiveSpace(self.polynomial_ring)
        self.plane_curve = self.projective_plane.curve(self.polynomial)
        self.standard_basis = self.polynomial_ring.gens()




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



    def instable_weight_vector_wrt(self, base_change_matrix):
        """
        Return an instable weight vector with respect to the basis given by
        base_change_matrix or 'None' if no such vector exist.

        INPUT:
            base_change_matrix - 3x3 invertible matrix over self.base_ring

        MATHEMATICAL INTERPRETATION:
            First, let

                K = self.base_ring,
                T = base_change_matrix,
                F = self.get_polynomial() .

            Furthermore, let

                (x0, x1, x2) = self.get_standard_basis()

            and

                G = F((x0, x1, x2) * T) .

            For a multi index set I subset NN^3 we can write

                G = sum_{i in I} a_i x0^i0 * x1^i1 * x2^i2

            with a_i != 0 for all i in I.

            Now, we want to find an ordered and balanced weight vector

                (w0, w1, w2) in QQ^3

            for G, i.e. w0 => w1 => w2, w0 + w1 + w2 = 0 and

                i0*w0 + i1*w1 + i2*w2 > 0

            for all i in I.

        REMARK. An ordered and balanced weight vector for G exists if

        and only if there exists a balanced weight vector (w0, w1, w2)

        for G such that

            -1 <= w0, w1, w2 <= 1.

        This follows from the fact, that any nonzero multiple of a

        balanced weight vector is again a balanced weight vector.

        Thus, we only have to maximize the function

            min(i0*w0 + i1*w1 + i2*w2 : i in I)

        under the constraints -1 <= w0, w1, w2 <= 1 and to check

        whether the maximum is > 0 or not.
        """

        if not base_change_matrix.is_invertible():
            raise ValueError

        G = self.polynomial(list(vector(self.standard_basis) * base_change_matrix))

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

        if values['maximum'] > 0:
            return (values['w0'], values['w1'], values['w2'])

        return None



    def _move_point_to_affine_origin( self, projective_point ):
        """
        Return a list consisting of an integer i with 0 <= i <=2 and the upper
        triangular unipotent matrix T with maximal number of zeros such that
        e_i * T = projective_point, where e_i is the i-th standard basis vector.

        INPUT:
            projective_point - [a, b, c]

        MATHEMATICAL INTERPRETATION:
            If a is not zero, we can assume a = 1. Then i = 0 and the matrix

                T = [[1, b, c],
                     [0, 1, 0],
                     [0, 0, 1]]

            satisfies (1, 0, 0) * T = (1, b, c). If a = 0, but b is not zero, we can

            assume a = 0 and b = 1. Then i = 1 and the matrix

                T = [[1, 0, 0],
                     [0, 1, c],
                     [0, 0, 1]]

            satisfies (0, 1, 0) * T = (a, b, c). Finally, if a = b = 0, we can assume

            c = 1. Then i = 2 and T is the identity matrix.
        """

        T = [[1,0,0],[0,1,0],[0,0,1]]

        # Find the minimal index, say i_min, with projective_point[i_min] != 0 
        L_proj_point = list(projective_point)
        i_min = L_proj_point.index(next(x for x in L_proj_point if x != 0))

        # Normalize projective_point
        for i in range(3):
            L_proj_point[i] = L_proj_point[i] / L_proj_point[i_min]

        T[i_min] = L_proj_point

        return [i_min, matrix(self.base_ring, T)]



    def _move_affine_line_to_coordinate_axis(self, line, affine_patch):
        """
        Return an upper triangular unipotent matrix with maximal number of zeros which describes
        a projective transformation moving the affine line, 'line', to one of the affine coordinate
        axes in the affine patch 'affine_patch'.

        INPUT:
            line         - A*x + B*y with A, B elements of self.base_ring
            affine_patch - 0 or 1 or 2

        OUTPUT:
            matrix over self.base_ring

        MATHEMATICAL INTERPRETATION:
            First, let

                k = self.base_ring,
                L = line .

            Then we have

                L = A*x + B*y in k[x,y] .

            If A = 0 or B = 0, the affine line L is already equal to one of the

            two coordinate axes. Thus, we can assume that A != 0 and B != 0. In

            that case the matrix

                M = [[1, -A/B],
                     [0,    1]]

            describes an affine transformation moving L to A*x. In fact, we have

                (x,y) * M = (x, -A / B * x + y)

            and therefore

                L((x,y) * M) = A*x + B*(-A / B * x + y) = B*y .

            Furthermore, depending on whether

                affine_patch = 0 or affine_patch = 1 or affine_patch = 2,

            the corresponding projective transformation if given by

                [1, 0,    0]
                [0, 1, -A/B]
                [0, 0,    1]

            or

                [1, 0, -A/B]
                [0, 1,    0]
                [0, 0,    1]

            or

                [1, -A/B, 0]
                [0,    1, 0]
                [0,    0, 1]

            respectively.

        REMARK: This function will be used in combination with the

        function _move_point_to_affine_origin(self, (a:b:c)) which

        returns one of the following matrices

            [1, b, c]        [1, 0, 0]        [1, 0, 0]
            [0, 1, 0]   or   [0, 1, c]   or   [0, 1, 0]
            [0, 0, 1]        [0, 0, 1]        [0, 0, 1]

        depending on whether

            a != 0 or b != 0 or c != 0 .

        Thus, the composition of the matrices given by the functions

            '_move_point_to_affine_origin'

        and

            '_move_projective_line_to_coordinate_axis'

        will be again an upper triangular unipotent matrix.
        """

        T = [[1,0,0],[0,1,0],[0,0,1]]
        line_coefficients = line.coefficients()

        if len(line_coefficients) == 1:
            return matrix(self.base_ring, T)

        if affine_patch == 0:
            T[1][2] = -line_coefficients[0] / line_coefficients[1]

        elif affine_patch == 1:
            T[0][2] = -line_coefficients[0] / line_coefficients[1]

        else:
            T[0][1] = -line_coefficients[0] / line_coefficients[1]

        return matrix(self.base_ring, T)



    def _move_projective_line_to_coordinate_axis(self, projective_line):
        """
        Return the list of tupels, t, where t[1] run over all upper triangular
        unipotent matrices with maximal number of zeros which transform
        'projective_line' to a coordinate axis and t[0] is the index of the
        corresponding coordinate axis.

        MATHEMATICAL INTERPRETATION:
            First, let

                L = projective_line .

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

                [[1, -A/B,  0],
                 [0,   1, -B/C],
                 [0,   0,   1]]

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

        if projective_line.degree() != 1:
            raise ValueError

        if projective_line.parent() != self.polynomial_ring:
            raise ValueError

        x0, x1, x2 = self.standard_basis
        A = projective_line.monomial_coefficient(x0)
        B = projective_line.monomial_coefficient(x1)
        C = projective_line.monomial_coefficient(x2)
        L = []

        if C != 0:
            coordinate_axis_index = 2
            T = [[1, 0, -A/C], [0, 1, -B/C], [0, 0, 1]]
            T = matrix(self.base_ring, T)
            L.append((coordinate_axis_index, T))
            if B != 0 and A != 0:
                T = [[1, -A/B, 0], [0, 1, -B/C], [0, 0, 1]]
                T = matrix(self.base_ring, T)
                L.append((coordinate_axis_index, T))

        elif B != 0:
            coordinate_axis_index = 1
            T = [[1, -A/B, 0], [0, 1, 0], [0, 0, 1]]
            T = matrix(self.base_ring, T)
            L.append((coordinate_axis_index, T))

        else:
            coordinate_axis_index = 0
            T = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            T = matrix(self.base_ring, T)
            L.append((coordinate_axis_index, T))

        return L



    def instabilities(self):
        """
        Return a list of all instabilities, given by unipotent upper triangular
        matrices with maximal number of zeros, of self.

        MATHEMATICAL INTERPRETATION:
            We will explain how to compute all such instabilities and why there are
            only finitely many of them.

        """

        list_of_instabilities = []


        # CASES (b) and (d)
        for P_initial, m in self.points_with_high_multiplicity():
            if m > 2 * self.degree / 3:  # CASE (b)
                affine_patch, T1 = self._move_point_to_affine_origin(P_initial)
                weight_vector = [1, 1, 1]
                weight_vector[affine_patch] = -2
                weight_vector = tuple(weight_vector)
                list_of_instabilities.append(PPC_Instability(T1, weight_vector, (P_initial, m)))

            elif m > self.degree / 2:  # CASE (d) 1/2
                affine_patch, T1 = self._move_point_to_affine_origin(P_initial)

                for line, line_multiplicity in PPC_TangentCone(self, P_initial).get_lines():
                    if line_multiplicity > m / 2:
                        T2 = self._move_affine_line_to_coordinate_axis(line, affine_patch)

                        weight_vector = self.instable_weight_vector_wrt(T2*T1)
                        if not weight_vector == None:
                            list_of_instabilities.append(PPC_Instability(T2 * T1, weight_vector, 'Case (d)'))

        # CASES (a) and (c)
        for line, line_multiplicity, G in self.lines_with_high_multiplicity():
            if line_multiplicity > self.degree / Integer(3):  # CASE (a)
                for coordinate_axis_index, T in self._move_projective_line_to_coordinate_axis(line):
                    weight_vector = [-1, -1, -1]
                    weight_vector[coordinate_axis_index] = 2
                    weight_vector = tuple(weight_vector)
                    list_of_instabilities.append(PPC_Instability(T, weight_vector, (line, line_multiplicity)))

            else: # CASE (c) 1/2
                line_curve = self.projective_plane.curve(line)
                G_curve = self.projective_plane.curve(G)
                for P_initial in line_curve.intersection_points(G_curve):
                    if line_curve.intersection_multiplicity(G_curve, P_initial) > (self.degree - line_multiplicity) / Integer(2):
                        affine_patch, T1 = self._move_point_to_affine_origin(P_initial)
                        R = PolynomialRing(self.base_ring, ['x', 'y'])
                        x, y = R.gens()
                        dehomogenization = [x,y]
                        dehomogenization.insert(affine_patch,R(1))
                        dehomogenization = vector(dehomogenization)
                        affine_line = line( list(dehomogenization * T1) )

                        T2 = self._move_affine_line_to_coordinate_axis(affine_line, affine_patch)

                        weight_vector = self.instable_weight_vector_wrt(T2 * T1)
                        if weight_vector != None:
                            geometric_type = "Case (c) with " + str(line) + " and " + str(P_initial)
                            list_of_instabilities.append(PPC_Instability(T2 * T1, weight_vector, geometric_type))

        return list_of_instabilities





class PPC_TangentCone:

    def __init__(self, projective_plane_curve, point_on_PPC):
        """
        At the moment only rational points are allowed
        """

        if projective_plane_curve.get_polynomial()(list(point_on_PPC)) != 0:
            return ValueError

        self.projective_plane_curve = projective_plane_curve
        self.base_ring = projective_plane_curve.get_base_ring()
        self.polynomial_ring = PolynomialRing(self.base_ring, ['x', 'y'])
        self.gen1, self.gen2 = self.polynomial_ring.gens()

        # Convert to list
        point_on_PPC = list(point_on_PPC)

        # Coerce coordinates to self.base_ring
        for i in range(3):
            point_on_PPC[i] = self.base_ring(point_on_PPC[i])

        # Find the minimal index, say i_min, with point_on_PPC[i_min] != 0 
        i_min = point_on_PPC.index(next(x for x in point_on_PPC if x != 0))

        # Normalize point_on_PPC
        for i in range(3):
            point_on_PPC[i] = point_on_PPC[i] / point_on_PPC[i_min]

        self.point_on_PPC = point_on_PPC
        self.affine_patch = i_min


    def get_polynomial(self):

        PPC_equation = self.projective_plane_curve.get_polynomial()
        dehomogenization = [self.gen1, self.gen2]
        dehomogenization.insert(self.affine_patch, self.polynomial_ring(0))
        dehomogenization_translated = []

        for i in range(3):
            dehomogenization_translated.append(dehomogenization[i] + self.point_on_PPC[i])

        f = PPC_equation(dehomogenization_translated)

        f_homo_comp_dict = f.homogeneous_components()
        minimal_degree = min(f_homo_comp_dict.keys())
        tangent_cone_polynomial = f_homo_comp_dict[minimal_degree]

        return tangent_cone_polynomial



    def get_lines(self):

        L = []
        tangent_cone_polynomial = self.get_polynomial()
        for factor, factor_multiplicity in list(tangent_cone_polynomial.factor()):
            if factor.degree() == 1:
                L.append((factor, factor_multiplicity))

        return L
                
                
        



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
