
# ****************************************************************************
#       Copyright (C) 2025 Kletus Stern <sternwork@gmx.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************



from sage.all import *
import LinearValuations as LV



class StabilityFunction:
    def __init__( self, homogeneous_form, base_ring_valuation ):
        """
        INPUT:
            homogeneous_form - homogeneous polynomial in K[x_0,...,x_n]
            base_ring_valuation - discrete valuation on K
        """

        if homogeneous_form.base_ring() != base_ring_valuation.domain():
            raise ValueError(f"Base rings of {homogeneous_form} and {base_ring_valuation} are not equal")

        self.homogeneous_form       = homogeneous_form
        self.base_ring_valuation    = base_ring_valuation

        self.standard_basis     = self.homogeneous_form.parent().gens()
        self.polynomial_ring    = self.homogeneous_form.parent()
        self.base_ring          = self.polynomial_ring.base_ring()
        self.dimension          = len(self.standard_basis) - 1


    def __repr__(self):
        return f"Stability Function of {self.homogeneous_form} over {self.base_ring} with {self.base_ring_valuation}"


    def get_homogeneous_form(self):
        return self.homogeneous_form


    def get_base_ring_valuation(self):
        return self.base_ring_valuation


    def get_dimension(self):
        return self.dimension


    def get_standard_basis(self):
        return self.standard_basis


    def get_polynomial_ring(self):
        return self.polynomial_ring


    def get_base_ring(self):
        return self.base_ring


    def graded_reduction_at(self, point_on_BTB):

        A = point_on_BTB.get_base_change_matrix()
        u = point_on_BTB.weight_vector()
        linear_valuation = LV.LinearValuation(self.polynomial_ring, self.base_ring_valuation, A, u)
        return linear_valuation.graded_reduction_of(self.homogeneous_form)


    def has_semistable_reduction_at(self, point_on_BTB):
        if self.graded_reduction_at(point_on_BTB).is_semistable():
            return True
        return False


    def affine_functions_on_apartment(self, base_change_matrix, affine_patch = None):
        """
        Return the stability function restricted to the apartment given by 'self.standard_basis*base_change_matrix.inverse()'.

        INPUT:
            base_change_matrix - invertible matrix in GL_{self.dimension + 1}(self.base_ring)
            affine_patch - integer between 0 and self.dimension

        OUTPUT:
            [affine_function_1, affine_function_2, ...]

        MATHEMATICS and IMPLEMENTATION:
            We will now explain the mathematics and its implementation in Sage. First, let
                v_K = self.base_ring_valuation,
                A   = base_change_matrix,
                B   = base_change_matrix.inverse(),
                E_0 = (x_0,...,x_n) = self.standard_basis,
                F   = self.homogeneous_form .
            Thus, F is a homogeneous polynomial in K[x_0,...,x_n]. Fruther, we call E_0 the standard
            basis and consider A and B as linear transformations, with respect to E_0, i.e.
                A(x_j) = sum_{i=0}^n a_{ij}*x_i  and  B(x_j) = sum_{i=0}^n b_{ij}*x_i .
            Then,
                E := (y_0,...,y_n) := ( B(x_0),...,B(x_n) )
            is a new basis of K[x_0,...,x_n]. Now if we view E_0 = (x_0,...,x_n) as a vector in Sage,
            we get
                (y_0,...,y_n) = (sum_{i=0}^n b_{i,0}*x_i,...,sum_{i=0}^n b_{i,n}*x_i)
                              = (x_0,...,x_n)*B
            and therefore
                F(x_0,...,x_n) = F( (y_0,...,y_n)*B^{-1} ) = F( (y_0,...,y_n)*A ) .
            Thus, the homogeneous polynomial
                G(y_0,...,y_n) := F( (y_0,...,y_n)*A ) in K[y_0,...,y_n]
            describes F with respect to the basis (y_0,...,y_n) and A describes the base change.
            Thus, for the valuation v_{E,w} we obtain
                v_{E,w}(F) = min( v_K(a_i) + <i,w> : i in I ) with G = sum_{i in I} a_i y^i,
            where i is a multi-index, i.e. I is a subset of NN^{n+1}. Moreover, the we have
                omega(v_{E,w}) = 1/(n+1) * ( w_0 + ... + w_n - v_K( det(E) ) .
            Note that per definition det(E) = det(B). Furthermore,
                v_K( det(B) ) = v_K( det(A^{-1}) ) = v_K( det(A)^{-1} ) = -v_K( det(A) )
            and therefore
                omega(v_{E,w}) = 1/(n+1) * ( w_0 + ... + w_n + v_K( det(A) ) .
            Now let N = n + 1. It follows, that 
                phi_E(w) = v_{E,w}(F) - d*omega(v_{E,w})
                         = min(v_K(a_i) - d/N*v_K(det(A)) + sum_{j=0}^n (i_j - d/N)*w_j : i in I) .
            Finally, we set w_{affine_patch} = 0, if affine_patch != None.
        IMPORTANT REMARK:
            We do not want to introduce new variables and hence we will consider
                G(x_0,...,x_n) = F( (x_0,...,x_n)*T^{-1} )
            and not G(y_0,...,y_n). But for mathematical clarity, we will still refer to (y_0,...,y_n) in the comments below.
        """

        if not base_change_matrix.is_invertible():
            raise ValueError

        # Set up variables
        d = Integer(self.homogeneous_form.degree())
        N = Integer(self.dimension + 1)   # N = n + 1

        # Compute G(x_0,...,x_n) = F( (x_0,...,x_n)*A )
        G = self.homogeneous_form( list( vector( self.standard_basis )*base_change_matrix ) )

        # Now create variables for affine functions
        w = list( PolynomialRing( QQ, N, 'w' ).gens() ) # w = [w_0,...,w_n] since N = n + 1

        # Now set w_j = 0 with j = affine_patch
        if not affine_patch == None:
            if affine_patch < 0 or N - 1 < affine_patch:
                raise ValueError
            else:
                w[affine_patch] = 0

        # Compute d/N*v_K( det(A) )
        const_A = d/N*self.base_ring_valuation(base_change_matrix.det())

        # Compute v_{E,w}(F) - d*omega( v_{E,w} )
        affine_functions = []
        for multi_index, G_coefficient in G.dict().items():
            affine_function = self.base_ring_valuation( G_coefficient ) - const_A
            for j in range(N):
                affine_function = affine_function + ( multi_index[j] - d/N )*w[j]
            affine_functions.append( affine_function )

        return affine_functions


    def active_functions_at(self, base_change_matrix, weight_vector):
        return RestrictedStabilityFunction(self, base_change_matrix).active_functions(weight_vector)


    def _maximum_on_apartment(self, base_change_matrix, affine_patch):
        """
        Return the maximum of the stability function on the apartment given by 'self.standard_basis*base_change_matrix.inverse()'.

        INPUT:
            base_change_matrix - invertible matrix in GL_{self.dimension + 1}(self.base_ring)
            affine_patch - integer between 0 and self.dimension

        OUTPUT:
            rational number, which equals the maximum of self on the apartment given by
            the basis 'self.standard_basis*base_change_matrix.inverse()'
        """

        affine_functions = self.affine_functions_on_apartment(base_change_matrix, affine_patch)

        # Note that by definition of the method 'affine_functions_on_apartment()' the elements
        # of 'affine_functions', i.e. affine_function_i's are linear polynomials in QQ[w_0,...,w_n]
        # and hence of the form q_0*w_0 + ... + q_n*w_n + q_constant with q_{ self.affine_patch } = 0
        affine_function_variables = list( affine_functions[0].parent().gens() ) # [w_0,...,w_n]
        N = self.dimension + 1

        MILP = MixedIntegerLinearProgram(solver='PPL')  # we need solver='PPL' for an exact rational solution
        v = MILP.new_variable()
        t = v['maximum']    
        MILP.set_objective(t)   # t will be maximized the under constraints below
 
        # We have affine_function = q_0*w_0 + ... q_n*w_n + q_constant. The next for-loop will replace w_i's by MILP variables u_i's.
        for affine_function in affine_functions:
            # affine_function = q_0*w_0 + ... + q_n*w_n + q_constant and hence affine_function.constant_coefficient() = q_constant
            affine_function_in_MILP_variables = affine_function.constant_coefficient()  # the constant term of affine_function
            for i in range(N):
                # w[i] = w_i is already a monomial in affine_function = q_0*w_0 + ... + q_n*w_n + q_constant
                # and the corresponding monomial coefficient is q_i.
                affine_function_in_MILP_variables = affine_function_in_MILP_variables + affine_function.monomial_coefficient( affine_function_variables[i] )*v["u"+str(i)]
            MILP.add_constraint( t <= affine_function_in_MILP_variables )
        MILP.solve()

        return MILP.get_values(v)


    def maximum_on_apartment(self, base_change_matrix, affine_patch = None):
        """
        Return the maximum on the apartment given by base_change_matrix and
        the point where the maximum is reached
        """

        if affine_patch != None:
            return self._maximum_on_apartment(base_change_matrix, affine_patch)

        solution_dict = self.maximum_on_apartment(base_change_matrix, 0)
        maximum = solution_dict['maximum']
        weight_vector = [solution_dict['u'+str(i)] for i in range(self.dimension + 1)]
        point_on_BTB = BTB_Point(self.base_ring_valuation, base_change_matrix, weight_vector)
        return [maximum, point_on_BTB]


    def ascent_directions_at(self, point_on_BTB, matrix_form = 'uut'):
        """
        Return a list of matrices which describe base changes, fixing 'point_on_BTB', to apartment,
        where the stability function can be maximized further.

        INPUT:
            point_on_BTB - object in the class 'BTB_Point' such that 'point_on_BTB.get_base_ring()'
                           equals 'self.base_ring'
            matrix_form  - one of the strings 'ult', 'uut', 'integral'

        OUTPUT:
            list of invertible matrix in GL_{self.dimension + 1}(self.base_ring)

        MATHEMATICAL INTERPRETATION:
            First, let
                A   = point_on_BTB.get_base_change_matrix(),
                u   = point_on_BTB.weight_vector(),
                B   = A.inverse(),
                K   = self.base_ring
                n   = self.dimension
                E_0 = self.standard_basis .
            Then the matrix A lies in GL_n(K) and point_on_BTB is represented by the valuation v_{E,u}.
            Now we view E_0 = (x_0,...,x_n) as a vector in Sage and define the basis
                E_1 = (x_0,...,x_n)*B .
            Furthermore, let phi be the stability function, 'self'. We want to find a matrix T in GL_n(K)
            such that for the basis E_2 = E_1*T^(-1) we get
                phi(point_on_BTB) < max phi_{E_2},
            where phi_{E_2} is the stability function phi restricted to the apartment given by E_2.
            Note that if n = 2, we can compute E_2 by computing a graded instability (E_2, w) of the graded
            reduction of 'self.homogeneous_form' with respect to a valuation representing 'point_on_BTB'.
        """

        if self.dimension != 2:
            raise NotImplementedError
        if self.base_ring != point_on_BTB.get_base_ring():
            raise ValueError

        # (1) Compute the graded reduction, say f, of self.homogeneous_form at point_on_BTB.
        # (2) Compute the list, say L, of all graded instabilities of f.
        # (3) Find all rational graded instability in L and return their lifted matrices. Otherwise return None.

        graded_reduction = self.graded_reduction_at(point_on_BTB)

        ascent_directions = []
        for rational_graded_instability in graded_reduction.rational_graded_instabilities(matrix_form):
            ascent_directions.append(rational_graded_instability.lift_matrix())
        return ascent_directions


    def ascent_direction_at(self, point_on_BTB, matrix_form = 'integral'):
        """
        Return a matrix which describes a base change, fixing 'point_on_BTB', to an apartment,
        where the stability function can be maximized further.

        INPUT:
            point_on_BTB - object in the class 'BTB_Point' such that 'point_on_BTB.get_base_ring()'
                           equals 'self.base_ring'
            matrix_form  - one of the strings 'ult', 'uut', 'integral'

        OUTPUT:
            invertible matrix in GL_{self.dimension + 1}(self.base_ring) 
        """

        ascent_directions = self.ascent_directions_at(point_on_BTB, matrix_form)

        if len(ascent_directions) == 0:
            return None
        else:
            return ascent_directions[0]


    def maximize(self, matrix_form = 'integral'):
        """
        Return the maximum and the point where the stability function takes it

        INPUT:
            matrix_form - one of the strings 'ult', 'uut', 'integral'
        """

        global_trafo_matrix = identity_matrix(self.base_ring, self.dimension + 1)

        # find maximum on standard apartment
        solution_dict = self.maximum_on_apartment(global_trafo_matrix, 0)
        maximum = solution_dict['maximum']
        weight_vector = [solution_dict['u'+str(i)] for i in range(self.dimension + 1)]

        point_on_BTB = BTB_Point(self.base_ring_valuation, global_trafo_matrix, weight_vector)

        if point_on_BTB.minimal_simplex_dimension() == self.dimension:
            return [maximum, point_on_BTB]

        affine_patch = point_on_BTB.affine_patch()
        local_trafo_matrix = self.ascent_direction_at(point_on_BTB, matrix_form)

        if local_trafo_matrix == None:
            return [maximum, point_on_BTB]

        while True:
            global_trafo_matrix = local_trafo_matrix * global_trafo_matrix

            # find maximum on the new apartment
            solution_dict = self.maximum_on_apartment(global_trafo_matrix, 0)
            maximum = solution_dict['maximum']
            weight_vector = [solution_dict['u'+str(i)] for i in range(self.dimension + 1)]

            point_on_BTB = BTB_Point(self.base_ring_valuation, global_trafo_matrix, weight_vector)

            if point_on_BTB.minimal_simplex_dimension() == self.dimension:
                return [maximum, point_on_BTB]

            affine_patch = point_on_BTB.affine_patch()
            local_trafo_matrix = self.ascent_direction_at(point_on_BTB, matrix_form)

            if local_trafo_matrix == None:
                return [maximum, point_on_BTB]

#            if local_trafo_matrix.is_diagonal():
#                return [maximum, BTB_Point(self.base_ring_valuation, global_trafo_matrix, weight_vector)]



    def evaluate_at(self, point_on_BTB):
        """
        Return the value of self at 'point_on_BTB'
        """

        T = point_on_BTB.get_base_change_matrix()
        w = point_on_BTB.weight_vector()

        return min(affine_function(w) for affine_function in self.affine_functions_on_apartment(T))



    def show( self ):
        print( self.affine_functions()[0].parent().gens(), "|---> min" + str( tuple( self.affine_functions() ) ) )



class RestrictedStabilityFunction:
    def __init__(self, stability_function, base_change_matrix):
        r"""
        Return ...
        """

        self.homogeneous_form = stability_function.get_homogeneous_form()
        self.base_ring_valuation = stability_function.get_base_ring_valuation()
        self.base_change_matrix = base_change_matrix


    def __repr__(self, ):
        n = len(self.homogeneous_form.parent().gens())
        return f"Stability Function of {self.homogeneous_form} over {self.base_ring} with {self.base_ring_valuation} restricted to RR^{n}"


    def maximum(self):
        raise NotImplementedError


    def active_functions(self, w, flag = True):
        r"""
        Return the set of active functions a w
        """

        d = self.homogeneous_form.degree()
        N = len(self.homogeneous_form.parent().gens())
        # Compute d/N*v_K( det(A) )
        const_A = d/N*self.base_ring_valuation(self.base_change_matrix.det())
        affine_functions_values = dict()
        F = _apply_matrix(self.base_change_matrix, self.homogeneous_form)
        for multi_index, coefficient in F.dict().items():
            value_at_w = self.base_ring_valuation(coefficient) - const_A
            for j in range(N):
                value_at_w = value_at_w + multi_index[j] * w[j]
            affine_functions_values[multi_index] = value_at_w
        min_value = min(affine_functions_values.values())
        return [key for key, value in affine_functions_values.items() if value == min_value]


    def embedding_matrix(self):
        r"""
        Return the matrix T such that...
        """
        return self.base_change_matrix



class BTB_Point:
    def __init__(self, base_ring_valuation, base_change_matrix, weight_vector):

        self.base_ring_valuation = base_ring_valuation
        self.base_change_matrix = base_change_matrix

        # convert all entries to rationals
        weight_vector_qq = [QQ(w) for w in weight_vector]

        normalized_weight_vector = []
        w_min = min(weight_vector_qq)
        for w in weight_vector_qq:
            normalized_weight_vector.append(w - w_min)

        self._weight_vector = normalized_weight_vector


    def __repr__(self):
        return f"Point on the Bruhat-Tits Building of SL({len(self._weight_vector)}) over {self.base_ring_valuation.domain()} with {self.base_ring_valuation}"


    def get_base_change_matrix(self):
        return self.base_change_matrix


    def weight_vector(self):
        return self._weight_vector


    def affine_patch(self):
        return self._weight_vector.index(0)


    def get_base_ring(self):
        return self.base_change_matrix.base_ring()


    def minimal_simplex_dimension(self):
        """
        Return the dimension of the minimal simplex containing self

        MATHEMATICAL INTERPRETATION (ToDo. But at this point we just give an example.): 
            If we start with the point (0, 1/2, 3/2, 1/3), we first move it to
            (0, 1/2, 1/2, 1/3). Then we see that {(0,0,0,0), (0, 1, 1, 0), (0, 1, 1, 1)}
            is the minimal simplex containing (0, 1/2, 1/2, 1/3). This is because we have
            to consider the lines <(0, 1, 0, 0)>, <(0, 0, 1, 0)>, <(0, 0, 0, 1)>,
            <(0, 1, 1, 0)>, <(0, 1, 0, 1)>, <(0, 0, 1, 1)>, <(0, 1, 1, 1)> and to find the
            minimal direct sum, containing (0, 1/2, 1/2, 1/3). In our case this is the sum
                <(0, 1, 1, 0)> oplus <(0, 0, 0, 1> = <(0,1,1,0), (0, 0, 0, 1)>.
            Since the simplices need to have a cyclic basis, we must choose
            (0, 1, 1, 0), (0, 1, 1, 1) as a basis of the direct sum above.
            In general, we first have to normalize the valuation to have the value group ZZ,
            i.e. to divide the vector 'self.weight_vector()' by v_K(pi_K). Let w_normalized be
            this normalized vector. Now we have to move w_normalized inside the cube [0,1]^N.
            Let w be this translated vector. Now remove all zeros from w. Let w_without_zeros be
            this modified vector. Now compute the cardinality of the set set(w_without_zeros).
            This cardinality is exactly the dimension of the minimal simplex containing self.
            Note that this is the same as the cardinality of the set set(w) minus 1, i.e.
            len(set(w)) - 1.
        """

        # normalize the value grout to be ZZ
        value_groug_generator = self.base_ring_valuation.value_group().gen()
        norm_weight_vector = []
        for i in self._weight_vector:
            norm_weight_vector.append(i/value_groug_generator)

        # translate inside the unit cube
        trans_norm_weight_vector = []
        for i in norm_weight_vector:
            trans_norm_weight_vector.append(i - floor(i))

        return len(set(trans_norm_weight_vector)) - 1





class BTB_apartment_point:
    def __init__(self, base_ring_valuation, base_change_matrix, affine_patch, weight_vector):
        """
        ToDo
        """



# ===================== helper functions =====================

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
