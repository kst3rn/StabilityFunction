
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
import re


class LinearValuation:
    
    def __init__(self, polynomial_ring, base_ring_valuation, base_change_matrix, weight_vector):
        """
        Return the valuation on polynomial_ring diagonalized by the basis
        vector( polynomial_ring.gens() )*base_change_matrix.inverse()
        with respect to weight vector weight_vector.

        INPUT:
            polynomial_ring - polynomial ring
            base_ring_valuation - valuation on polynomial_ring.base_ring()
            base_change_matrix - invertible matrix over polynomial_ring.base_ring()
            weight_vector - vector over polynomial_ring.base_ring()

        OUTPUT:
            v_{E,u} - linear valuation on polynomial_ring with, where:
                        w = weight_vector,
                        E = vector( polynomial_ring.gens() )*base_change_matrix.inverse()

        INTERPRETATION:
            First, let
                K   = polynomial_ring.base_ring(),
                E_0 = (x_0,...,x_n) = polynomial_ring.gens(),
                v_K = base_ring_valuation,
                A   = base_change_matrix,
                T   = base_change_matrix.inverse(),
                w   = weight_vector .
            Thus, we can write
                polynomial_ring = K[x_0,...,x_n] .
            We call E_0 the standard basis and consider A and T as linear transformations,
            with respect to E_0, i.e.
                A(x_j) = sum_{i=0}^n a_{ij}*x_i  and  T(x_j) = sum_{i=0}^n t_{ij}*x_i .
            Then E = (y_0,...,y_n) = ( T(x_0),...,T(x_n) ) is a new basis of K[x_0,...,x_n].
            and we set
                LinearValuation( K[x_0,...,x_n], v_K, A, w ) = v_{E,u} .
            Now let F be a polynomial in K[x_0,...,x_n]. To compute v_{E,u}(F) we have to
            write F with respect to the basis E. This will be done as follows. If we view
            E_0 = (x_0,...,x_n) as a vector in Sage, we get
                (x_0,...,x_n)*T = (sum_{i=0}^n t_{i,0}*x_i,...,sum_{i=0}^n t_{i,n}*x_i)
                                = (y_0,...,y_n)
            and therefore
                F(x_0,...,x_n) = F( (y_0,...,y_n)*T^{-1} ) = F( (y_0,...,y_n)*A ) .
            Thus, the polynomial
                G(y_0,...,y_n) = F( (y_0,...,y_n)*A ) in K[y_0,...,y_n]
            describes F with respect to the basis (y_0,...,y_n) and A describes the base change.
            Thus, with the notation
                <i,u> = i_0*u_0 + ... + i_n*u_n,
            we obtain
                v_{E,u}(F) = min( v_K(a_i) + <i,u> : i in I ) with G = sum_{i in I} a_i y^i,
            where i is a multi-index, i.e. I is a subset of NN^{n+1}.
            Remark. Note that the valuations
                LinearValuation( K[x_0,...,x_n], v_K, A, w )
            and
                LinearValuation( K[y_0,...,y_n], v_K, identity_matrix, w )
            are equal. We will use this fact for example to compute graded reduction.
        IMPORTANT REMARK:
            We do not want to introduce new variables and hence we will consider
                G(x_0,...,x_n) = F( (x_0,...,x_n)*T^{-1} )
            and not G(y_0,...,y_n). But for mathematical clarity, we will still
            refer to (y_0,...,y_n) in the comments below.
        """

        # convert all entries to rationals
        weight_vector_qq = [QQ(w) for w in weight_vector]
    
        self.polynomial_ring = polynomial_ring
        self.base_ring_valuation = base_ring_valuation
        self.base_change_matrix = base_change_matrix
        self.weight_vector = weight_vector_qq


    def __repr__(self):
        return "Linear valuation"


    def __call__(self, polynomial):
        return self.evaluate_at(polynomial)


    def get_polynomial_ring(self):
        return self.polynomial_ring


    def get_base_field(self):
        return self.polynomial_ring.base_ring()


    def get_base_ring_valuation(self):
        return self.base_ring_valuation


    def get_dimension(self):
        return len(self.weight_vector)


    def standard_basis(self):
        return self.polynomial_ring.gens()


    def get_weight_vector(self):
        return self.weight_vector


    def get_base_change_matrix(self):
        return self.base_change_matrix


    def evaluate_at( self, polynomial ):
        """
        Return the value of self at polynomial, i.e. the valuation of polynomial

        INPUT:
            polynomial - element of self.polynomial_ring

        OUTPUT:
            rational number

        MATHEMATICAL INTERPRETATION:
            First, let
                K   = self.polynomial_ring.base_ring(),
                E_0 = (x_0,...,x_n) = self.polynomial_ring.gens(),
                v_K = self.base_ring_valuation,
                A   = self.base_change_matrix,
                T   = self.base_change_matrix.inverse(),
                u   = self.weight_vector .
            Thus, we can write
                self.polynomial_ring = K[x_0,...,x_n] .
            We call E_0 the standard basis and consider A and T as linear transformations,
            with respect to E_0, i.e.
                A(x_j) = sum_{i=0}^n a_{ij}*x_i  and  T(x_j) = sum_{i=0}^n t_{ij}*x_i .
            Then E = (y_0,...,y_n) = ( T(x_0),...,T(x_n) ) is a new basis of K[x_0,...,x_n].
            and
                LinearValuation( K[x_0,...,x_n], v_K, A, u ) = v_{E,u} .
            Let F be a polynomial in K[x_0,...,x_n]. To evaluate v_{E,u} at F we have to
            write F with respect to the basis E. This will be done as follows. If we view
            E_0 = (x_0,...,x_n) as a vector in Sage, we get
                (x_0,...,x_n)*T = (sum_{i=0}^n T_{i,0}*x_i,...,sum_{i=0}^n T_{i,n}*x_i)
                                = (y_0,...,y_n)
            and therefore
                F(x_0,...,x_n) = F( (y_0,...,y_n)*T^{-1} ) = F( (y_0,...,y_n)*A ) .
            Thus, the polynomial
                G(y_0,...,y_n) = F( (y_0,...,y_n)*A ) in K[y_0,...,y_n]
            describes F with respect to the basis (y_0,...,y_n) and A describes the base change.
            For a multi index set I subset NN^(n+1) we can write
                G = sum_{i in I} a_i y_0^i_0 * ... * y_n^i_n
            Thus, with the notation
                <i,u> = i_0*u_0 + ... + i_n*u_n,
            we obtain
                v_{E,u}(F) = min( v_K(a_i) + <i,u> : i in I ).
            Remark. Note that the valuations
                LinearValuation( K[x_0,...,x_n], v_K, A, u )
            and
                LinearValuation( K[y_0,...,y_n], v_K, identity_matrix, u )
            are equal. We will use this fact for example to compute graded reduction.

        IMPORTANT REMARK:
            We do not want to introduce new variables and hence we will consider
                G(x_0,...,x_n) = F( (x_0,...,x_n)*T^{-1} )
            and not G(y_0,...,y_n). But for mathematical clarity, we will still
            refer to (y_0,...,y_n) in the comments below.
        """
        if not polynomial.parent() == self.polynomial_ring:
            raise TypeError(f"Polynomial must be in {self.polynomial_ring}, but got parent {polynomial.parent()}")

        if polynomial == 0:
            return +Infinity
        
        N = self.get_dimension()
        G = polynomial( list( vector(self.standard_basis()) * self.base_change_matrix ) )
        values = []
        for multi_index, G_coefficient in G.dict().items():
            value = self.base_ring_valuation( G_coefficient )
            for j in range(N):
                value = value + multi_index[j]*self.weight_vector[j]
            values.append( value )
        return min(values)


    def initial_form(self, polynomial):
        r"""
        Return the initial form of polynomial with respect to self
        INPUT:
            polynomial - element of self.polynomial_ring
        OUTPUT:
            string
        MATHEMATICAL INTERPRETATION:
            First, let
                K   = self.polynomial_ring.base_ring(),
                E_0 = (x_0,...,x_n) = self.polynomial_ring.gens(),
                v_K = self.base_ring_valuation,
                A   = self.base_change_matrix,
                T   = self.base_change_matrix.inverse(),
                u   = self.weight_vector .
            Thus, we can write
                self.polynomial_ring = K[x_0,...,x_n] .
            We call E_0 the standard basis and consider A and T as linear transformations,
            with respect to E_0, i.e.
                A(x_j) = sum_{i=0}^n a_{ij}*x_i  and  T(x_j) = sum_{i=0}^n t_{ij}*x_i .
            Then E = (y_0,...,y_n) = ( T(x_0),...,T(x_n) ) is a new basis of K[x_0,...,x_n].
            and
                LinearValuation( K[x_0,...,x_n], v_K, A, u ) = v_{E,u} .
            Let F be a polynomial in K[x_0,...,x_n]. To evaluate v_{E,u} at F we have to
            write F with respect to the basis E. This will be done as follows. If we view
            E_0 = (x_0,...,x_n) as a vector in Sage, we get
                (x_0,...,x_n)*T = (sum_{i=0}^n T_{i,0}*x_i,...,sum_{i=0}^n T_{i,n}*x_i)
                                = (y_0,...,y_n)
            and therefore
                F(x_0,...,x_n) = F( (y_0,...,y_n)*T^{-1} ) = F( (y_0,...,y_n)*A ) .
            Thus, the polynomial
                G(y_0,...,y_n) = F( (y_0,...,y_n)*A ) in K[y_0,...,y_n]
            describes F with respect to the basis (y_0,...,y_n) and A describes the base change.
            For a multi index set I subset NN^(n+1) we can write
                G = sum_{i in I} a_i y_0^i_0 * ... * y_n^i_n
            Thus, with the notation
                <i,u> = i_0*u_0 + ... + i_n*u_n,
            we obtain
                v_{E,u}(F) = min( v_K(a_i) + <i,u> : i in I ).
            Remark. Note that the valuations
                LinearValuation( K[x_0,...,x_n], v_K, A, u )
            and
                LinearValuation( K[y_0,...,y_n], v_K, identity_matrix, u )
            are equal. The initial form of F with respect to v_{E, u} is the polynomial
                sum_{i in I_0} a_i y_0^i_0 * ... * y_n^i_n
            for
                I_0 = {i \in I : v_K(a_i) + <i,u> = v_{E, u}(F)}.
        """
        polynomial_value = self(polynomial)
        G = _apply_matrix(self.base_change_matrix, polynomial)
        if G.is_zero():
            return "0"

        parent_ring = G.parent()
        gens = parent_ring.gens()
        G_initial = parent_ring(0)
        N = self.get_dimension()
        for multi_index, coefficient in G.dict().items():
            value = self.base_ring_valuation(coefficient)
            value += sum(multi_index[j] * self.weight_vector[j] for j in range(N))
            if value == polynomial_value:
                G_initial += coefficient * prod(gens[j]**multi_index[j] for j in range(N))

        s = str(G_initial)

        var_names = [str(g) for g in gens]
        sorted_var_names = sorted(var_names, key=len, reverse=True)
        for v_name in sorted_var_names:
            pattern = r'\b' + re.escape(v_name) + r'\b'
            replacement = f"T({v_name})"
            s = re.sub(pattern, replacement, s)
        return s


    def graded_reduction_of( self, polynomial ):
        """
        Return the graded reduction of polynomial with respect to self.

        INPUT:
            polynomial  - element of self.polynomial_ring

        OUTPUT:
            object in the class GradedReduction

        MATHEMATICS AND IMPLEMENTATION:
            First, let
                F   = polynomial,
                K   = self.polynomial_ring.base_ring(),
                E_0 = (x_0,...,x_n) = self.polynomial_ring.gens(),
                v_K = self.base_ring_valuation,
                k   = v_K.residue_field()
                p_K = v_K.uniformizer()
                A   = self.base_change_matrix,
                T   = self.base_change_matrix.inverse(),
                u   = self.weight_vector .
            Thus, we can write
                self.polynomial_ring = K[x_0,...,x_n] .
            We call E_0 the standard basis and consider A and T as linear transformations,
            with respect to E_0, i.e.
                A(x_j) = sum_{i=0}^n a_{ij}*x_i  and  T(x_j) = sum_{i=0}^n t_{ij}*x_i .
            Then E = (y_0,...,y_n) = ( T(x_0),...,T(x_n) ) is a new basis of K[x_0,...,x_n].
            and
                LinearValuation( K[x_0,...,x_n], v_K, A, u ) = v_{E,u} .
            The graded reduction ring of v_{E,u} is the graded ring
                R = k[t,t^(-1)][Y_0,...,Y_n]
            where
                t = [p_K], Y_0 = [y_0],..., Y_n = [y_n]
            are the graded reduction of
                p_K, y_0,..., y_n,
            respectively. Moreover, the grading of R is given by the degrees of
            the generators
                |t| = v_K(p_K), |Y_0| = v_{E,u}(y_0),..., |Y_n| = v_{E,u}(y_n) .
            To compute the graded reduction  [F] of F we first have to write F with
            respect to the basis E. This will be done as follows. If we view
            E_0 = (x_0,...,x_n) as a vector in Sage, we get
                (x_0,...,x_n)*T = (sum_{i=0}^n T_{i,0}*x_i,...,sum_{i=0}^n T_{i,n}*x_i)
                                = (y_0,...,y_n)
            and therefore
                F(x_0,...,x_n) = F( (y_0,...,y_n)*T^(-1) ) = F( (y_0,...,y_n)*A ) .
            Since A = T^(-1), the polynomial
                G(y_0,...,y_n) = F( (y_0,...,y_n)*A ) in K[y_0,...,y_n]
            describes F with respect to the basis (y_0,...,y_n) and A describes the base change.
            Note that mathematically the rings
                K[x_0,...,x_n] and K[y_0,...,y_n]
            are the same as well as the polynomials F and G are the same. We only perform a change
            of coordinates. So, for a multi index set I subset NN^(n+1) let
                F = G = sum_{i in I} g_i y_0^i_0 * ... * y_n^i_n .
            Now, with the notation
                <i,u> = i_0*u_0 + ... + i_n*u_n
            and
                J = { i in I : v_{E,u}(F) = v_{E,u}(g_i y_0^i_0 * ... * y_n^i_n) },
            we have
                 v_{E,u}(F) = min( v_K(g_i) + <i,u> : i in I )
            and therefore the graded reduction [F] is given by
                sum_{i in J} [g_i / p_K^(v_{E,u}(F)-<i,u>)]*t^(v_{E,u}(F)-<i,u>) * Y_0^i_1*...*Y_n^i_n,
            where
                [g_i / p_K^(v_{E,u}(F)-<i,u>)]
            is the reduction of
                g_i / p_K^(v_{E,u}(F)-<i,u>)
            with respect to v_K.

            REMARK 1. Note that the valuations
                v1 = LinearValuation( K[x_0,...,x_n], v_K, A, u )
            and
                v2 = LinearValuation( K[y_0,...,y_n], v_K, identity_matrix, u )
            are both (mathematically) equal to v_{E,u}. As well as F and G are
            both (mathematically) equal. In particular, (mathematically) we have
                v_{E,u}(F) = v_{E,u}(G) .
            But in SageMath we have
                v1 != v2
            as well as
                F != G .
            But by construction of v1 and v2 as well as F and G we have
                v1.evaluate_at(F) = v2.evaluate_at(G) .
            Moreover, to compute
                v_{E,u}(g_j y_0^j_0 * ... * y_n^j_n)
            in SageMath, we have to compute
                v2.evaluate_at(g_j y_0^j_0 * ... * y_n^j_n),
            since the monomial
                g_j y_0^j_0 * ... * y_n^j_n
            is already written with respect to the basis E = (y_0,...,y_n).

            REMARK 2. We do not want to introduce new variables and hence we will consider
                G(x_0,...,x_n) = F( (x_0,...,x_n)*T^(-1) )
            and not G(y_0,...,y_n). But for mathematical clarity, we will still
            refer to (y_0,...,y_n) in the comments.
        """

        if not polynomial.parent() == self.polynomial_ring:
            raise TypeError

        prime_element = self.base_ring_valuation.uniformizer()

        BaseRing = LaurentPolynomialRing( self.base_ring_valuation.residue_field(), 't' )
        BaseRing_generator = BaseRing.gen()

        N = self.get_dimension()
        GradedReductionRing = PolynomialRing( BaseRing, 'y', N )

        E_matrix = identity_matrix( self.polynomial_ring.base_ring(), N )
        G = polynomial( list( vector(self.standard_basis()) * self.base_change_matrix ) )
        v_help = LinearValuation( self.polynomial_ring, self.base_ring_valuation, E_matrix, self.weight_vector )

        # By REMARK 1 we have self.evaluate_at(polynomial) = v_help.evaluate_at(G)
        polynomial_valuation = v_help.evaluate_at( G )

        f = 0
        for monom_coefficient, monom in list(G):
            if v_help.evaluate_at( monom_coefficient*monom ) == polynomial_valuation:
                prime_element_exponent = self.base_ring_valuation( monom_coefficient ) / self.base_ring_valuation( prime_element )
                monom_coefficient_reduction = self.base_ring_valuation.reduce( monom_coefficient / prime_element**prime_element_exponent )
                monomial_term_reduction = monom_coefficient_reduction
                for j, variable_degree in enumerate(monom.degrees()):
                    monomial_term_reduction = monomial_term_reduction * GradedReductionRing.gen(j)**variable_degree
                f = f + BaseRing_generator**prime_element_exponent * monomial_term_reduction

        return GradedReduction(self, f)





class GradedReduction:

    def __init__(self, linear_valuation, gr_polynomial):
        """
        Return ...

        INPUT:
            linear_valuation    - linear valuation on a polynomial ring K[x_0,...,x_n]
                                  over a field K with discrete valuation v_K
            gr_polynomial       - graded reduction of F in K[x_0,...,x_n] with respect
                                  to 'linear_valuation'

        OUTPUT:
            graded reduction of F with respect to 'linear_valuation' as object of the
            class 'GradedReduction'

        MATHEMATICS and IMPLEMENTATION:
            First, let
                f = gr_polynomial,
                s = linear_valuation.get_base_ring_valuation().value_group().gen(),
                u = linear_valuation.get_weight_vector() .
            Then s is a rational number and u = (u_0,...,u_n) a trupel of rational
            numbers. Moreover, f is the graded reduction of a polynomial
                F in K[x_0,...,x_n]
            with respect to a linear valuation v. The graded reduction ring of v is
            given by the polynomial ring
                R = k[t,t^(-1)][y_0,...,y_n]
            over the Laurent polynomial ring k[t,t^(-1)]. Here y_i is the reduction
            of x_i. Thus, we can write f as
                sum_{j,i_0,...,i_n} a_{j,i_0,...,i_n} * t^j * y_0^i_0 * ... * y_n^i_n .
            Furthermore, R has a grading given by
                |a| = 0 for all a in k, |t| = s,
            and
                |y_0| = v(x_0) = u_0, ..., |y_n| = v(x_n) = u_n .
            With respect to this grading, f is a homogeneous element. With the notation
                <i,u> = i_0 * u_0 + ... + i_n * u_n
            we have
                |f| = j*s + <i,u>
            for all j and i = (i_0,...,i_n) with nonzero a_{j,i_0,...,i_n}.
            Not let
                d_0 = (u_0/s).denominator(), ..., d_n = (u_n/s).denominator()
            and
                d = lcm(d_0,...,d_n) .
            In an algebraic closure of the fraction field of k[t,t^(-1)] let T be an
            element with
                T^d = t .
            Then, in the integral extension of rings
                k[t,t^(-1)] ---> k[T,T^(-1)]
            we can write
                T^(d/d_i) = (t^(1/d))^(d/d_i) = t^(1/d * d/d_i) = t^(1/d_i) .
            Thus, over k[T,T^(-1)] we can introduce the variables
                z_0 = t^(-u_0/s) * y_0, ..., z_n = t^(-u_n/s) * y_n .
            Note that
                |z_i| = -u_i / s * |t| + |y_i| = -u_i + u_i = 0
            for any i. Further, using
                y_i = t^(u_i/s) * z_i
            we get
                t^j*y_0^i_0*...*y_n^i_n = t^(j+1/s*<i,u>)*z_0^i_0*...*z_n^i_n .
            As explained above, we have
                |f| = j*u + <i,u>
            and therefore
                |f|/s = j + 1/s * <i,u>
            for all j and i = (i_0,...,i_n) with nonzero a_{j,i_0,...,i_n}.
            Thus, we have
                t^j*y_0^i_0*...*y_n^i_n = t^(|f|/s)*z_0^i_0*...*z_n^i_n
            and therefore f is equal to
                t^(|f|/s)*sum_{j,i_0,...,i_n} a_{j,i_0,...,i_n}*z_0^i_0*...*z_n^i_n .
            Now we view
                g = sum_{j,i_0,...,i_n} a_{j,i_0,...,i_n}*z_0^i_0*...*z_n^i_n
            as an element of the polynomial ring
                k[z_0,...,z_n]
            over the field k.
            For computational reasons we will work with g and transfere the
            results back to f with the inverse coordinate change.
        """

        self.linear_valuation = linear_valuation
        self.polynomial_ring = linear_valuation.get_polynomial_ring()
        self.gr_polynomial = gr_polynomial
        self.base_ring_grading = linear_valuation.get_base_ring_valuation().value_group().gen()
        self.GRR_grading = linear_valuation.get_weight_vector()

        self.graded_reduction_ring = gr_polynomial.parent()
        self.residue_field = self.graded_reduction_ring.base_ring().base_ring() # linear_valuation.get_base_ring_valuation.residue_field()

        self.scaled_GRR_grading = [u / self.base_ring_grading for u in self.GRR_grading]

        N = len(self.GRR_grading)
        self.reduction_ring = PolynomialRing(self.residue_field, N, 'z')
        reduction_ring_basis = self.reduction_ring.gens()
        self.r_polynomial = 0
        for multi_index, coefficient in self.gr_polynomial.dict().items():
            monom = coefficient.coefficients()[0]
            for j, degree in enumerate(multi_index):
                monom = monom * reduction_ring_basis[j]**degree
            self.r_polynomial = self.r_polynomial + monom


    def __repr__(self):
        return str(self.gr_polynomial)


    def get_GRR(self):
        return self.graded_reduction_ring


    def GRR_standard_basis(self):
        return self.graded_reduction_ring.gens()


    def get_base_ring_grading(self):
        return self.base_ring_grading


    def get_GRR_grading(self):
        return self.GRR_grading


    def RR_standard_basis(self):
        return self.reduction_ring.gens()


    def get_scaled_GRR_grading(self):
        return self.scaled_GRR_grading


    def get_gr_polynomial(self):
        return self.gr_polynomial


    def get_r_polynomial(self):
        return self.r_polynomial


    def get_initial_polynomial_ring(self):
        return self.polynomial_ring


    def get_initial_base_ring(self):
        return self.polynomial_ring.base_ring()


    def get_valuation(self):
        return self.linear_valuation


    def graded_instabilities(self, matrix_form = 'uut'):
        """
        Return a list of graded instabilities of self

        INPUT:
            matrix_form - one of the strings 'ult', 'uut', 'integral'

        MATHEMATICAL INTERPRETATION:
            First, let
                (y_0,...,y_n) = self.GRR_standard_basis(),
                (z_0,...,z_n) = self.RR_standard_basis(),
                f = self.get_gr_polynomial(),
                g = self.get_r_polynomial(),
            compare the description of the '__init__' method.
            Let (E,w) with E = (e_0,...,e_n) be an instability of g over k, i.e.
            there exists an invertible matrix
                M in GL_n(k)
            which describes the base change from E to (z_0,...,z_n). So, if we
            view E and (z_0,...,z_n) as a vectors in SageMath, we can write
                E * M = (z_0,...,z_n) .
            Note that M, as a matrix over k, desribes a k-linear transforamtion
            of degree zero, i.e. it does not change the grading.
            Now let D be the invertible diagonal matrix
                D = diag(t^(u_0/s), ..., t^(u_n/s)) in GL_3(k[T,T^(-1)]) .
            We define the basis
                E_new = (e_0_new,...,e_n_new)
            by
                E_new * D^(-1) = E .
            Note that D describes a k-linear transforamtion which changes the degree
            of homogeneous elements. Furthermore, we obtain
                E_new * D^(-1) * M * D = E * M * D = (z_0,...,z_n) * D
                                       = (y_0,...,y_n) .
            In particular, (E_new,w) is an instability of f over k[T,T^(-1)]. Note also
            that D^(-1) * M * D describes a k-linear isomorphism of degree zero and
            hence induces an isomorphism of graded rings,
                k[T,T^(-1)][y_0,...,y_n] ---> k[T,T^(-1)][y_0,...,y_n] .
            REMARK. Of course we have
                E * M * D = (z_0,...,z_n) * D = (y_0,...,y_n)
            and therefore (E,w) is already an instability of f over k[T,T^(-1)]. But the
            matrix M*D describes a k-linear transforamtion, which shifts degrees. For
            reasons of computational clarity, we want to avoid transformations that
            lead to degree shifts. Thus, we only consider transformations of the form
                D^(-1) * M * D
            as explained above. The (i,j)-th entry of such a matrix is equal to
                m_{ij} * t^{(u_j - u_i) / s} .
        """

        if len(self.RR_standard_basis()) != 3:
            raise NotImplementedError
        if matrix_form == 'integral':
            matrix_form = self.GRR_grading

        import ProjectivePlaneCurves as PPC
        reduced_curve = PPC.ProjectivePlaneCurve(self.r_polynomial)
        Graded_Instabilities = []

        for instability in reduced_curve.instabilities():
            for T in instability.base_change_matrices(matrix_form):
                Graded_Instabilities.append(GradedInstability(self, T))

        return Graded_Instabilities


    def is_semistable(self):
        if len(self.graded_instabilities()) == 0:
            return True
        return False


    def rational_graded_instabilities(self, matrix_form = 'uut'):
        """
        Return a sublist of 'self.graded_instabilities()' consisting only of
        retional graded instabilities.

        INPUT:
            matrix_form - one of the strings 'ult', 'uut', 'integral'
        """

        Rational_Graded_Instabilities = []
        for Graded_Instability in self.graded_instabilities(matrix_form):
            if Graded_Instability.is_rational():
                Rational_Graded_Instabilities.append(Graded_Instability)

        return Rational_Graded_Instabilities





class GradedInstability:

    def __init__(self, graded_reduction, instability_matrix):
        """
        INPUT:
            graded_reduction    - object in the class 'GradedReduction'
            instability_matrix  - invertible matrix over the degree zero subring of the graded reduction ring


        MATHEMATICAL INTERPRETATION:
            First, let
                f in k[t,t^(-1)][y_0,...,y_n]
            and
                g in k[z_0,...,z_n]
            as defined in the mathematical interpretation of the '__init__' method
            in the 'GradedReduction' class. Moreover, let
                M = (m_{ij}) in GL_n(k)
            and
                D = diag(t^(u_0/s), ..., t^(u_n/s)) in GL_3(k[T,T^(-1)]) .
            as in the 'instabilities' method in the 'GradedReduction' class.
            ToDo...
        """

        self.graded_reduction = graded_reduction
        self.polynomial_ring = graded_reduction.get_initial_polynomial_ring()
        self.graded_reduction_ring = graded_reduction.get_GRR()
        self.base_ring_grading = graded_reduction.get_base_ring_grading()
        self.GRR_grading = graded_reduction.get_GRR_grading()
        self.instability_matrix = instability_matrix


    def __repr__(self):
        return f"Graded Instability of {self.graded_reduction}"


    def is_rational(self):
        """
        Return True if self.instability_matrix is defined over the graded_reduction_ring
        and False otherwise.

        MATHEMATICAL INTERPRETATION:
            To understand the notation, first read the mathematical interpretation in
            the '__init__' method of self. The (i,j)-th entry of the matrix
                self.instability_matrix
            is of the form
                m_{ij} * t^{(u_j - u_i) / s}
            with m_{ij} in k. Thus, we have
                self.instability_matrix[i][j] in k[t,t^(-1)]
            if and only if
                m_{ij} = 0 or s.divides(u_j - u_i) .
        """

        for i, row in enumerate(self.instability_matrix):
            for j, entry in enumerate(row):
                if entry != 0:
                    t_exponent = (self.GRR_grading[j] - self.GRR_grading[i]) / self.base_ring_grading
                    if not t_exponent.is_integral() :
                        return False
        return True


    def get_matrix(self):
        t = self.graded_reduction_ring.base_ring().gen()
        N = len(self.GRR_grading)
        graded_trafo_matrix = [[0 for i in range(N)] for j in range(N)]

        for i, row_of_instability_matrix in enumerate(self.instability_matrix):
            for j, entry_in_row in enumerate(row_of_instability_matrix):
                if entry_in_row != 0:
                    t_exponent = (self.GRR_grading[j] - self.GRR_grading[i]) / self.base_ring_grading
                    if t_exponent.is_integral():
                        graded_trafo_matrix[i][j] = entry_in_row * t**t_exponent
                    else:
                        return None
        return matrix(self.graded_reduction_ring, graded_trafo_matrix)


    def lift_matrix(self):
        if not self.is_rational():
            return None

        base_ring = self.graded_reduction.get_initial_base_ring()
        _valuation_ = self.graded_reduction.get_valuation()
        base_ring_valuation = _valuation_.get_base_ring_valuation()
        prime_element = base_ring_valuation.uniformizer()
        N = len(self.GRR_grading)
        lifted_graded_trafo_matrix = [[0 for i in range(N)] for j in range(N)]

        for i, row_of_instability_matrix in enumerate(self.instability_matrix):
            for j, entry_in_row in enumerate(row_of_instability_matrix):
                if entry_in_row != 0:
                    exponent = (self.GRR_grading[j] - self.GRR_grading[i]) / self.base_ring_grading
                    if exponent.is_integral():
                        lifted_graded_trafo_matrix[i][j] = base_ring_valuation.lift(entry_in_row) * prime_element**exponent

        return matrix(base_ring, lifted_graded_trafo_matrix)


    def print_matrix(self):
        """
        Just to see the matrix in case if not self.is_rational(). 
        """

        if self.is_rational():
            print(self.get_matrix())
        else:
            N = len(self.GRR_grading)
            t = var('t')
            array_2d = [[0 for i in range(N)] for j in range(N)]
            for i in range(N):
                for j in range(N):
                    t_exponent = (self.GRR_grading[j] - self.GRR_grading[i]) / self.base_ring_grading
                    if self.instability_matrix[i][j] == 0:
                        array_2d[i][j] = 0
                    elif t_exponent == 0:
                        array_2d[i][j] = self.instability_matrix[i][j]
                    else:
                        array_2d[i][j] = "(" + str(self.instability_matrix[i][j]) + ")*" + str(t) + "^(" + str(t_exponent) + ")"

            v_K = self.graded_reduction.get_valuation().get_base_ring_valuation()
            print("Let t be the reduction of the uniformizer " + str(v_K.uniformizer()) + ". Then the base change to instability is given by the following matrix:")

            # Calculate the maximum width of each column
            col_widths = [max( len(str(item)) for item in col ) for col in zip(*array_2d)]

            # Print the array with aligned columns
            for row in array_2d:
                for i, item in enumerate(row):
                    print(str(item).ljust(col_widths[i] + 2), end="")  # Add spacing
                print()



# ================== helper functions ==================

def _apply_matrix(T, F, affine_patch = None):
    r"""
    Return F((x_0,...,x_n) * T) or its dehomogenization
    at affine_patch, i.e. x_{affine_patch} = 1 if
    affine_patch != None

    INPUT:
        T            - (n+1)x(n+1) matrix over K, where n+1 is the number of variables for F.
        F            - polynomial in K[x_0,...,x_n]
        affine_patch - integer between 0 and n, or None

    OUTPUT:
        Polynomial resulting from F((x_0,...,x_n) * T), where x_{affine_patch}
        is set to 1 *before* the matrix multiplication if affine_patch is not None.

    MATHEMATICAL INTERPRETATION:
        Let X = (x_0, ..., x_n) be the coordinate variables of the polynomial ring of F.
        This function computes the polynomial H(X) = F(X * T).
        If affine_patch = k is specified, it then returns H(x_0, ..., x_{k-1}, 1, x_{k+1}, ..., x_n),
        effectively performing a substitution x_k = 1 into the expression H(X).
    """
    R = F.parent()
    num_gens = R.ngens()

    if not (T.nrows() == num_gens and T.ncols() == num_gens):
        raise ValueError(f"Matrix T must be {num_gens}x{num_gens}")

    generators = list(R.gens())

    if affine_patch is not None:
        if not (0 <= affine_patch < num_gens):
            raise ValueError(f"affine_patch must be between 0 and {num_gens-1}")
        generators[affine_patch] = R(1)

    transformed_vars = list(vector(generators) * T)        
    return F(transformed_vars)

