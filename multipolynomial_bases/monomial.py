# -*- coding: utf-8 -*-
r"""
Monomial basis

This module implements the Monomial Basis of multivariate polynomials.

"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer
from sage.structure.unique_representation import UniqueRepresentation

from basis import PolynomialRingWithBasis, FiniteRankPolynomialRingWithBasis, MonomialKeyWrapper
class MonomialBasis(PolynomialRingWithBasis):
    r"""
    This class implements the monomial basis. Polynomials are seen as
    sum of monomials. The monomials are indexed by integer vectors that are
    the exponents of the variables.

    The class is called by
    ``multivariate_polinomials.MultivariatePolynomialAlgebra.monomial_basis``

    It is a representation of ``multivariate_polinomials.MultivariatePolynomialAlgebra``

    The number of variables is not set, this class is a facade for
    ``FiniteMonomialBasis``.

    INPUT:

    -``abstract_polynomial_ring``, The facade abstract polynomial ring of
    which ``self`` is a representation
    -``basis_repr`` (optional), the string representating the monomials,
     by default it is ``abstract_polynomial_ring._main_repr_var``


    EXAMPLES::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
        sage: x
        The Multivariate polynomial algebra on x over Rational Field on the monomial basis
        sage: x.an_element()
        2*x[1, 0, 0] + 3*x[0, 1, 0] + x[0, 0, 0] + x[1, 2, 3]
        sage: pol = x[2,2,1] + x[3,2]; pol
        x[2, 2, 1] + x[3, 2, 0]



    TESTS::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A = MultivariatePolynomialAlgebra(QQ)
        sage: m = A.monomial_basis()
        sage: TestSuite(m).run()

    """
    def __init__(self, abstract_polynomial_ring, basis_repr = None):
        r"""
        TESTS::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: x = A.monomial_basis()
        """
        if(basis_repr is None): basis_repr = abstract_polynomial_ring._main_repr_var
        PolynomialRingWithBasis.__init__(
            self,
            abstract_polynomial_ring,
            "monomial basis",
            1,
            basis_repr
        )

    def equivalent_basis(self, abstract_polynomial_ring):
        r"""
            Returns the monomial basis of another abstract polynomial ring.

            INPUT:
            - ``abstract_polynomial_ring``, an abstract polynomial ring
            which is not the abstract polynomial ring of ``self``

            OUTPUT:
            The monomial basis of ``abstract_polynomial_ring``

            EXAMPLES::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis();m
                The Multivariate polynomial algebra on x over Rational Field on the monomial basis
                sage: B = MultivariatePolynomialAlgebra(ZZ)
                sage: mb = m.equivalent_basis(B); mb
                The Multivariate polynomial algebra on x over Integer Ring on the monomial basis

        """
        return abstract_polynomial_ring.monomial_basis()

    def _finite_rank_basis_instance(self, n):
        r"""
        INPUT:
        - ``n`` : the number of variables

        OUTPUT:
            - the monomial basis on ``n`` variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: m = A.monomial_basis()
            sage: m._finite_rank_basis_instance(3)
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis

        """
        return self.abstract_algebra().algebra_finite_nb_variables(n).monomial_basis(self.basis_repr())


    def _to_expr_on_basis(self, vector, alphabet = None):
        r"""
        Transforms a vector key of a finite basis into a symbolic
        expression factor.

        It uses formal variables with names : basis_repr + i
        (by default, 'x1' to 'xn')

        INPUT :
            - ``vector`` a key element of a finite basis
            ``CominatorialFreeModule``
            - ``alphabet`` an obtional alphabet of formal variables, if let
            to ``None``, it uses basis_repr + i (by default, 'x1' to 'xn')

        OUTPUT:
            - a symbolic expression corresponding to the monomial

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: pol = x[2,2,1] + x[3,2]; pol
            x[2, 2, 1] + x[3, 2, 0]
            sage: key = list(pol)[0][0]; key
            [2, 2, 1]
            sage: x._to_expr_on_basis(key)
            x1^2*x2^2*x3
            sage: var('a,b,c')
            (a, b, c)
            sage: alphabet = [a,b,c]
            sage: x._to_expr_on_basis(key, alphabet=alphabet)
            a^2*b^2*c

        """
        from sage.calculus.var import var
        key = list(vector)
        if alphabet is None:
            alphabet = [var(self.basis_repr() + str(i)) for i in xrange(1,len(key)+1)]
        expr = 1
        for i in xrange(len(key)):
            expr*=alphabet[i]**key[i]

        return expr

    def _maxDiffDiv_on_basis(self, key):
        r"""
        Apply the maximal divided difference on a monomial key.

        The algorithm just sorts the key and counts the number of inversions.
        Then the result is `(-1)^{inv}Schub_{sortedKey - [0, ...,n-1]}`.

        The result is given as a Schubert poynomial which is symmetric and
        can be understood as a symmetric function.

        This method is then wrapped in a module morphism to give the
        maximal divided difference of a polynomial.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: m = A.monomial_basis()
            sage: key = list(m[1,2,3])[0][0];key
            [1, 2, 3]
            sage: m._maxDiffDiv_on_basis(key)
            -Y[1, 1, 1]

        """
        result_key = list(key)
        Schub = self.abstract_algebra().schubert_basis("A")
        sign = -1
        l = len(result_key) -1
        flag = True
        while flag:
            flag = False
            for i in xrange(l):
                if(result_key[i] == result_key[i+1]): return Schub(0)
                if(result_key[i] > result_key[i+1]):
                    result_key[i], result_key[i+1] = result_key[i+1], result_key[i]
                    flag = True
                    sign*=-1
            l = l-1
        for i in xrange(len(result_key)):
            result_key[i]-=i
        return sign * Schub( result_key )

    class _divided_difference_wrapper(PolynomialRingWithBasis._divided_difference_wrapper):
        r"""
        This inner class is used to wrap the divided difference on basis
        methods. The morphism module on basis method of
        ``CombinatorialFreeModule`` only takes the element key as a
        parameter. But we also need to know the type, number of variables
        and number of the divided difference.

        """
        def __init__(self, module, i, otype = "A", t1=1, t2=2):
            r"""
            INPUT:
            - ``module``, a ``CombinatorialFreeModule``, the
            ``FiniteMonomialBasis`` on which the morphism apply
            - ``i`` the number of the wanted divided difference
            - ``otype`` the type of wanted divided difference, default is
            ``A``
            - ``t1`` and ``t2``, parameters for the hecke algebra generator

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
            """
            PolynomialRingWithBasis._divided_difference_wrapper.__init__(self,module,i)
            self._type = otype
            self._t1 = t1
            self._t2 = t2

        def _generic_on_basis(self, key, method):
            r"""
            Generic code to apply the on basis divided difference.

            The actual divided difference is done in
            :class: ``ambient_space_basis.PolynomialRingWithBasisFromAmbientSpace``

            This method computes the needed divided difference on the rigth
            number of variables and gives the result in the monomial basis.

            As an example, if the second type ``B`` divided difference is
            needed on a 3 variables polynomial. The ambient space basis of
            size 2 is created and applies the divided difference on the two
            first variables of the polynomial The third variables is added
            afterward and the polynomial is reconverted into the monomial
            basis.

            INPUT:

            - ``key``, the element key to apply the divided difference
            - ``method`` the name of the method to use in the ambient space
            basis.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
                sage: key = list(m[1,2,1])[0][0];key
                [1, 2, 1]
                sage: wrapp._generic_on_basis(key,"divided_difference")
                -x[1, 1, 1]

            """
            i = self._i
            otype  = self._type
            if(otype == "A"):
                monomial_basis_with_type = self._module.abstract_algebra().monomial_basis_with_type(otype)
                morph = monomial_basis_with_type.get_morphism(i,method=method,t1=self._t1,t2 = self._t2)
                return self._module( morph(monomial_basis_with_type( key )))

            monomial_basis_with_type = self._module.basis_tower().abstract_algebra().algebra_finite_nb_variables(i).monomial_basis_with_type(otype)
            morph = monomial_basis_with_type.get_morphism(i,method=method,t1=self._t1,t2 = self._t2)
            key1 = key[:i]
            key2 = [0 for j in xrange(i)]
            key2.extend(key[i:])
            return self._module( morph(monomial_basis_with_type( key1 )).change_nb_variables(len(key))) * self._module( key2 )

        def divided_difference_on_basis(self,key):
            r"""
            On basis method for the divided difference, see
            ``FiniteMonomialBasis.divided_difference_morphism`` method
            for more details on what it does.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
                sage: key = list(m[1,2,1])[0][0];key
                [1, 2, 1]
                sage: wrapp.divided_difference_on_basis(key)
                -x[1, 1, 1]

            """
            return self._generic_on_basis(key, "divided_difference")

        def isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis method for the isobar divided difference, see
            ``FiniteMonomialBasis.divided_difference_isobar_morphism`` method
            for more details on what it does.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
                sage: key = list(m[2,1,1])[0][0];key
                [2, 1, 1]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                x[1, 2, 1] + x[2, 1, 1]

            """
            return self._generic_on_basis(key, "isobaric_divided_difference")

        def hat_isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis method for the isobar divided difference, see
            ``FiniteMonomialBasis.divided_difference_isobar_hat_morphism`` method
            for more details on what it does.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
                sage: key = list(m[2,1,1])[0][0];key
                [2, 1, 1]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                x[1, 2, 1]

            """
            return self._generic_on_basis(key, "hat_isobaric_divided_difference")

        def hecke_generator_on_basis(self, key):
            r"""
            On basis method for the hecke generator, see
            ``FiniteMonomialBasis.hecke_generator_morphism`` method
            for more details on what it does.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
                sage: key = list(m[1,2,1])[0][0];key
                [1, 2, 1]
                sage: wrapp.hecke_generator_on_basis(key)
                -2*x[2, 1, 1]

            """
            return self._generic_on_basis(key, "hecke_generator")

        def product_variable_on_basis(self, key):
            r"""
            On basis method for a product by a variable.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
                sage: key = list(m[2,1,1])[0][0];key
                [2, 1, 1]
                sage: wrapp.product_variable_on_basis(key)
                x[3, 1, 1]
            """
            key = list(key)
            i = self._i
            key[i-1] = key[i-1]+1
            return self._module(key)

        def si_on_basis(self, key):
            r"""
            On basis method for the si action

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: m3 = m.finite_rank_basis(3)
                sage: wrapp = m._divided_difference_wrapper(m3,1,"A")
                sage: key = list(m[2,1,1])[0][0];key
                [2, 1, 1]
                sage: wrapp.si_on_basis(key)
                x[1, 2, 1]
                sage: wrapp = m._divided_difference_wrapper(m3,1,"B")
                sage: wrapp.si_on_basis(key)
                x[-2, 1, 1]
                sage: wrapp = m._divided_difference_wrapper(m3,1,"C")
                sage: wrapp.si_on_basis(key)
                x[-2, 1, 1]
                sage: wrapp = m._divided_difference_wrapper(m3,2,"D")
                sage: wrapp.si_on_basis(key)
                x[-1, -2, 1]
                sage: wrapp = m._divided_difference_wrapper(m3,2,"E")
                sage: wrapp.si_on_basis(key)
                Traceback (most recent call last):
                ...
                ValueError: Unknown type E

            """
            key = list(key)
            i = self._i -1
            otype = self._type
            if(otype=="A"):
                key[i],key[i+1] = key[i+1], key[i]
            elif(otype=="B" or otype == "C"):
                key[i] = -key[i]
            elif(otype=="D"):
                key[i],key[i-1] = -key[i-1],-key[i]
            else:
                raise ValueError, "Unknown type %s"%(otype)
            return self._module(key)

class FiniteMonomialBasis(FiniteRankPolynomialRingWithBasis):
    r"""
    This class implements the monomial basis on a given number of variables.
    It is obtained automatically by MonomialBasis when a polynomial is created,
    see :class: MonomialBasis for more details.

    It is a realization of both the monomial basis on an unset number of
    variables and the finite abstract polynomial ring on n variables.

    INPUT:
    - ``abstract_polynomial_ring`` the realization of Abstract polynomial
    ring on ``n`` variables. The number of variables is contained in this
    argument.
    -``basis_repr`` (optional), the string representating the monomials,
     by default it is ``abstract_polynomial_ring._main_repr_var``

    EXAMPLES::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A = MultivariatePolynomialAlgebra(QQ)
        sage: m = A.monomial_basis()
        sage: m3 = m.finite_rank_basis(3)
        sage: m3
        The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis


    TESTS::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A = MultivariatePolynomialAlgebra(QQ)
        sage: m = A.monomial_basis()
        sage: m3 = m.finite_rank_basis(3)
        sage: m3
        The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
        sage: TestSuite(m3).run()
    """
    def __init__(self, abstract_algebra, basis_repr = None, extra_category = None):
        r"""
        TESTS::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: m = A.monomial_basis()
            sage: m3 = m.finite_rank_basis(3)

        """
        if(basis_repr is None): basis_repr = abstract_algebra._main_repr_var
        basis_keys = MonomialKeyWrapper(abstract_algebra.nb_variables())
        FiniteRankPolynomialRingWithBasis.__init__(
            self,
            abstract_algebra,
            abstract_algebra.polynomial_ring_tower().monomial_basis(),
            basis_keys,
            "monomial basis",
            basis_repr,
            extra_category = extra_category
        )

    def __getitem__(self, c, *rest):
        r"""
        Allows the creation of elements with [ ].

        TESTS::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: m = A.monomial_basis()
            sage: m3 = m.finite_rank_basis(3)
            sage: m3[2,1,3]
            x[2, 1, 3]

        """
        if len(rest) > 0 or type(c) is int or type(c) is Integer:
            c = tuple([c])+tuple(rest)
        return self.term( self._basis_keys( c ) )


    def product_on_basis(self, vector1, vector2):
        r"""
        Return the element of ``self`` which is the product of
        basis element indexed by ``vector1`` and ``vector2``.

        This method is wrapped by ``CombinatorialFreeModule`` to compute
        the product

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: m = A.monomial_basis()
            sage: m3 = m.finite_rank_basis(3)
            sage: m3[2,3,1] * m3[1,1,2]
            x[3, 4, 3]

        """
        return self.term( self._basis_keys(tuple(vector1[i]+vector2[i] for i in xrange(len(vector1)))) )


    class Element(FiniteRankPolynomialRingWithBasis.Element):


        def __invert__(self):
            r"""
                Inverts ``self`` if ``self`` contains only one monomial.
                 Otherwise, it raises a ``ValueError`` exception

                EXAMPLES::

                    sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                    sage: A = MultivariatePolynomialAlgebra(QQ);
                    sage: M = A.monomial_basis()
                    sage: pol = M[1,2,3];pol
                    x[1, 2, 3]
                    sage: pol^-1
                    x[-1, -2, -3]
                    sage: pol2 = M[1,1,1];pol2
                    x[1, 1, 1]
                    sage: (pol+pol2)^-1
                    Traceback (most recent call last):
                    ...
                    ValueError: x[1, 1, 1] + x[1, 2, 3] is not invertible in The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis

            """
            if( len(self) ==1 ):
                l = list(self)
                vect = l[0][0]
                coef = l[0][1]
                return coef**-1 * self.parent()( tuple(-v for v in vect) )
            raise ValueError,"%s is not invertible in %s"%(self, self.parent())

        def _right_number_of_variables_pol(self, i, otype= "A"):
            r"""
            Gives the polynomial in a big enought number of variables to
            compute an action ``i`` of type ``otype``

            INPUT:
            -``i``, the number of the action
            - ``otype``, the type of the action

            OUTPUT:
            - ``self`` in a big enought number of variables to
            compute an action ``i`` of type ``otype``

            EXAMPLES::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: m = A.monomial_basis()
                sage: pol = m[2,1,1]
                sage: pol._right_number_of_variables_pol(3,"A")
                x[2, 1, 1, 0]
                sage: pol._right_number_of_variables_pol(4,"B")
                x[2, 1, 1, 0]

            """
            min = self.parent().basis_tower()._right_number_of_variables(i, otype)
            if(self.nb_variables() < min):
                return self.change_nb_variables(min)
            else:
                return self




