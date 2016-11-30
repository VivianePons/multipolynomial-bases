# -*- coding: utf-8 -*-
"""
Monomial basis with type

This module implements the Monomial Basis of multivariate polynomials associated
with a group type that determines the actions on monomials.

"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.combinat.root_system.root_system import RootSystem
from sage.rings.integer import Integer
from sage.structure.unique_representation import UniqueRepresentation

from basis import PolynomialRingWithBasis, FiniteRankPolynomialRingWithBasis, MonomialKeyWrapper

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

class PolynomialRingWithBasisFromAmbientSpace(PolynomialRingWithBasis):
    r"""
    This class implements the ambient space basis. This is really close to
    the monomial basis as polynomials are also seen as sum of monomials. But
    now, each monomial is indexed by an element of the ambient space basis
    and so has a group type embeded in it.

    The class is called by
    ``multivariate_polinomials.MultivariatePolynomialAlgebra.monomial_basis_with_type``

    It is a representation of ``multivariate_polinomials.MultivariatePolynomialAlgebra``


    As the number of variable is not set, then the ambient space basis is
    not directly created, but the group type is kept.

    INPUT:

    -``abstract_polynomial_ring``, The facade abstract polynomial ring of
    which ``self`` is a representation
    - ``group_type`` the group type of the ambient space bases (``A``, ``B``,
    ``C`` or ``D``)
    - ``basis_name``, the name of the basis
    -``basis_repr`` (optional), the string representating the monomials,
     by default it is ``abstract_polynomial_ring._main_repr_var``

    EXAMPLES::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A = MultivariatePolynomialAlgebra(QQ)
        sage: ma = A.monomial_basis_with_type("A"); ma
        The Multivariate polynomial algebra on x over Rational Field on the Ambient space basis of type A

    TESTS::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A = MultivariatePolynomialAlgebra(QQ)
        sage: ma = A.monomial_basis_with_type("A")
        sage: TestSuite(ma).run()

    """

    def __init__(self, abstract_polynomial_ring, group_type, basis_name, basis_repr = None):
        r"""
        TESTS::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: ma = A.monomial_basis_with_type("A")
        """
        if(basis_repr is None): basis_repr = abstract_polynomial_ring._main_repr_var
        self._group_type = group_type
        PolynomialRingWithBasis.__init__(
            self,
            abstract_polynomial_ring,
            basis_name,
            1,
            basis_repr
        )

    def equivalent_basis(self, abstract_polynomial_ring):
        r"""
            Returns the ambien space basis of another abstract polynomial ring.

            INPUT:
            - ``abstract_polynomial_ring``, an abstract polynomial ring
            which is not the abstract polynomial ring of ``self``

            OUTPUT:
            The ambient space basis of ``abstract_polynomial_ring``

            EXAMPLES::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A = MultivariatePolynomialAlgebra(QQ)
                sage: ma = A.monomial_basis_with_type("A"); ma
                The Multivariate polynomial algebra on x over Rational Field on the Ambient space basis of type A
                sage: AZ = MultivariatePolynomialAlgebra(ZZ)
                sage: maZ = ma.equivalent_basis(AZ); maZ
                The Multivariate polynomial algebra on x over Integer Ring on the Ambient space basis of type A

        """
        return abstract_polynomial_ring.monomial_basis_with_type(self._group_type)


    def _finite_rank_basis_instance(self,n):
        r"""
        INPUT:
        - ``n`` : the number of variables

        OUTPUT:
            - the ambient space basis on ``n`` variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: ma = A.monomial_basis_with_type("A")
            sage: ma._finite_rank_basis_instance(3)
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the Monomial basis of type A

        """
        return self.abstract_algebra().algebra_finite_nb_variables(n).monomial_basis_with_type(self._group_type, self.basis_repr())


    def group_type(self):
        r"""
        Return the group type of ``self``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._group_type

    class _divided_difference_wrapper(PolynomialRingWithBasis._divided_difference_wrapper):

        def __init__(self, module, i, t1=1, t2=2):
            r"""
            TESTS::

                sage: # Fix a nice test
            """
            PolynomialRingWithBasis._divided_difference_wrapper.__init__(self,module,i)
            self._t1 = t1
            self._t2 = t2

        @cached_method
        def divided_difference_on_basis(self,key):
            r"""
            Returns ...

            EXAMPLES::

                sage: # Fix a nice example
            """
            i = self._i
            ambient_space = key.parent().ambient_space()
            as_key = key.ambient_space_element()
            n = as_key.scalar(ambient_space.simple_coroot(i))
            if n >= 0:
                return self._module.sum_of_monomials((key.parent()(as_key-(j)*ambient_space.simple_root(i)-ambient_space.basis()[i-1]) for j in xrange(n)))
            else:
                return -self.divided_difference_on_basis(key.parent()(ambient_space.simple_reflection(i)(as_key)))

        @cached_method
        def isobaric_divided_difference_on_basis(self, key):
            r"""
            Returns ...

            EXAMPLES::

                sage: # Fix a nice example
            """
            i = self._i
            ambient_space = key.parent().ambient_space()
            as_key = key.ambient_space_element()
            n = as_key.scalar(ambient_space.simple_coroot(i))
            if n >= -1:
                return self._module.sum_of_monomials((key.parent()(as_key-j*ambient_space.simple_root(i)) for j in range(n+1)))
            else:
                return -self._module.sum_of_monomials((key.parent()(as_key+(j+1)*ambient_space.simple_root(i)) for j in range(-n-1)))

        @cached_method
        def hat_isobaric_divided_difference_on_basis(self, key):
            res1 = self.isobaric_divided_difference_on_basis(key)
            res2 = self._module(key)
            return res1 - res2

        def si_on_basis(self, key):
            return self._module(key.parent()(key.ambient_space_element().weyl_action([self._i])))

        @cached_method
        def hecke_generator_on_basis(self, key):
            if(self._module.group_type() != "A"):
                raise NotImplementedError, "The hecke algebra operator is only implemented in type A"%()
            pi = self.isobaric_divided_difference_on_basis
            res1 = pi(key) * (self._t1 + self._t2)
            res2 = self._module(key)
            res2 = res2.si(self._i) * self._t2
            return res1 - res2

class FiniteRankPolynomialRingWithBasisFromAmbientSpace(FiniteRankPolynomialRingWithBasis):
    r"""
    This class implements the ambient basis on a given number of variables
    it is obtained automatically by PolynomialRingWithBasisFromAmbientSpace when a polynomial is created
    see AbastractPolynomialRing.monomial_basis_with_type for more information

    EXAMPLES::

        sage: # Fix a nice example

    TESTS::

        sage: # Fix a nice test
    """
    def __init__(self, abstract_polynomial_ring, group_code, group_type, basis_name, basis_repr = None, extra_category = None):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        if(basis_repr is None): basis_repr = abstract_polynomial_ring._main_repr_var
        self._root_system = RootSystem(group_code)
        self._group_type = group_type
        self._basis_keys = MonomialKeyWrapper(self._root_system)
        FiniteRankPolynomialRingWithBasis.__init__(
            self,
            abstract_polynomial_ring,
            abstract_polynomial_ring.polynomial_ring_tower().monomial_basis_with_type(group_type),
            self._basis_keys,
            basis_name,
            basis_repr,
            extra_category = extra_category
        )

        monomial_basis = self.abstract_algebra().monomial_basis()
        self._to_monomial_morphism = self._module_morphism(
            self._to_monomial_on_basis,
            codomain = monomial_basis
        )
        self._from_monomial_morphism = monomial_basis._module_morphism(
            self._from_monomial_on_basis,
            codomain = self
        )

        #temp#
        self._to_monomial_morphism.register_as_coercion()
        self._from_monomial_morphism.register_as_coercion()


    def product_on_basis(self, key1, key2):
        r"""
        Returns the element of ``self`` which is the product of basis
        element indexed by ``key1`` and ``key2``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self.term( self._basis_keys(tuple(key1[i]+key2[i] for i in xrange(len(key1)))) )

    def _to_monomial_on_basis(self, key):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        monomial_basis = self.abstract_algebra().monomial_basis()
        return monomial_basis(key)

    def _from_monomial_on_basis(self, key):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        return self.term( self._basis_keys( key ) )

    def to_monomial_morphism(self):
        """
        Returns the module morphism from this basis to the monomial basis

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._to_monomial_morphism

    def from_monomial_morphism(self):
        """
        Returns the module morphism from monomial basis to this basis

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._from_monomial_morphism

    def weyl_group(self):
        r"""
        Returns the weyl group acting on ``self``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._basis_keys.ambient_space().weyl_group()

    def group_type(self):
        r"""
        Returns the group type acting on ``self``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._group_type


    def __getitem__(self, c, *rest):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        if len(rest) > 0 or type(c) is int or type(c) is Integer:
            c = tuple([c])+tuple(rest)
        return self.term( self._basis_keys( list(c) ) )





    class Element(FiniteRankPolynomialRingWithBasis.Element):

        def __invert__(self):
            """
                Inverts self if self is composed by only one element. Otherwise,
                it raises a valueError exception

                EXAMPLES::

                    sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                    sage: A = MultivariatePolynomialAlgebra(QQ);
                    sage: MA = A.monomial_basis_with_type("A")
                    sage: pol = MA[2,2,3]; pol
                    xA[2, 2, 3]
                    sage: pol^-1
                    xA[-2, -2, -3]
                    sage: pol2 = MA[1,1,1]; pol2
                    xA[1, 1, 1]
                    sage: (pol+pol2)^-1
                    Traceback (most recent call last):
                    ...
                    ValueError: xA[2, 2, 3] + xA[1, 1, 1] is not invertible in The Multivariate polynomial algebra on x over Rational Field with 3 variables on the Monomial basis of type A

            """
            if( len(self) ==1 ):
                l = list(self)
                vect = l[0][0]
                coef = l[0][1]
                return coef**-1 * self.parent()( tuple(-v for v in vect) )
            raise ValueError,"%s is not invertible in %s"%(self, self.parent())

