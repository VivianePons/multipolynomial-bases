# -*- coding: utf-8 -*-
r"""
This module implements the ring of abstract polynomials over a double set of variables
This ring is actually an ``MultivariatePolynomialAlgebra`` over antoher ``MultivariatePolynomialAlgebra``
but the module provides methods that are specific to double multivariate polynomials to make it easy
to use

EXAMPLES::

    sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
    sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ); A
    The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
    sage: A.an_element()
    y[0]*x[0, 0, 0] + 3*y[0]*x[0, 1, 0] + 2*y[0]*x[1, 0, 0] + y[0]*x[1, 2, 3]

``x`` and ``y`` correspond to the monomial basis on the two sets of variables::

    sage: x
    The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field on the monomial basis
    sage: y
    The Multivariate polynomial algebra on y over Rational Field on the monomial basis
    sage: x.an_element()
    y[0]*x[0, 0, 0] + 3*y[0]*x[0, 1, 0] + 2*y[0]*x[1, 0, 0] + y[0]*x[1, 2, 3]
    sage: y.an_element()
    y[0, 0, 0] + 3*y[0, 1, 0] + 2*y[1, 0, 0] + y[1, 2, 3]

By default, all actions are done on the ``x`` variable set::

    sage: pol = x[1,2,3] + x[2,2,4]*y[1,2,3]; pol
    y[0]*x[1, 2, 3] + (y[1,2,3])*x[2, 2, 4]
    sage: pol.divided_difference(1)
    (-y[0])*x[1, 1, 3]


You can obtain the polynomial ring on the ``y`` variables which correspond
to the coeff ring::

    sage: Coeffs = A.coeffs_ring(); Coeffs
    The Multivariate polynomial algebra on y over Rational Field

You can change the bases of the ``x`` polynomial by the usual coercion system::

    sage: Schub = A.schubert_basis(); Schub
    The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field on the Schubert basis of type A
    sage: Schub(A.an_element())
    y[0]*Yx[0, 0, 0] + 3*y[0]*Yx[0, 1, 0] + (-y[0])*Yx[1, 0, 0] + y[0]*Yx[1, 2, 3] + (-y[0])*Yx[1, 3, 2] + (-y[0])*Yx[2, 1, 3] + y[0]*Yx[2, 3, 1] + y[0]*Yx[3, 1, 2] + (-y[0])*Yx[3, 2, 1] + y[0]*Yx[4, 1, 1]

You can also change the base for the ``y``::

    sage: YSchub = Coeffs.schubert_basis(); YSchub
    The Multivariate polynomial algebra on y over Rational Field on the Schubert basis of type A
    sage: pol = y[1,2,3] * x[2,2,4]; pol
    (y[1,2,3])*x[2, 2, 4]
    sage: pol.change_coeffs_bases(YSchub)
    (Yy[1,2,3]-Yy[1,3,2]-Yy[2,1,3]+Yy[2,3,1]+Yy[3,1,2]-Yy[3,2,1]+Yy[4,1,1])*x[2, 2, 4]


Also, you can obtain the ring where the roles of `x` and `y` are exchanged::

    sage: A2 =  A.inversed_ring(); A2
    The Multivariate polynomial algebra on y over The Multivariate polynomial algebra on x over Rational Field
    sage: A2.an_element()
    x[0]*y[0, 0, 0] + 3*x[0]*y[0, 1, 0] + 2*x[0]*y[1, 0, 0] + x[0]*y[1, 2, 3]

There is a coercion between ``A`` and ``A2``::

    sage: pol = x[2,2,4]*y[1,2,3]; pol
    (y[1,2,3])*x[2, 2, 4]
    sage: A2(pol)
    (x[2,2,4])*y[1, 2, 3]


But this coercion doesn't allow operations including polynomials from ``A`` and ``A2`` as the coercion
only exists between abstract polynomial rings but not between the concrete bases : sage coercion system
doesn't look for a parent where the coercion could be made::

    sage: A.an_element() + A2.an_element()
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +: 'The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field with 3 variables on the monomial basis' and 'The Multivariate polynomial algebra on y over The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis'
    sage: A2(A.an_element()) + A2.an_element()
    (2*x[0,0,0]+3*x[0,1,0]+2*x[1,0,0]+x[1,2,3])*y[0, 0, 0] + 3*x[0]*y[0, 1, 0] + 2*x[0]*y[1, 0, 0] + x[0]*y[1, 2, 3]

Some special bases have been implemented for double polynomials::

    sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
    sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
    sage: DSchub = A.double_schubert_basis()
    sage: DSchub[0,1]
    y[0]*YY[0, 1]
    sage: DSchub[0,1].expand()
    (-yA[0,1]-yA[1,0])*xA[0, 0] + y[0]*xA[0, 1] + y[0]*xA[1, 0]
    sage: DSchub[0,1]^2
    (yA[0,0,1]-yA[0,1,0])*YY[0, 1] + y[0]*YY[0, 2] + y[0]*YY[1, 1]
    sage: DSchub(x[0,1])
    (yA[0,1])*YY[0, 0] + y[0]*YY[0, 1] + (-y[0])*YY[1, 0]
    sage: DSchub(y[0,1]*x[0,1])
     (yA[0,2])*YY[0, 0] + (y[0,1])*YY[0, 1] + (-y[0,1])*YY[1, 0]
    sage: DSchub(y[0,1])
    (y[0,1])*YY[0, 0]
    sage: DGroth = A.double_grothendieck_basis()
    sage: DGroth
    The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field on the Double Grothendieck basis of type A
    sage: DGroth[0,1]
    y[0]*GG[0, 1]
    sage: DGroth[0,1].expand()
    (-yA[1,1])*xA[-1, -1] + y[0]*xA[0, 0]
    sage: DGroth[3,2,3].isobaric_divided_difference(1)
    y[0]*GG[2, 2, 3]
    sage: DGroth[3,2,3].hat_isobaric_divided_difference(1)
    y[0]*GG[2, 2, 3] + (-y[0])*GG[3, 2, 3]
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.category import Category
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from .multivariate_polynomials import MultivariatePolynomialAlgebra, MultivariatePolynomialAlgebra_generic
from .basis import PolynomialRingWithBasis
from .linear_basis_on_vectors import LinearBasisOnVectors, FiniteLinearBasisOnVectors


class DoubleMultivariatePolynomialAlgebra_generic(MultivariatePolynomialAlgebra_generic):
    r"""
    Implementation for double polynomial algebra

    INPUT:
        - ``R``: the base ring of the algebra
        - ``repr_var1``, a string representing the main variable set, default is `x`
        - ``repr_var2``, a string representing the secondary variable set, default is `y`
        - ``inversed_ring``, the ring where the roles of the two sets of variables are inversed.
          By default, nothing is sent and the ring is created, the field is then used to avoid infinite
          recursion

    OUTPUT:

        - The abstract ring of multivariate polynomials on ``repr_var1`` over te abstract ring
          of multivariate polynomials on ``repr_var2`` over ``R``

    TESTS::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
        sage: TestSuite(A).run()
    """

    def __init__(self, R, repr_var1 = 'x', repr_var2 = 'y', inversed_ring = None):
        self._coeffs_ring = MultivariatePolynomialAlgebra_generic(R, repr_var2, always_show_main_var = True)
        MultivariatePolynomialAlgebra_generic.__init__(
            self,
            self._coeffs_ring,
            repr_var1,
            always_show_main_var = True,
            extra_bases_category = Finite_rank_double_bases()
        )
        self._repr_var1 = repr_var1
        self._repr_var2 = repr_var2
        self._coeffs_base_ring = R
        if(inversed_ring is None):
            self._inversed_ring = DoubleMultivariatePolynomialAlgebra_generic(R, repr_var2, repr_var1, self)
        else:
            self._inversed_ring = inversed_ring

        m = SetMorphism( Hom(self, self.inversed_ring()), lambda x : x.swap_coeffs_elements())
        m.register_as_coercion()



    def coeffs_ring(self):
        """r
        returns the multivariate polynomial ring on the second set of variables
        used as coefficients of the main ring on the first set of variables

        OUPUT:

            - the ring of multivariate polynomials on the second set of variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ); A
            The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
            sage: A.coeffs_ring()
            The Multivariate polynomial algebra on y over Rational Field

        """
        return self._coeffs_ring

    def inversed_ring(self):
        """r
        returns the ring of multivariate polynomials where the roles of the two sets of variables
        are exchanged

        OUTPUT:

            - the ring of multivariate polynomials on the second set of variables
              over the ring of multivariate polynomials on the first set of variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ); A
            The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
            sage: A2 = A.inversed_ring(); A2
            The Multivariate polynomial algebra on y over The Multivariate polynomial algebra on x over Rational Field

            There is a coercion between `D` and `D2`::

            sage: pol = x.an_element() * y.an_element(); pol
            (y[0,0,0]+3*y[0,1,0]+2*y[1,0,0]+y[1,2,3])*x[0, 0, 0] + (3*y[0,0,0]+9*y[0,1,0]+6*y[1,0,0]+3*y[1,2,3])*x[0, 1, 0] + (2*y[0,0,0]+6*y[0,1,0]+4*y[1,0,0]+2*y[1,2,3])*x[1, 0, 0] + (y[0,0,0]+3*y[0,1,0]+2*y[1,0,0]+y[1,2,3])*x[1, 2, 3]
            sage: A2(pol)
            (x[0,0,0]+3*x[0,1,0]+2*x[1,0,0]+x[1,2,3])*y[0, 0, 0] + (3*x[0,0,0]+9*x[0,1,0]+6*x[1,0,0]+3*x[1,2,3])*y[0, 1, 0] + (2*x[0,0,0]+6*x[0,1,0]+4*x[1,0,0]+2*x[1,2,3])*y[1, 0, 0] + (x[0,0,0]+3*x[0,1,0]+2*x[1,0,0]+x[1,2,3])*y[1, 2, 3]

        But this coercion doesn't allow for operations including polynomials from both ``A`` and ``A2`` as the coercion
        only exists between abstract polynomial rings but not between the concrete bases : sage coercion system
        doesn't look for a parent where the coercion could be made
        """
        return self._inversed_ring

    def gens(self):
        r"""
        Return a tuple whose entries are the generators for this
        object.

        In the case of the multivariate polynomial algebra, the number
        of actual generators is potentatially infinite. So this method
        actually return a tuple containing the monomial basis which can
        be seen as a multivariate gerator.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ); A
            The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
            sage: x
            The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field on the monomial basis
            sage: y
            The Multivariate polynomial algebra on y over Rational Field on the monomial basis
            sage: A.gens()
            [The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field on the monomial basis,
             The Multivariate polynomial algebra on y over Rational Field on the monomial basis]
            sage: A.gens()[0] == x
            True
            sage: A.gens()[1] == y
            True

        """
        return [self.monomial_basis(), self.coeffs_ring().monomial_basis()]

    def double_schubert_basis(self, group_type = "A", basis_name = None, basis_repr = "YY"):
        if(basis_name is None):
            basis_name = "Double Schubert basis of type " +group_type
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return DoubleSchubertBasis(self, monomial_basis_with_type, basis_name, basis_repr)

    def double_grothendieck_basis(self, group_type = "A", basis_name = None, basis_repr = "GG"):
        if(basis_name is None):
            basis_name = "Double Grothendieck basis of type " +group_type
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return DoubleGrothendieckBasis(self, monomial_basis_with_type, basis_name, basis_repr)

class Finite_rank_double_bases(Category):
    """r
    This category is an extra category for bases of double multivariate
    polynomials. It is used to add extra methods on double polynomials.

    OUTPUT:
        The category of finite rank bases for double polynomials

    EXAMPLES::
        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
        sage: A._extra_bases_category
        The category of finite rank bases for double polynomials
        sage: A.an_element().parent().category()
        Join of Category of graded algebras with basis over The Multivariate polynomial algebra on y over Rational Field and The category of finite rank bases for double polynomials and Category of realizations of The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field with 3 variables

    TESTS::

        sage: from multipolynomial_bases.double_multivariate_polynomials import Finite_rank_double_bases
        sage: C = Finite_rank_double_bases()
        sage: C
        The category of finite rank bases for double polynomials
        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
        sage: TestSuite(A.an_element().parent().category()).run()

    """

    def _repr_(self):
        """
        TESTS::

            sage: from multipolynomial_bases.double_multivariate_polynomials import Finite_rank_double_bases
            sage: Finite_rank_double_bases()
            The category of finite rank bases for double polynomials
        """
        return "The category of finite rank bases for double polynomials"

    def super_categories(self):
        """
        TESTS::

            sage: from multipolynomial_bases.double_multivariate_polynomials import Finite_rank_double_bases
            sage: Finite_rank_double_bases().super_categories()
            []
        """
        return []

    class ElementMethods:

        def change_coeffs_bases(self, new_base):
            """r
            This method changes the base of the coefficients of a given polynomial which are polynomials
            on the second set of variables

            INPUT:
                - ``new_base`` : a base of :class:`MultivariatePolynomialAlgebra` on the second set of variables

            OUTPUT:
                - the polynomial where the base of the coefficients has been changed

            EXAMPLES::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
                sage: AY = A.coeffs_ring()
                sage: pol = x[1,0,0] * y[1,2,3]; pol
                (y[1,2,3])*x[1, 0, 0]
                sage: pol.change_coeffs_bases(AY.schubert_basis())
                (Yy[1,2,3]-Yy[1,3,2]-Yy[2,1,3]+Yy[2,3,1]+Yy[3,1,2]-Yy[3,2,1]+Yy[4,1,1])*x[1, 0, 0]

            """
            return sum( [ self.parent().term(x,new_base(y)) for (x,y) in self ] )

        def inversed_ring(self):
            """r
            Return the ring of polynomials where coefficients and values
            have been inverted.

            EXAMPLES::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ); A
                The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
                sage: A.an_element().inversed_ring()
                The Multivariate polynomial algebra on y over The Multivariate polynomial algebra on x over Rational Field
            """
            return self.parent().abstract_algebra().polynomial_ring_tower().inversed_ring()

        def swap_coeffs_elements(self):
            """
            Return a double polynomial whose coefficients and values have been
            swapped.

            EXAMPLES::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
                sage: pol = A.an_element(); pol
                y[0]*x[0, 0, 0] + 3*y[0]*x[0, 1, 0] + 2*y[0]*x[1, 0, 0] + y[0]*x[1, 2, 3]
                sage: pol.swap_coeffs_elements()
                (x[0,0,0]+3*x[0,1,0]+2*x[1,0,0]+x[1,2,3])*y[0]
            """
            inversed_ring = self.inversed_ring()
            coeffs_ring = inversed_ring.coeffs_ring()
            base_x = self.parent().basis_tower().equivalent_basis(coeffs_ring)
            base_y = None
            default_base_y = None
            result = None
            for (x,y) in self:
                if(base_y is None):
                    default_base_y = y.parent().basis_tower()
                    base_y = default_base_y.equivalent_basis(inversed_ring)
                    result = base_y.zero()
                y = default_base_y( y )
                for (key, coeff) in y:
                    result += base_y.term( key, base_x.term(x,coeff))
            return result



class DoubleSchubertBasis(LinearBasisOnVectors):
    r"""
    Explain this class

    EXAMPLES::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ); A
        The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
        sage: DSchub = A.double_schubert_basis(); DSchub
        The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field on the Double Schubert basis of type A
        sage: DSchub[0,1]
        y[0]*YY[0, 1]
        sage: DSchub[0,1].expand()
        (-yA[0,1]-yA[1,0])*xA[0, 0] + y[0]*xA[0, 1] + y[0]*xA[1, 0]
        sage: DSchub[0,1]^2
        (yA[0,0,1]-yA[0,1,0])*YY[0, 1] + y[0]*YY[0, 2] + y[0]*YY[1, 1]
        sage: DSchub(x[0,1])
        (yA[0,1])*YY[0, 0] + y[0]*YY[0, 1] + (-y[0])*YY[1, 0]
        sage: DSchub(y[0,1]*x[0,1])
        (yA[0,2])*YY[0, 0] + (y[0,1])*YY[0, 1] + (-y[0,1])*YY[1, 0]
        sage: DSchub(y[0,1])
        (y[0,1])*YY[0, 0]

    """

    def __init__(self, abstract_polynomial_ring, monomial_basis_with_type, basis_name, basis_repr):
        LinearBasisOnVectors.__init__(
            self,
            abstract_polynomial_ring,
            monomial_basis_with_type,
            basis_name,
            basis_repr,
            on_basis_method = self.on_basis_method,
            variables_auto_coerce = True,
            cmp = self.cmp
        )

    def cmp(self, key1, key2):
        l = len(key1)
        for i in range(l-1,-1,-1):
            if (key1[i]>key2[i]):
                return 1
            if (key1[i]<key2[i]):
                return -1
        return 0

    def on_basis_method(self, x, basis, call_back):
        for i in range( len( x ) - 1 ):
            if( x[i]<x[i+1] ):
                x[i], x[i+1] = x[i+1]+1, x[i]
                return call_back(x).divided_difference(i+1)
        basex = basis
        basey = basex.basis_tower().equivalent_basis(basex.base_ring())
        return prod( [basex.var(i+1) - basey.var(j+1) for i in range(len(x)) for j in range(x[i])], basex.one())



class DoubleGrothendieckBasis(LinearBasisOnVectors):

    r"""
    Explain this class

    EXAMPLES::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ); A
        The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
        sage: DGroth = A.double_grothendieck_basis()
        sage: DGroth
        The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field on the Double Grothendieck basis of type A
        sage: DGroth[0,1]
        y[0]*GG[0, 1]
        sage: DGroth[0,1].expand()
        (-yA[1,1])*xA[-1, -1] + y[0]*xA[0, 0]
        sage: DGroth[3,2,3].isobaric_divided_difference(1)
        y[0]*GG[2, 2, 3]
        sage: DGroth[3,2,3].hat_isobaric_divided_difference(1)
        y[0]*GG[2, 2, 3] + (-y[0])*GG[3, 2, 3]


    """

    def __init__(self, abstract_polynomial_ring, monomial_basis_with_type, basis_name, basis_repr):

        LinearBasisOnVectors.__init__(
            self,
            abstract_polynomial_ring,
            monomial_basis_with_type,
            basis_name,
            basis_repr,
            on_basis_method = self.on_basis_method,
            variables_auto_coerce = True,
            triangular = None
        )


    def on_basis_method(self, x, basis, call_back):
        for i in range( len( x ) - 1 ):
            if( x[i]<x[i+1] ):
                x[i], x[i+1] = x[i+1]+1, x[i]
                return call_back(x).isobaric_divided_difference(i+1)
        prod = basis.one()
        basex = basis
        basey = basex.basis_tower().equivalent_basis(basex.base_ring())
        for i in range(len(x)):
            inv_x_i = basex.var(i+1)**(-1)
            for j in range(x[i]):
                prod *= (basis.one() - basey.var(j+1) * inv_x_i)
        return prod

    class _divided_difference_wrapper(PolynomialRingWithBasis._divided_difference_wrapper):
        r"""
        This class is a wrapper for the on basis isobaric divided differences
        methods.
        It contains a optimized version of the divided difference
        for Grothendieck polynomials.

        TESTS::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
            sage: DG = A.double_grothendieck_basis()
            sage: DG3 = DG.finite_rank_basis(3)
            sage: wrapp = DG._divided_difference_wrapper(DG3,1)


        """
        def isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the isobaric divided difference on the
            Grothendieck basis.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
                sage: DG = A.double_grothendieck_basis()
                sage: DG3 = DG.finite_rank_basis(3)
                sage: wrapp = DG._divided_difference_wrapper(DG3,1)
                sage: pol = DG[3,2,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                y[0]*GG[2, 2, 3]
                sage: pol = DG[2,3,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.isobaric_divided_difference_on_basis(key)
                y[0]*GG[2, 3, 3]

            Test consistency..::

                sage: pol = DG[3,2,3]
                sage: a = pol.isobaric_divided_difference(1).expand()
                sage: b = pol.expand().isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.isobaric_divided_difference(2).expand()
                sage: b = pol.expand().isobaric_divided_difference(2)
                sage: a == b
                True
                sage: pol = DG[2,3,3]
                sage: a = pol.isobaric_divided_difference(1).expand()
                sage: b = pol.expand().isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.isobaric_divided_difference(2).expand()
                sage: b = pol.expand().isobaric_divided_difference(2)
                sage: a == b
                True
            """
            i = self._i
            if(key[i-1] > key[i]):
                key2 = [key[j] for j in range(self._module.nb_variables())]
                key2[i-1], key2[i] = key2[i], key2[i-1]-1
                return self._module(key2)
            else:
                return self._module(key)

        def hat_isobaric_divided_difference_on_basis(self, key):
            r"""
            On basis action for the hat isobaric divided difference on the
            Grothendieck basis.

            TESTS::

                sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
                sage: A.<x,y> = MultivariatePolynomialAlgebra(QQ)
                sage: DG = A.double_grothendieck_basis()
                sage: DG3 = DG.finite_rank_basis(3)
                sage: wrapp = DG._divided_difference_wrapper(DG3,1)
                sage: pol = DG[3,2,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                y[0]*GG[2, 2, 3] + (-y[0])*GG[3, 2, 3]
                sage: pol = DG[2,3,3]
                sage: key = list(pol)[0][0]
                sage: wrapp.hat_isobaric_divided_difference_on_basis(key)
                0

            Test consistency..::

                sage: pol = DG[3,2,3]
                sage: a = pol.hat_isobaric_divided_difference(1).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.hat_isobaric_divided_difference(2).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(2)
                sage: a == b
                True
                sage: pol = DG[2,3,3]
                sage: a = pol.hat_isobaric_divided_difference(1).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(1)
                sage: a == b
                True
                sage: a = pol.hat_isobaric_divided_difference(2).expand()
                sage: b = pol.expand().hat_isobaric_divided_difference(2)
                sage: a == b
                True
            """
            return self.isobaric_divided_difference_on_basis(key) - self._module(key)

