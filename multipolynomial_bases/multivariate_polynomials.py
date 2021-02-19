# -*- coding: utf-8 -*-
r"""
Multivariate Polynomial Algebra

This modules implements the polynomial ring seen as an algebra with multibases.
Especially, it implements bases such as the Schubert, Grothendieck, and
Key polynomials, and any basis based on a divided difference type operation.

In the monomial basis, a multivariate polynomial is seen as a linear combination
of vectors. Where each vector represents the exponents of the given monomial.

The number of variables is not set: the algebra can be understood as the projective limit
of all polynomial rings with a finite number of variables.

EXAMPLES::

    sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
    sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
    sage: A
    The Multivariate polynomial algebra on x over Rational Field
    sage: A.an_element()
    x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
    sage: x
    The Multivariate polynomial algebra on x over Rational Field on the monomial basis
    sage: x[1,1,2] + x[3,2,4]
    x[1, 1, 2] + x[3, 2, 4]
    sage: x[1,1,2] * x[3,2,4]
    x[4, 3, 6]

Here is how to access a single variable::

    sage: x1 = A.var(1)
    sage: x1
    x[1]
    sage: x2 = A.var(2)
    sage: x2
    x[0, 1]
    sage: x1 * x2
    x[1, 1]

Get back a symbolic expression::

    sage: pol = A.an_element(); pol
    x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
    sage: pol.to_expr()
    x1*x2^2*x3^3 + 2*x1 + 3*x2 + 1

You can apply many different actions on the polynomial::

    sage: pol = A.an_element(); pol
    x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
    sage: pol.divided_difference(2)
    3*x[0, 0, 0] - x[1, 2, 2]
    sage: pol.isobaric_divided_difference(2)
    x[0, 0, 0] + 3*x[0, 0, 1] + 3*x[0, 1, 0] + 2*x[1, 0, 0]
    sage: pol.hat_isobaric_divided_difference(2)
    3*x[0, 0, 1] - x[1, 2, 3]
    sage: pol.si(1)
    x[0, 0, 0] + 2*x[0, 1, 0] + 3*x[1, 0, 0] + x[2, 1, 3]

The main purpose of this module is to implement different bases of the
polynomial algebra based on these actions::

    sage: Schub = A.schubert_basis()
    sage: K = A.demazure_basis()
    sage: Khat = A.demazure_hat_basis()
    sage: Groth = A.grothendieck_positive_basis()
    sage: Schub(pol)
    Y[0, 0, 0] + 3*Y[0, 1, 0] - Y[1, 0, 0] + Y[1, 2, 3] - Y[1, 3, 2] - Y[2, 1, 3] + Y[2, 3, 1] + Y[3, 1, 2] - Y[3, 2, 1] + Y[4, 1, 1]
    sage: K(pol)
    K[0, 0, 0] + 3*K[0, 1, 0] - K[1, 0, 0] + K[1, 2, 3] - K[1, 3, 2] - K[2, 1, 3] + K[2, 3, 1] + K[3, 1, 2] - K[3, 2, 1]
    sage: Khat(pol)
    ^K[0, 0, 0] + 3*^K[0, 1, 0] + 2*^K[1, 0, 0] + ^K[1, 2, 3]
    sage: Groth(pol)
    G[0, 0, 0] + 3*G[0, 1, 0] - G[1, 0, 0] + 3*G[1, 1, 0] + G[1, 2, 3] - G[1, 3, 2] + G[1, 3, 3] - G[2, 1, 3] + G[2, 2, 3] + G[2, 3, 1] - 2*G[2, 3, 2] + G[2, 3, 3] + G[3, 1, 2] - G[3, 1, 3] - G[3, 2, 1] + G[3, 2, 2] - G[3, 3, 2] + G[3, 3, 3] + G[4, 1, 1] - G[4, 1, 3]
"""
from __future__ import absolute_import
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.graded_algebras import GradedAlgebras
from sage.categories.category import Category
from sage.categories.graded_modules import GradedModules
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.realizations import Realizations
from sage.categories.rings import Rings
from sage.categories.category_types import Category_over_base_ring
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.all import var

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

def MultivariatePolynomialAlgebra(base_ring, names = None, set_of_variables = 1):
    r"""
    Return the multivariate polynomial algebra.

    This is a polynomial ring in one or two inifinite sets of variables,
    interpreted as a multibases algebra.

    INPUT:

        - ``base_ring`` the base ring of the algebra
        - ``names`` a list of maximal size 2 of names for the set of variables
        - ``set_of_variables`` an integer default:1, if names is not set, then
          ``set_of_variables`` is used to know on how many set of variables
          the algebra is defined (max = 2)

    OUTPUT:

    The multivariate polynomial algebra.

    EXAMPLES::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
        sage: A
        The Multivariate polynomial algebra on x over Rational Field

    The generator ``x`` is not a single variable, it represents a infinite
    set of variables. More precisely, it is the monomial basis of the algebra.

    ::

        sage: x
        The Multivariate polynomial algebra on x over Rational Field on the monomial basis
        sage: x[1,1,2] + x[2,3,4]
        x[1, 1, 2] + x[2, 3, 4]
        sage: x == A.monomial_basis()
        True

    Here is how to access a single variable::

        sage: x1 = A.var(1)
        sage: x1
        x[1]
        sage: x2 = A.var(2)
        sage: x2
        x[0, 1]
        sage: x1 * x2
        x[1, 1]

    TESTS::

        sage: A = MultivariatePolynomialAlgebra(QQ); A
        The Multivariate polynomial algebra on x over Rational Field
        sage: A = MultivariatePolynomialAlgebra(QQ, names = ["x"]); A
        The Multivariate polynomial algebra on x over Rational Field
        sage: A = MultivariatePolynomialAlgebra(QQ, names = ["x", "y"]); A
        The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field
        sage: A = MultivariatePolynomialAlgebra(QQ, set_of_variables = 2); A
        The Multivariate polynomial algebra on x over The Multivariate polynomial algebra on y over Rational Field

    """
    assert(base_ring in Rings())
    if names is None:
        main_repr_var = 'x'
        main_repr_var2 = 'y'
    elif len(names) > 2 or set_of_variables > 2:
        raise ValueError("The multivariate polynomial algebra only allows two sets of variables")
    else:
        main_repr_var = names[0]
        if len(names) > 1:
            set_of_variables = 2
            main_repr_var2 = names[1]
    if not set_of_variables == 2:
        return MultivariatePolynomialAlgebra_generic(base_ring, main_repr_var)
    else:
        from .double_multivariate_polynomials import DoubleMultivariatePolynomialAlgebra_generic
        return DoubleMultivariatePolynomialAlgebra_generic(base_ring, main_repr_var, main_repr_var2)



class MultivariatePolynomialAlgebra_generic(UniqueRepresentation, Parent):
    r"""
    A class implementing the multivariate polynomial ring as multibases
    algebra.

    INPUT:

    - ``R``: the base ring of the algebra
    - ``main_repr_var``, the letter corresponding to the set of variables,
      it is used to represent several bases, default is ``x``
    - ``always_show_main_var``, if True ``main_repr_var`` will be displayed
      on elements of every basis, even the ones that don't use it directly
      (Schubert basis, Demazure basis, ...), false by default, used on
      ``DoubleMultivariatePolynomialAlgebra`` to differentiate the two sets of
      variables


    OUTPUT:

    - The Multivariate polynomial algebra on ``main_repr_var`` over ``R``

    EXAMPLES::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A.<x> = MultivariatePolynomialAlgebra(QQ); A
        The Multivariate polynomial algebra on x over Rational Field

    The monomial basis is given as the algebra generator::

        sage: x
        The Multivariate polynomial algebra on x over Rational Field on the monomial basis
        sage: A.monomial_basis() == x
        True

    You can use it to create polynomials.

    ::

        sage: pol = x[1,1,2] + x[3,4]; pol
        x[1, 1, 2] + x[3, 4, 0]
        sage: pol * x[2,3]
        x[3, 4, 2] + x[5, 7, 0]

    You can also access a single variable.

    ::

        sage: A.var(1)
        x[1]
        sage: A.var(2)
        x[0, 1]
        sage: A.var(3)
        x[0, 0, 1]

    The coercion between elements with a different number of variables is done
    automatically.

    ::

        sage: pol1 = x[1,1,2]
        sage: pol1.nb_variables()
        3
        sage: pol2 = x[2]
        sage: pol2.nb_variables()
        1
        sage: (pol1 * pol2).nb_variables()
        3

    TESTS::

        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: A = MultivariatePolynomialAlgebra(QQ)
        sage: TestSuite(A).run()
    """
    def __init__(self, R, main_repr_var, always_show_main_var = False, extra_category = None, extra_bases_category = None):
        r"""
        TESTS::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: A
            The Multivariate polynomial algebra on x over Rational Field
            sage: x
            The Multivariate polynomial algebra on x over Rational Field on the monomial basis
        """
        category = [GradedAlgebras(R).WithRealizations(), CommutativeAlgebras(R)]
        if extra_category is not None:
            category.add(extra_category)
        Parent.__init__(
            self,
            base = R,
            category = category
        )
        self._finite_rings = set([])
        self._main_repr_var = main_repr_var
        self._show_main_var = always_show_main_var
        self._extra_bases_category = extra_bases_category

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
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: x
            The Multivariate polynomial algebra on x over Rational Field on the monomial basis
            sage: x == A.gens()[0]
            True
            sage: x[1,2,3]
            x[1, 2, 3]
        """
        return [self.monomial_basis()]

    def facade_for(self):
        r"""
        Return all the parents ``self`` is a facade for

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(CC)
            sage: A.facade_for()
            []

        The facade_for parents are added dynamically.

        """
        l = [r for r in self.realizations()]
        l.extend(self._finite_rings)
        for r in self.realizations():
            l.extend(r.facade_for())
        return l

    def _repr_(self):
        r"""
        Return the string representation of ``self``

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: MultivariatePolynomialAlgebra(QQ)
            The Multivariate polynomial algebra on x over Rational Field
        """
        return "The Multivariate polynomial algebra on %s over %s"%(
            self._main_repr_var,
            self.base_ring()
        )


    def a_realization(self):
        r"""
        Return a default realization of ``self`` : the monomial basis

        EXAMPLES ::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: A.a_realization()
            The Multivariate polynomial algebra on x over Rational Field on the monomial basis

        """
        return self.monomial_basis()

    def an_element(self):
        r"""
        Return an element of ``self``. By default, this element lies in the
        monomial basis.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: A.an_element()
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
        """
        return self.monomial_basis().an_element()

    def _element_constructor_(self, element):
        r"""
        Construct an element of ``self``.

        As ``self`` is an abstract algebra, this method will
        just check if ``element`` belongs to ``self``

        INPUT:

        -``element`` the element to be contructed from

        OUTPUT:

        The element itself if it belongs to ``self``, if not,
        it raises a TypeError Exception

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: p = A.an_element(); p
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
            sage: A._element_constructor_(p)
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
            sage: A._element_constructor_(1)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make '1' an element of 'The Multivariate polynomial algebra on x over Rational Field'
        """
        if any (parent.is_parent_of(element) for parent in self.facade_for()):
            return element
        raise TypeError("do not know how to make '%s' an element of '%s'"%(element, self))

    def var(self, i, nb_variables = 0):
        r"""
        Return the i_th variable as a monomial base element

        INPUT:

        - ``i``: the index of the variable to return
        - ``nb_variables``: the number of variables of the result,
          default is ``i``, if ``nb_variables`` is lower than ``i`` it is
          ignored and changed to ``i``

        OUTPUT:

        - the ith variable as a monomial element

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ);
            sage: A.var(1)
            x[1]
            sage: A.var(3)
            x[0, 0, 1]
            sage: A.var(1,3)
            x[1, 0, 0]
            sage: A.var(4,3)
            x[0, 0, 0, 1]

        """
        if(nb_variables ==0 or nb_variables < i ): nb_variables =i
        vect = [0 for j in range(nb_variables)]
        vect[i-1] = 1
        return self.monomial_basis()( vect )

    # TODO : rewrite !
    #def from_expr(self, expr, alphabet = None, second_alphabet = None):
        #r"""
        #Construct a polynomial from a symbolic expression.

        #INPUT:
            #- ``expr`` a symbolic expression, it must be a polynomial
              #in all variables of ``variables``
            #- ``alphabet`` (optional), a list of symbolic variables.
              #If not set, it takes ``expr.variables()``. The variables are matched
              #to the vector key of the monomials by the order of the list.
            #- ``second_alphabet`` (optional) a list of symbolic variables
              #when working on a polynomial on 2 sets of variables

        #EXAMPLES::

            #sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            #sage: A = MultivariatePolynomialAlgebra(QQ)
            #sage: var('x1,x2,x3')
            #(x1, x2, x3)
            #sage: expr = 3*x3 + x2^2 - x1*x3
            #sage: A.from_expr(expr)
            #-x[1, 0, 1] + x[0, 2, 0] + 3*x[0, 0, 1]
            #sage: var('t1,t2')
            #(t1, t2)
            #sage: K.<t1,t2> = QQ[]
            #sage: K = K.fraction_field()
            #sage: A = MultivariatePolynomialAlgebra(K)
            #sage: expr = t1*t2*x2 - x3^4*x1*4*t2^2
            #sage: A.from_expr(expr,[x1,x2,x3])
            #(-4*t2^2)*x[1, 0, 4] + t1*t2*x[0, 1, 0]

        #TODO Works with polynomials in two sets of variables

            #sage: D = DoubleMultivariatePolynomialAlgebra(QQ)
            #sage: var('x1,x2,x3,y1,y2,y3')
            #(x1, x2, x3, y1, y2, y3)
            #sage: expr = x1*y1 +(y2*y3^2 - y1)*x3*x1^4
            #sage: D.from_expr(expr,[x1,x2,x3],[y1,y2,y3])
            #(y[1,0,0])*x[1, 0, 0] + (-y[1,0,0]+y[0,1,2])*x[4, 0, 1]

        #"""
        #return self.monomial_basis().from_expr(expr, alphabet, second_alphabet)

    def algebra_finite_nb_variables(self, nb_variables, basis_repr = None):
        r"""
        Return the realization of ``self`` in a given number of variables.

        INPUT:

        - ``nb_variables``: the number of variables
        - ``basis_repr``, the representation letter for the elements of the base,
          by default, it is the main representation for the set of variable :
          ``self._main_repr_var``

        OUTPUT:

        The multivariate polynomial algebra in ``nb_variables`` variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: A3 = A.algebra_finite_nb_variables(3); A3
            The Multivariate polynomial algebra on x over Rational Field with 3 variables

        The finite number of variables ring contains method to obtain the algebra
        bases on a finite number of variables::

            sage: m3 = A3.monomial_basis(); m3
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
            sage: ma3 = A3.monomial_basis_with_type("A"); ma3
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the Monomial basis of type A

        Coercions between rings with a different number of variables are created
        dynamically::

            sage: x = A.monomial_basis()
            sage: pol1 = x[1,2,3]; pol1
            x[1, 2, 3]
            sage: pol1.parent()
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
            sage: pol2 = x[1,1]; pol2
            x[1, 1]
            sage: pol2.parent()
            The Multivariate polynomial algebra on x over Rational Field with 2 variables on the monomial basis
            sage: pol1 + pol2
            x[1, 1, 0] + x[1, 2, 3]
            sage: (pol1 + pol2).parent()
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        F = FiniteRankMultivariatePolynomialAlgebra(self,nb_variables, basis_repr, extra_bases_category = self._extra_bases_category)
        return F

    def _register_finite_ring(self,F):
        r"""
        Add ``F`` as one of ``self`` finite variables realizations. It is called
        by the init function of ``F``

        INPUT:

        - ``F`` -- a realization of ``self`` in a given number of variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: var('t')
            t
            sage: A = MultivariatePolynomialAlgebra(QQ[t])
            sage: A._finite_rings
            set()
            sage: F3 = A.algebra_finite_nb_variables(3)
            sage: A._finite_rings
            {The Multivariate polynomial algebra on x over Univariate Polynomial Ring in t over Rational Field with 3 variables}
        """
        if(not(F in self._finite_rings)):
            for ring in self._finite_rings:
                if( ring.nb_variables() != F.nb_variables() ): self._create_morphism(ring, F)
            self._finite_rings.add(F)

    def monomial_basis(self, basis_repr = None):
        r"""
        Return the monomial basis of ``self``.

        INPUT:

        - ``basis_repr``, the representation letter for the elements of the base,
          by default, it is the main representation for the set of variable :
          ``self._main_repr_var``

        OUTPUT:

        The monomial basis of the multivariate polynomial algebra.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: x = A.monomial_basis();x
            The Multivariate polynomial algebra on x over Rational Field on the monomial basis
            sage: x[1,2,3]
            x[1, 2, 3]
            sage: x( [2,2] )
            x[2, 2]
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        from .monomial import MonomialBasis
        return MonomialBasis(self, basis_repr)

    def monomial_basis_with_type(self, group_type, basis_repr = None):
        r"""
        Return a typed monomial basis.

        Monomials are indexed by a root system lattice. They embed a group
        type and all divided difference are done within this group type.

        INPUT:

        - ``group_type``: -- the letter that represents type of the weyl group,
          can be either ``"A"``, ``"B"``, ``"C"``, or ``"D"``
        - ``basis_repr`` -- (optional) the representation letter for the elements of the base,
          by default, it is using both ``self._main_repr_var`` and ``group_type``

        OUTPUT:

        The monomial basis with type ``group_type``.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: xA = A.monomial_basis_with_type("A"); xA
            The Multivariate polynomial algebra on x over Rational Field on the Ambient space basis of type A
            sage: xB = A.monomial_basis_with_type("B")
            sage: xB
            The Multivariate polynomial algebra on x over Rational Field on the Ambient space basis of type B
            sage: xA == x
            False
            sage: xA == xB
            False
            sage: xA.group_type()
            'A'


        Default coercions are created between the typed and untyped basis::

            sage: xA( x[1,2,3])
            xA[1, 2, 3]
            sage: x( xA[2,4] )
            x[2, 4]
            sage: xB( x[1,2,3])
            xB[1, 2, 3]
            sage: x( xB[2,4])
            x[2, 4]

        """
        if(basis_repr is None): basis_repr = self._main_repr_var + group_type
        from .ambient_space_basis import PolynomialRingWithBasisFromAmbientSpace
        return PolynomialRingWithBasisFromAmbientSpace(self,group_type,"Ambient space basis of type " + group_type, basis_repr)

    def from_morphism_basis(self, neutral_nb_variables, morphism_to_basis, get_basis_keys, get_morphism_on_basis, basis_name, basis_repr, variables_auto_coerce =False, **keywords):
        r"""
        Create a basis defined by its morphism to another basis

        INPUT:

        - ``neutral_nb_variables`` -- the default number of variables to get the
            one element
        - ``morphism_to_basis`` -- the basis of the polynomial algebra on
            which the morphism will be defined
        - ``get_basis_keys`` -- a function with :
                input:
                    - ``nb_variables`` -- the number of variables
                output:
                    - the set of indexes that will be used to index elements
                      of the basis on the given
                      number of variables
        - ``get_morphism_on_basis`` -- a function with :
                input:
                    -``nb_variables``, the number of variables
                output:
                    - the function that will be used to create the module morphims on
                      basis on the given number of variables
        - ``basis_name`` -- the name of the basis (used in repr)
        - ``basis_repr``-- the basis representation for elements (exemple ``x``)
        - ``variables_auto_coerce`` -- if set to ``True``, a coercion will be
          created between elements of the basis indexed by vectors of size
          n to basis on m>n variables by extending the vectors with zeros
          (example: x[2,2,1] -> x[2,2,1,0,0]. Default is ``False``.
        - ``**keywords`` -- the keywords sent to the ``CombinatorialFreeModule``
          morphism.

        OUTPUT:

        - the basis of which elements are indexed by the sets return
          by ``get_basis_keys`` and can be coerced on
          the ``morphims_to_basis`` basis

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: m = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space(QQ)
            sage: def get_morphism_on_basis(n): return lambda key: m( [key[i] for i in range(n)])
            sage: MyBasis = A.from_morphism_basis(1,m,get_basis_keys,get_morphism_on_basis,"My Basis", "x"); MyBasis
            The Multivariate polynomial algebra on x over Rational Field on the My Basis
            sage: MyBasis.an_element()
            x(2, 2, 3)
            sage: m( MyBasis.an_element() )
            x[2, 2, 3]

        We have recreated the typed basis.
        """
        from .basis import PolynomialRingWithBasisFromMorphism
        return PolynomialRingWithBasisFromMorphism(self, neutral_nb_variables, morphism_to_basis, basis_name, basis_repr, get_basis_keys, get_morphism_on_basis,variables_auto_coerce, **keywords)

    def linear_basis_on_vectors(self, group_type, basis_name, basis_repr, on_basis_method, extra_parameters = (), **keywords):
        r"""
        Create a linear basis on objects inedexed by vectors based on an operation
        to convert each object (through its vector) into a typed polynomial.

        INPUT:

        - ``group_type`` -- the letter that represents the type of the weyl group that
           will be used for the ambient space basis
        - ``basis_name`` -- the name of the basis (used in repr)
        - ``basis_repr``-- the basis representation for elements
        - ``on_basis_method`` -- a method that takes a vector (python list) and returns
          the converted polynomial associated with it
          The ``on_basis_method`` should have the following signature :
          Input:

          - ``v`` a python list representing the vector
          - ``basis`` the ambient space basis used to make the conversion
          - ``call_back`` a call_back method to use the conversion recursively
          - ``**keywords`` extra parameters that could be used for convertion

          Output:

          - a polynomial expanded into the sent ``basis`` and corresponding to
            the objected indexed by the sent vector ``v``

        - ``extra_parameters``: (default : empty) a tuple containing the extra parameters
          to be sent to the ``on_basis_method`` as tuples ``(key,val)``
        - ``**keyword`` : parameters used to create the morphism to the ambient space basis,
          sent to ``CombinatorialFreeModule.module_morphism``. By default, ``triangular`` is
          set to ``upper`` : change it explicitly to ``None`` if you're basis is not.
          A default ``cmp`` method is also used to order keys : it compares degrees first and
          then orders vectors by inversed lexical order. This comparasion method is the classic
          one, used by Schubert, Demazure and Grothendieck polynomials.

        OUTPUT :

            - a basis named ``basis_name`` and defined by its conversion to the ambient space basis of type ``group_type``

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: def schubert_on_basis(v, basis, call_back):
            ....:    for i in range(len(v)-1):
            ....:        if(v[i]<v[i+1]):
            ....:            v[i], v[i+1] = v[i+1] + 1, v[i]
            ....:            return call_back(v).divided_difference(i+1)
            ....:    return basis(v)
            sage: myBasis = A.linear_basis_on_vectors("A","MySchub","Y",schubert_on_basis)
            sage: pol = myBasis[2,1,3];pol
            Y[2, 1, 3]
            sage: pol.expand()
            xA[2, 1, 3] + xA[2, 2, 2] + xA[2, 3, 1] + xA[3, 1, 2] + xA[3, 2, 1] + xA[4, 1, 1]
            sage: myBasis(A.an_element())
            Y[0, 0, 0] + 3*Y[0, 1, 0] - Y[1, 0, 0] + Y[1, 2, 3] - Y[1, 3, 2] - Y[2, 1, 3] + Y[2, 3, 1] + Y[3, 1, 2] - Y[3, 2, 1] + Y[4, 1, 1]

        We have recreated the Schubert basis. Let's see an example with parameters::

            sage: def t_inverse(v, basis, call_back, t1=1, t2=1): return t1 * (t2*basis(v))**(-1)
            sage: tInverse = A.linear_basis_on_vectors("A","tInverse","T",t_inverse, (("t1",2), ("t2",4)), triangular = None)
            sage: pol = tInverse[1,1,2]; pol
            T[1, 1, 2]
            sage: pol.expand()
            1/2*xA[-1, -1, -2]

        """
        from .linear_basis_on_vectors import LinearBasisOnVectors
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return LinearBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, on_basis_method, extra_parameters, **keywords)

    def schubert_basis(self, basis_name = None, basis_repr = "Y"):
        r"""
        Creates the simple Schubert basis where schubert polynomials are indexed
        by vectors.

        Here is the definition we use. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define

        $Y_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.

        Otherwise, we have for $ v_i > v_{i+1}$

        $Y_{\cdots v_{i+1} v_i-1 \cdots} = Y_v \partial_i$ where $\partial_i$ is the ith divided difference.

        The vectors indexing the Schubert polynomials can as well been seen as
        lehmer codes.

        INPUT:

        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (defaul: ``Y``) the basis representation for elements

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the Schubert basis
          of type ``group_type`` index by vectors where ``R`` is the algebra base
          ring defined in the abstract algebra.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: Schub = A.schubert_basis(); Schub
            The Multivariate polynomial algebra on x over Rational Field on the Schubert basis of type A
            sage: Schub.an_element()
            Y[0, 0, 0] + 3*Y[0, 1, 0] + 2*Y[1, 0, 0] + Y[2, 2, 3]
            sage: Schub[1,2,3]
            Y[1, 2, 3]

        Let us see the coercions::

            sage: y = Schub.an_element(); y
            Y[0, 0, 0] + 3*Y[0, 1, 0] + 2*Y[1, 0, 0] + Y[2, 2, 3]
            sage: y.expand()
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 5*xA[1, 0, 0] + xA[2, 2, 3] + xA[2, 3, 2] + xA[3, 2, 2]
            sage: xA = A.monomial_basis_with_type("A")
            sage: xA(y)
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 5*xA[1, 0, 0] + xA[2, 2, 3] + xA[2, 3, 2] + xA[3, 2, 2]
            sage: x = A.monomial_basis()
            sage: x (y )
            x[0, 0, 0] + 3*x[0, 1, 0] + 5*x[1, 0, 0] + x[2, 2, 3] + x[2, 3, 2] + x[3, 2, 2]
            sage: xA.an_element(); Schub( xA.an_element())
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 2*xA[1, 0, 0] + xA[2, 2, 3]
            Y[0, 0, 0] + 3*Y[0, 1, 0] - Y[1, 0, 0] + Y[2, 2, 3] - Y[2, 3, 2]
            sage: x.an_element(); Schub(x.an_element())
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
            Y[0, 0, 0] + 3*Y[0, 1, 0] - Y[1, 0, 0] + Y[1, 2, 3] - Y[1, 3, 2] - Y[2, 1, 3] + Y[2, 3, 1] + Y[3, 1, 2] - Y[3, 2, 1] + Y[4, 1, 1]

        Let us see some operations::

            sage: Schub[1,2] + Schub[3,0,0]
             Y[1, 2, 0] + Y[3, 0, 0]
            sage: Schub.an_element() * Schub.an_element()
            Y[0, 0, 0] + 6*Y[0, 1, 0] + 9*Y[0, 2, 0] + 4*Y[1, 0, 0] + 21*Y[1, 1, 0] + 16*Y[2, 0, 0] + 2*Y[2, 2, 3] + 6*Y[2, 3, 3] + 6*Y[2, 4, 2] + 4*Y[3, 2, 3] + Y[4, 4, 6] + Y[4, 5, 5]
            sage: Schub[1,2] * Schub[3,0,0]
            Y[4, 2, 0] + Y[5, 1, 0]
        """
        from .linear_basis_on_vectors import SchubertBasisOnVectors
        if(basis_name is None):
            basis_name = "Schubert basis of type A"
        if(self._show_main_var): basis_repr+= self._main_repr_var
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return SchubertBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr)

    def demazure_basis(self, group_type ="A", basis_name = None, basis_repr = "K"):
        r"""
        Creates the Demazure basis where demazure / key polynomials are indexed
        by vectors.

        Here is the definition we use for type A. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$,
        we define

        $K_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.

        Otherwise, we have for $ v_i > v_{i+1}$

        $K_{\cdots v_{i+1} v_i \cdots} = K_v \pi_i$ where $\pi_i$ is the ith isobar divided difference.

        The vectors indexing the key polynomials can as well been seen as lehmer codes.

        INPUT:

        - ``group_type``: (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``K``) the basis representation for elements

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the Demazure basis
          of type ``group_type`` indexed by vectors where ``R`` is the algebra
          base ring defined in the abstract algebra.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: Dem = A.demazure_basis("A"); Dem
            The Multivariate polynomial algebra on x over Rational Field on the Demazure basis of type A
            sage: Dem.an_element()
            K[0, 0, 0] + 3*K[0, 1, 0] + 2*K[1, 0, 0] + K[2, 2, 3]
            sage: Dem[1,2,3]
            K[1, 2, 3]

        Let us see some coercions::

            sage: k = Dem.an_element(); k
            K[0, 0, 0] + 3*K[0, 1, 0] + 2*K[1, 0, 0] + K[2, 2, 3]
            sage: k.expand()
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 5*xA[1, 0, 0] + xA[2, 2, 3] + xA[2, 3, 2] + xA[3, 2, 2]
            sage: xA = A.monomial_basis_with_type("A")
            sage: xA(k)
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 5*xA[1, 0, 0] + xA[2, 2, 3] + xA[2, 3, 2] + xA[3, 2, 2]
            sage: x = A.monomial_basis()
            sage: x(k)
            x[0, 0, 0] + 3*x[0, 1, 0] + 5*x[1, 0, 0] + x[2, 2, 3] + x[2, 3, 2] + x[3, 2, 2]
            sage: xA.an_element(); Dem(xA.an_element())
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 2*xA[1, 0, 0] + xA[2, 2, 3]
            K[0, 0, 0] + 3*K[0, 1, 0] - K[1, 0, 0] + K[2, 2, 3] - K[2, 3, 2]
            sage: x.an_element(); Dem(x.an_element())
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
            K[0, 0, 0] + 3*K[0, 1, 0] - K[1, 0, 0] + K[1, 2, 3] - K[1, 3, 2] - K[2, 1, 3] + K[2, 3, 1] + K[3, 1, 2] - K[3, 2, 1]

        Let us see some operations::

            sage: Dem[1,2] + Dem[3,0,0]
            K[1, 2, 0] + K[3, 0, 0]
            sage: Dem.an_element() * Dem.an_element()
            K[0, 0, 0] + 6*K[0, 1, 0] + 9*K[0, 2, 0] + 4*K[1, 0, 0] + 21*K[1, 1, 0] + 16*K[2, 0, 0] + 2*K[2, 2, 3] + 6*K[2, 3, 3] + 6*K[2, 4, 2] + 4*K[3, 2, 3] + 4*K[4, 2, 2] + K[4, 4, 6] + K[4, 5, 5]
            sage: Dem[1,2] * Dem[3,0,0]
            K[4, 2, 0] + K[5, 1, 0]

        We can also have type B, C or D key polynomials::

            sage: DemB = A.demazure_basis("B"); DemB
            The Multivariate polynomial algebra on x over Rational Field on the Demazure basis of type B
            sage: pol = DemB[2,1,-2]; pol
            K[2, 1, -2]
            sage: pol.expand()
            xB[2, 1, -2] + xB[2, 1, -1] + xB[2, 1, 0] + xB[2, 1, 1] + xB[2, 1, 2] + xB[2, 2, -1] + xB[2, 2, 0] + xB[2, 2, 1]


        """
        from .linear_basis_on_vectors import DemazureBasisOnVectors
        if(basis_name is None):
            basis_name = "Demazure basis of type " +group_type
        if(self._show_main_var): basis_repr+= self._main_repr_var
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return DemazureBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, "isobaric_divided_difference")

    def demazure_hat_basis(self, group_type ="A", basis_name = None, basis_repr = "^K"):
        r"""
        Creates the Demazure hat basis where demazure polynomials are indexed
        by vectors.

        Here is the definition we use for type A. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define

        $K_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.

        Otherwise, we have for $ v_i > v_{i+1}$

        $K_{\cdots v_{i+1} v_i \cdots} = K_v \pi_i$ where $\pi_i$ is the ith isobar hat divided difference.

        The vectors indexing the key polynomials can as well been seen as lehmer codes.

        INPUT:

        - ``group_type``: (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``^K``) the basis representation for elements

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the Demazure hat basis
          of type ``group_type`` indexed by vectors where ``R`` is the algebra
          base ring defined in the abstract algebra.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ);
            sage: Demh = A.demazure_hat_basis("A"); Demh
            The Multivariate polynomial algebra on x over Rational Field on the Demazure hat basis of type A
            sage: Demh.an_element()
            ^K[0, 0, 0] + 3*^K[0, 1, 0] + 2*^K[1, 0, 0] + ^K[2, 2, 3]

        Let us see some coercions::

            sage: kh = Demh[1,2,4]; kh
            ^K[1, 2, 4]
            sage: kh.expand()
            xA[1, 2, 4] + xA[1, 3, 3] + xA[2, 2, 3]
            sage: xA = A.monomial_basis_with_type("A")
            sage: xA(kh)
            xA[1, 2, 4] + xA[1, 3, 3] + xA[2, 2, 3]
            sage: x = A.monomial_basis()
            sage: x(kh)
            x[1, 2, 4] + x[1, 3, 3] + x[2, 2, 3]
            sage: Demh(xA[1,2,4])
            ^K[1, 2, 4] - ^K[1, 3, 3] + ^K[2, 3, 2]
            sage: Demh(x[1,2,4])
            ^K[1, 2, 4] - ^K[1, 3, 3] + ^K[2, 3, 2]

        Let us see some operations::

            sage: Demh[1,2] + Demh[3,0,0]
            ^K[1, 2, 0] + ^K[3, 0, 0]
            sage: Demh.an_element() * Demh.an_element()
            ^K[0, 0, 0] + 6*^K[0, 1, 0] + 9*^K[0, 2, 0] + 4*^K[1, 0, 0] + 3*^K[1, 1, 0] + 4*^K[2, 0, 0] + 2*^K[2, 2, 3] + 6*^K[2, 3, 3] + 4*^K[3, 2, 3] + ^K[4, 4, 6] - ^K[4, 5, 5] - ^K[5, 4, 5]
            sage: Demh[1,2] * Demh[3,0,0]
            ^K[4, 2, 0]

        We can also have type B, C or D hat key polynomials::

            sage: DemhB = A.demazure_hat_basis("B"); DemhB
            The Multivariate polynomial algebra on x over Rational Field on the Demazure hat basis of type B
            sage: pol = DemhB[2,1,-2]; pol
            ^K[2, 1, -2]
            sage: pol.expand()
            xB[2, 1, -2] + xB[2, 1, -1] + xB[2, 1, 0] + xB[2, 1, 1]

        """
        from .linear_basis_on_vectors import DemazureBasisOnVectors
        if(basis_name is None):
            basis_name = "Demazure hat basis of type " +group_type
        if(self._show_main_var): basis_repr+= self._main_repr_var
        monomial_basis_with_type = self.monomial_basis_with_type(group_type)
        return DemazureBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, "hat_isobaric_divided_difference")


    def grothendieck_negative_basis(self, basis_name = None, basis_repr = "G"):
        r"""
        Creates the simple Grothendieck basis where Grothendieck polynomials are indexed
        by vectors. The Negative stands for that we use a definition of the basis where
        the variable have negative exposants.

        Here is the definition we use. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define

        $G_v = \prod_{1 \cdots n} (1 - \frac{1}{x_i})^{v_i}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.

        Otherwise, we have for $ v_i > v_{i+1}$

        $G_{\cdots v_{i+1} v_i-1 \cdots} = G_v \pi_i$ where $\pi_i$ is the ith isobar divided difference.

        The vectors indexing the Grothendieck polynomials can as well been seen as lehmer codes.



        INPUT:

        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (defaul: ``G``) the basis representation for elements

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the
          Grothendieck basis of type ``group_type`` with negative exposants
          indexd by vectors where ``R`` is the algebra base ring defined in
          the abstract algebra.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: Groth = A.grothendieck_negative_basis(); Groth
            The Multivariate polynomial algebra on x over Rational Field on the Grothendieck basis of type A with negative exposants
            sage: Groth.an_element()
            G[0, 0, 0] + 3*G[0, 1, 0] + 2*G[1, 0, 0] + G[2, 2, 3]
            sage: Groth[1,2,3]
            G[1, 2, 3]

        We can convert a Grothendieck polynomial into the Monomial or
        Ambient Space basis but not the other way around as Grothendieck
        polynomials are with negative exposants. Note that conversion
        from monomials with negative exposants into Grothendieck polynomials
        is NOT implemented ::

            sage: g = Groth[0,1]; g
            G[0, 1]
            sage: g.expand()
            -xA[-1, -1] + xA[0, 0]
            sage: xA = A.monomial_basis_with_type("A")
            sage: xA(g)
            -xA[-1, -1] + xA[0, 0]
            sage: x = A.monomial_basis()
            sage: x(g)
            -x[-1, -1] + x[0, 0]
            sage: pol = x[0,0] - x[-1,-1]
            sage: Groth( pol)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= -x[-1, -1] + x[0, 0]) an element of self (=The Multivariate polynomial algebra on x over Rational Field with 2 variables on the Grothendieck basis of type A with negative exposants)

        We can add Grothendieck polynomials but not multiply them as this
        would use conversion from monomials into Grothendieck polynomials ::

            sage: Groth[1,2] + Groth[3,0,0]
            G[1, 2, 0] + G[3, 0, 0]
            sage: Groth[1,2] * Groth[3,0,0]
            Traceback (most recent call last):
            ...
            NotImplementedError: The product is not implemented for this basis
        """
        from .linear_basis_on_vectors import GrothendieckNegativeBasisOnVectors
        if(basis_name is None):
            basis_name = "Grothendieck basis of type A with negative exposants"
        if(self._show_main_var): basis_repr+= self._main_repr_var
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return GrothendieckNegativeBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr)

    def grothendieck_positive_basis(self, basis_name = None, basis_repr = "G"):
        r"""
        Creates the simple Grothendieck basis where Grothendieck polynomials
        are indexed by vectors. The positive stands for that we use a definition
        of the basis where variables have positive exposants.

        It corresponds to the basis given by ``grothendieck_negative_basis``
        by a change of variables :

        $x_i = 1 - \fract{1}{x_i}$

        Here is the definition we use. For $v = (v_1, \cdots, v_n) \in \mathbb{N}^n$, we define

        $G_v = x_1^{v_1}x_2^{v_2}\cdotsx_n^{v_n}$ if $v$ is a partition, i.e, if $v_1 \geq v_2 \geq \cdots \geq v_n$.

        Otherwise, we have for $ v_i > v_{i+1}$

        $G_{\cdots v_{i+1} v_i-1 \cdots} = \left( G_v \frac{1 - x_{i+1}}{x_i} \right) \pi_i$ where $\pi_i$ is the ith isobar divided difference.

        The vectors indexing the Grothendieck polynomials can as well been seen as lehmer codes.

        INPUT:

        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``G``) the basis representation for elements

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the Grothendieck basis
          of type ``group_type`` with positive exposants indexed by vectors where ``R``
          is the algebra base ring defined in the abstract algebra.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: Groth = A.grothendieck_positive_basis(); Groth
            The Multivariate polynomial algebra on x over Rational Field on the Grothendieck basis of type A, with positive exposants
            sage: Groth.an_element()
            G[0, 0, 0] + 3*G[0, 1, 0] + 2*G[1, 0, 0] + G[2, 2, 3]
            sage: Groth[1,2,3]
            G[1, 2, 3]

        Let us see some coercions::

            sage: g = Groth.an_element(); g
            G[0, 0, 0] + 3*G[0, 1, 0] + 2*G[1, 0, 0] + G[2, 2, 3]
            sage: g.expand()
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 5*xA[1, 0, 0] - 3*xA[1, 1, 0] + xA[2, 2, 3] + xA[2, 3, 2] - xA[2, 3, 3] + xA[3, 2, 2] - xA[3, 2, 3] - xA[3, 3, 2] + xA[3, 3, 3]
            sage: xA = A.monomial_basis_with_type("A")
            sage: xA(g)
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 5*xA[1, 0, 0] - 3*xA[1, 1, 0] + xA[2, 2, 3] + xA[2, 3, 2] - xA[2, 3, 3] + xA[3, 2, 2] - xA[3, 2, 3] - xA[3, 3, 2] + xA[3, 3, 3]
            sage: x = A.monomial_basis()
            sage: x(g)
            x[0, 0, 0] + 3*x[0, 1, 0] + 5*x[1, 0, 0] - 3*x[1, 1, 0] + x[2, 2, 3] + x[2, 3, 2] - x[2, 3, 3] + x[3, 2, 2] - x[3, 2, 3] - x[3, 3, 2] + x[3, 3, 3]
            sage: xA.an_element(); Groth(xA.an_element())
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 2*xA[1, 0, 0] + xA[2, 2, 3]
            G[0, 0, 0] + 3*G[0, 1, 0] - G[1, 0, 0] + 3*G[1, 1, 0] + G[2, 2, 3] - G[2, 3, 2] + G[2, 3, 3] - G[3, 3, 2] + G[3, 3, 3]
            sage: x.an_element(); Groth(x.an_element())
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
            G[0, 0, 0] + 3*G[0, 1, 0] - G[1, 0, 0] + 3*G[1, 1, 0] + G[1, 2, 3] - G[1, 3, 2] + G[1, 3, 3] - G[2, 1, 3] + G[2, 2, 3] + G[2, 3, 1] - 2*G[2, 3, 2] + G[2, 3, 3] + G[3, 1, 2] - G[3, 1, 3] - G[3, 2, 1] + G[3, 2, 2] - G[3, 3, 2] + G[3, 3, 3] + G[4, 1, 1] - G[4, 1, 3]


        Let us see some operations::

            sage: Groth[1,2] + Groth[3,0,0]
            G[1, 2, 0] + G[3, 0, 0]
            sage: Groth.an_element() * Groth.an_element()
            G[0, 0, 0] + 6*G[0, 1, 0] + 9*G[0, 2, 0] + 4*G[1, 0, 0] + 21*G[1, 1, 0] - 9*G[1, 2, 0] + 16*G[2, 0, 0] - 12*G[2, 1, 0] + 2*G[2, 2, 3] + 6*G[2, 3, 3] + 6*G[2, 4, 2] - 6*G[2, 4, 3] + 4*G[3, 2, 3] + G[4, 4, 6] + G[4, 5, 5] - G[4, 5, 6]
            sage: Groth[1,2] * Groth[3,0,0]
            G[4, 2, 0] + G[5, 1, 0] - G[5, 2, 0]

        """

        from .linear_basis_on_vectors import GrothendieckPositiveBasisOnVectors
        if(basis_name is None):
            basis_name = "Grothendieck basis of type A, with positive exposants"
        if(self._show_main_var): basis_repr+= self._main_repr_var
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return GrothendieckPositiveBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr)


    def macdonald_basis_on_vectors(self, t1 =None , t2=None, q=None , basis_name = None, basis_repr = "M"):
        r"""
        Creates the the basis of non symmetric Macdonald polynomials indexed by vectors.

        INPUT:

        - ``t1``: (default: symbolic variable t1) the first parameter for the Hecke algebra operator
        - ``t2``: (default: symbolic variable t2) the second parameter for the Hecke algebra operator
        - ``q``: (default: symbolic variable q) the specific q parmater of the polynomials
        - ``basis_name``: (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``: (default: ``M``) the basis representation for elements

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the Macdolald basis
          of type A indexed by vectors where ``R`` is the algebra base ring defined
          in the abstract algebra.

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: K.<t1,t2,q> = QQ[]
            sage: K = K.fraction_field()
            sage: A = MultivariatePolynomialAlgebra(K)
            sage: Mac = A.macdonald_basis_on_vectors(); Mac
            The Multivariate polynomial algebra on x over Fraction Field of Multivariate Polynomial Ring in t1, t2, q over Rational Field on the Macdonald basis of type A (indexed by vectors)
            sage: Mac.an_element()
            M[0, 0, 0] + 3*M[0, 1, 0] + 2*M[1, 0, 0] + M[2, 2, 3]
            sage: Mac[1,2]
            M[1, 2]

        Let us see some coercions::

            sage: pol = Mac[1,2];pol
            M[1, 2]
            sage: pol.expand()
            t2^3*xA[0, 0] + ((t2^2*q+t2^2)/q)*xA[0, 1] + t2/q*xA[0, 2] + t2^2*xA[1, 0] + ((t2*q+t2)/q)*xA[1, 1] + 1/q*xA[1, 2]
            sage: xA = A.monomial_basis_with_type("A")
            sage: xA(pol)
            t2^3*xA[0, 0] + ((t2^2*q+t2^2)/q)*xA[0, 1] + t2/q*xA[0, 2] + t2^2*xA[1, 0] + ((t2*q+t2)/q)*xA[1, 1] + 1/q*xA[1, 2]
            sage: x = A.monomial_basis()
            sage: x(pol)
            t2^3*x[0, 0] + ((t2^2*q+t2^2)/q)*x[0, 1] + t2/q*x[0, 2] + t2^2*x[1, 0] + ((t2*q+t2)/q)*x[1, 1] + 1/q*x[1, 2]
            sage: Mac(x[1,0] + x[0,1])
            (t1-t2)*M[0, 0] + ((t1*q-t1)/(t1*q+t2))*M[0, 1] + (1/(-t2))*M[1, 0]


        Let us see some operations::

            sage: Mac[1,2] + Mac[1,0]
            M[1, 0] + M[1, 2]
            sage: Mac[1,2] * Mac[1,0]
            ((t1^2*t2*q^2-t1^2*t2*q-t2^3*q+t2^3)/(-t1*q-t2))*M[1, 2] + ((t1*q^2+t2*q^2)/(t1*q+t2))*M[1, 3] + ((t1^2*t2*q^3-t1^2*t2*q^2-t2^3*q^2+t2^3*q)/(-t1^2*q^2-2*t1*t2*q-t2^2))*M[2, 2]
            """
        if(t1 is None): t1 = self.base_ring()(var('t1'))
        if(t2 is None): t2 = self.base_ring()(var('t2'))
        if(q is None): q = self.base_ring()(var('q'))
        from .linear_basis_on_vectors import MacdonaldBasisOnVectors
        if(basis_name is None):
            basis_name = "Macdonald basis of type A (indexed by vectors) "
        if(self._show_main_var): basis_repr+= self._main_repr_var
        monomial_basis_with_type = self.monomial_basis_with_type("A")
        return MacdonaldBasisOnVectors(self, monomial_basis_with_type, basis_name, basis_repr, t1, t2,q)

    def _create_morphism(self,f1,f2):
        r"""
        Creates a morphism between two `FiniteRankMultivariatePolynomialAlgebra` on their `FiniteMonomialBasis`
        by adding extra variables to the elements of the one with least variables.
        The morphism is then registered as a coercion.

        This method is called by `finite_polynomial_ring` each time a new `FinitePolynomialRing` is created

        INPUT:
            - ``f1`` a `FiniteRankMultivariatePolynomialAlgebra`
            - ``f2`` another `FiniteRankMultivariatePolynomialAlgebra` with a different number of variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: F1 = A.algebra_finite_nb_variables(1); F1
            The Multivariate polynomial algebra on x over Rational Field with 1 variable
            sage: F2 = A.algebra_finite_nb_variables(2); F2
            The Multivariate polynomial algebra on x over Rational Field with 2 variables
            sage: M1 = F1.monomial_basis();M1
            The Multivariate polynomial algebra on x over Rational Field with 1 variable on the monomial basis
            sage: M2 = F2.monomial_basis(); M2
            The Multivariate polynomial algebra on x over Rational Field with 2 variables on the monomial basis
            sage: M2(M1.an_element())
            x[0, 0] + 3*x[1, 0] + 3*x[2, 0]

            by creating ``F1`` and ``F2`` through ``A`` a coercion between their monomial basis ``M1`` and ``M2``
            has been created

        """
        if(f1.nb_variables() > f2.nb_variables()):
            temp = f1
            f1 = f2
            f2 = temp



    def change_nb_variables(self, pol, nb_variables):
        r"""
        Forcing the addition of variables

        INPUT:

        - ``pol``: a polynomial expressed in any basis
        - ``nb_variables``: the new number of variables

        OUTPUT:

        - the polynomial seen as a polynomial of ``nb_variables`` variables

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: pol = x.an_element(); pol
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
            sage: A.change_nb_variables(pol, 5)
            x[0, 0, 0, 0, 0] + 3*x[0, 1, 0, 0, 0] + 2*x[1, 0, 0, 0, 0] + x[1, 2, 3, 0, 0]
            sage: xA = A.monomial_basis_with_type("A")
            sage: pol = xA.an_element(); pol
            xA[0, 0, 0] + 3*xA[0, 1, 0] + 2*xA[1, 0, 0] + xA[2, 2, 3]
            sage: A.change_nb_variables(pol, 5)
            xA[0, 0, 0, 0, 0] + 3*xA[0, 1, 0, 0, 0] + 2*xA[1, 0, 0, 0, 0] + xA[2, 2, 3, 0, 0]
        """
        if(nb_variables < pol.nb_variables()):
            raise NotImplementedError("This method doesn't reduce the number of variables, use reduce_nb_variables")
        basis = pol.parent().basis_tower().finite_rank_basis(nb_variables)
        return basis( pol )

    def reduce_nb_variables(self, pol):
        """
        Creates a polynomial by removing all last variables with exposant 0 of the given polynomial

        INPUT:

        - ``pol``: the polynomial to be reduced

        OUTPUT:

        - a polynomial equal to ``pol`` and without all the last variables
          with exposant 0

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: pol = x[1,2,2,0,0] + x[2,5];pol
            x[1, 2, 2, 0, 0] + x[2, 5, 0, 0, 0]
            sage: pol.nb_variables()
            5
            sage: red = A.reduce_nb_variables(pol); red
            x[1, 2, 2] + x[2, 5, 0]
            sage: red.nb_variables()
            3
            sage: red == pol
            True
        """
        max = 1
        for ind, coeff in pol:
            for i in range(pol.nb_variables()-1,-1,-1):
                if(ind[i]!=0):
                    if(i+1>max): max = i+1
                    break
        if(max==pol.nb_variables()): return pol
        codomain = pol.parent().basis_tower().finite_rank_basis(max)
        return sum( [coeff * codomain([ind[i] for i in range(0,max)]) for ind, coeff in pol] )

    def maxDiffDiv(self, pol):
        """
        Apply the maximum divided difference to the polynomial.
        As the result is a symmetrical function, it is writen as a symmetrical
        Schubert polynomial (same as Schur functions). Giving the result in this
        basis makes the algorithm faster and the result compact.

        INPUT:

        - ``pol``: the polynomial to apply the maximum divided difference on

        OUTPUT:

        - the result polynomial in Schubert basis after applying maximum
          divided difference

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: pol = x[1,2,3]
            sage: A.maxDiffDiv(pol)
            -Y[1, 1, 1]
        """
        return pol.maxDiffDiv()

    def maxPi(self, pol):
        """
        Apply the maximum isobaric divided difference to the polynomial.
        As the result is a symmetrical function, it is writen as a symmetrical
        Schubert polynomial (same as Schur functions). Giving the result in this
        basis makes the algorithm faster and the result compact.

        INPUT:

        - ``pol``: the polynomial to apply the maximum isobaric divided
          difference on

        OUTPUT:

        - the result polynomial in Schubert basis after applying maximum
          isobaric divided difference

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
            sage: pol = x[2,1,3]
            sage: A.maxPi(pol)
            -Y[2, 2, 2]
            """
        return pol.maxPi()

class FiniteRankMultivariatePolynomialAlgebra(UniqueRepresentation, Parent):
    r"""
    This class implements the polynomial algebra in a given number of variables.

    INPUT:

    - ``polynomial_ring_tower`` -- the class of the polynomial algebra in an
      unset number of variables. from which the ``FiniteRankMultivariatePolynomialAlgebra``
      comes from. A ``FiniteRankMultivariatePolynomialAlgebra`` always comes from a
      ``MultivariatePolynomialAlgebra`` which contains general informations like
      the base ring.
    - ``nb_variables`` -- the number of variables
    - ``main_repr_var`` -- the letter corresponding to the set of variables,
      it is used to represent several bases, default is ``x``


    TESTS::

        sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
        sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
        sage: var('t')
        t
        sage: A = MultivariatePolynomialAlgebra(QQ[t])
        sage: B = FiniteRankMultivariatePolynomialAlgebra(A,4)
        sage: B
        The Multivariate polynomial algebra on x over Univariate Polynomial Ring in t over Rational Field with 4 variables
        sage: TestSuite(B).run()
    """
    def __init__(self, polynomial_ring_tower, nb_variables, main_repr_var = 'x', extra_bases_category = None):
        r"""
        TESTS::

            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteRankMultivariatePolynomialAlgebra(A,4)
        """
        self._polynomial_ring_tower = polynomial_ring_tower
        self._nb_variables = nb_variables
        Parent.__init__(
            self,
            base = polynomial_ring_tower.base_ring(),
            category = GradedAlgebras(polynomial_ring_tower.base_ring()).WithRealizations()
        )
        self._main_repr_var = main_repr_var
        self._polynomial_ring_tower._register_finite_ring(self)
        self._extra_bases_category = extra_bases_category
        m = SetMorphism( Hom(self, polynomial_ring_tower), lambda x: x)
        m.register_as_coercion()

    def _repr_(self):
        r"""
        TESTS::

            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: FiniteRankMultivariatePolynomialAlgebra(A,4)
            The Multivariate polynomial algebra on x over Rational Field with 4 variables
        """
        if(self.nb_variables()>1):
            variables_str = "variables"
        else:
            variables_str = "variable"
        return "%s with %s %s"%(self.polynomial_ring_tower(),self.nb_variables(),variables_str)


    def nb_variables(self):
        r"""
        Return the number of variables of ``self``.

        EXAMPLES::

            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteRankMultivariatePolynomialAlgebra(A,8)
            sage: B.nb_variables()
            8
            sage: B = FiniteRankMultivariatePolynomialAlgebra(A,0)
            sage: B.nb_variables()
            0
        """
        return self._nb_variables

    def polynomial_ring_tower(self):
        r"""
        Return the polynomial ring tower given to define ``self``.

        EXAMPLES::

            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: K = CyclotomicField(3)
            sage: A = MultivariatePolynomialAlgebra(K)
            sage: B = FiniteRankMultivariatePolynomialAlgebra(A,3)
            sage: B.polynomial_ring_tower()
            The Multivariate polynomial algebra on x over Cyclotomic Field of order 3 and degree 2
        """
        return self._polynomial_ring_tower


    def _element_constructor_(self, element):
        r"""
        Construct an element of ``self``.

        As ``self`` is an abstract algebra, this method will
        just check if ``element`` belongs to ``self``

        INPUT:

        -``element`` the element to be contructed from

        OUTPUT:

        The element itself if it belongs to ``self``, if not,
        it raises a TypeError Exception

        TESTS::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteRankMultivariatePolynomialAlgebra(A,3)
            sage: p = B.an_element()
            sage: B._element_constructor_(p)
            x[0, 0, 0] + 3*x[0, 1, 0] + 2*x[1, 0, 0] + x[1, 2, 3]
            sage: B._element_constructor_(1)
            Traceback (most recent call last):
            ...
            ValueError: '1' is not an element of 'The Multivariate polynomial algebra on x over Rational Field with 3 variables'
        """
        if self.is_parent_of(element):
            return element
        raise ValueError("'%s' is not an element of '%s'"%(element, self))


    def an_element(self):
        r"""
        Return an element of ``self``. By default, this element lies in
        the monomial basis.

        EXAMPLES::

            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteRankMultivariatePolynomialAlgebra(A,7)
            sage: B.an_element()
            x[0, 0, 0, 0, 0, 0, 0] + 3*x[0, 1, 0, 0, 0, 0, 0] + 2*x[1, 0, 0, 0, 0, 0, 0] + x[1, 2, 3, 4, 5, 6, 7]

        """
        return self.a_realization().an_element()

    def a_realization(self):
        r"""
        Returns a default realization of ``self``, the monomial basis

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: F = A.algebra_finite_nb_variables(3)
            sage: F.a_realization()
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the monomial basis
        """
        return self.monomial_basis()


    def monomial_basis(self, basis_repr = None):
        r"""
        Returns the algebra ``self`` view in the monomials basis.

        INPUT:

        - ``basis_repr``, the representation letter for the elements of the base, by default, it is the main representation for
          the set of variable : ``self._main_repr_var``

        EXAMPLES::

            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: B = FiniteRankMultivariatePolynomialAlgebra(A,5)
            sage: B.monomial_basis()
            The Multivariate polynomial algebra on x over Rational Field with 5 variables on the monomial basis
        """
        if(basis_repr is None): basis_repr = self._main_repr_var
        from .monomial import FiniteMonomialBasis
        return FiniteMonomialBasis(self, basis_repr,  extra_category = self._extra_bases_category)

    def monomial_basis_with_type(self, letter, basis_repr = None):
        r"""
        Return the algebra ``self`` view in the proper ambient space of the
        root system design by ``letter``.

        EXAMPLES::

            sage: from multipolynomial_bases.multivariate_polynomials import FiniteRankMultivariatePolynomialAlgebra
            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: F = FiniteRankMultivariatePolynomialAlgebra(A,2)
            sage: F.monomial_basis_with_type("B")
            The Multivariate polynomial algebra on x over Rational Field with 2 variables on the Monomial basis of type B
        """
        if(basis_repr is None): basis_repr = self._main_repr_var + letter
        from .ambient_space_basis import FiniteRankPolynomialRingWithBasisFromAmbientSpace
        if(letter == "A"): number = self.nb_variables()-1
        else: number = self.nb_variables()
        code = str(letter) + str(number)
        basis = FiniteRankPolynomialRingWithBasisFromAmbientSpace(self,code,letter,"Monomial basis of type " + letter, basis_repr,  extra_category = self._extra_bases_category)
        return basis

    def from_morphism_basis(self, polynomial_ring_tower, basis_name, basis_repr):
        r"""
        Creates a basis defined by its morphism to another basis

        INPUT:

        - ``polynomial_ring_tower``: the basis of ``AbsractPolynomialRing`` which is a facade to this basis and represents
          it on a undefined number of variables. It must have a `get_morphism_on_basis` method and a `get_basis_keys` method
          as well as a ``morphism_to_basis``
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements (exemple "x")

        OUTPUT:

        - the basis of which elements are indexed by the set return
          by ``polynomial_ring_tower.get_basis_keys(self.nb_variables())`` and can be coerced on
          the ``morphims_to_basis`` basis of the ``polynomial_ring_tower`` on the right number of variables (``self.nb_variables()``)

        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: M = A.monomial_basis()
            sage: def get_basis_keys(n): code = "A" + str(n-1); return RootSystem(code).ambient_space(QQ)
            sage: def get_morphism_on_basis(n): return lambda key: M( [key[i] for i in range(n)])
            sage: MyBasis = A.from_morphism_basis(1,M,get_basis_keys,get_morphism_on_basis,"My Basis", "X"); MyBasis
            The Multivariate polynomial algebra on x over Rational Field on the My Basis
            sage: from multipolynomial_bases.multivariate_polynomials  import FiniteRankMultivariatePolynomialAlgebra
            sage: F2 = FiniteRankMultivariatePolynomialAlgebra(A,2)
            sage: MyFiniteBasis = F2.from_morphism_basis(MyBasis,"MyBasis", "X"); MyFiniteBasis
            The Multivariate polynomial algebra on x over Rational Field with 2 variables on the MyBasis
            sage: MyFiniteBasis.an_element()
            X(2, 2)
            sage: M( MyFiniteBasis.an_element())
            x[2, 2]


        We have recreated the basis on ambient space.

        """

        from .basis import FiniteRankPolynomialRingWithBasisFromMorphism
        return FiniteRankPolynomialRingWithBasisFromMorphism(self, polynomial_ring_tower, basis_name, basis_repr, extra_category = self._extra_bases_category)


    def linear_basis_on_vectors(self, polynomial_ring_tower, basis_name, basis_repr, **keywords):
        r"""
        Creates a linear basis on objects inedexed by vectors based on an operation
        to convert each object (through its vector) into a ambient space basis polynomial.
        The type of the ambient space basis and the method of conversion are all contained into the
        ``polynomial_ring_tower``

        - ``polynomial_ring_tower``: the basis of ``AbsractPolynomialRing`` which is a facade to this basis and
          represents it on a undefined number of variables. It should inherit from a basis.LinearBasisOnVectors
        - ``basis_name``: the name of the basis (used in repr)
        - ``basis_repr``: the basis representation for elements
        - ``**keyword`` : parameters used to create the morphism to the ambient space basis,
          sent to ``CombinatorialFreeModule.module_morphism``.

        OUTPUT:

        - a basis named ``basis_name`` and defined by its conversion to an ambient space basis,
          the type of the ambient space basis and the method of conversion are all contained into the
          ``polynomial_ring_tower``


        EXAMPLES::

            sage: from multipolynomial_bases import MultivariatePolynomialAlgebra
            sage: A = MultivariatePolynomialAlgebra(QQ)
            sage: def schubert_on_basis(v, basis, call_back):
            ....:    for i in range(len(v)-1):
            ....:        if(v[i]<v[i+1]):
            ....:            v[i], v[i+1] = v[i+1] + 1, v[i]
            ....:            return call_back(v).divided_difference(i+1)
            ....:    return basis(v)
            sage: myBasis = A.linear_basis_on_vectors("A","MySchub","Y",schubert_on_basis)
            sage: F3 = A.algebra_finite_nb_variables(3)
            sage: myFiniteBasis = F3.linear_basis_on_vectors(myBasis,"MySchub","Y")
            sage: myFiniteBasis
            The Multivariate polynomial algebra on x over Rational Field with 3 variables on the MySchub

        """
        from .linear_basis_on_vectors import FiniteLinearBasisOnVectors
        return FiniteLinearBasisOnVectors(self, polynomial_ring_tower, basis_name, basis_repr, extra_category = self._extra_bases_category, **keywords)



