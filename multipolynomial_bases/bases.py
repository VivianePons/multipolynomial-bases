# -*- coding: utf-8 -*-
r"""
Shortcuts for most used polynomial bases
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from multivariate_polynomials import MultivariatePolynomialAlgebra

def SchubertPolynomials(R, basis_name = None, basis_repr= "Y", **keywords):
    r"""
    The Schubert polynomial ring indexed by vectors as a basis of
    multivariate polynomials on different bases, see :class:`MultivariatePolynomialAlgebra`

    For double Schubert polynomials, see :class:`DoubleMultivariatePolynomialAlgebra`

    Here is the definition we use. For `v = (v_1, \cdots, v_n) \in
    \mathbb{N}^n`, we define

    `Y_v = x_1^{v_1}x_2^{v_2}\cdots x_n^{v_n}` if `v` is a partition,
    i.e, if `v_1 \geq v_2 \geq \cdots \geq v_n`.

    Otherwise, we have for ` v_i > v_{i+1}`

    `Y_{\cdots v_{i+1} v_i-1 \cdots} = Y_v \partial_i` where `\partial_i`
    is the ith divided difference.

    The vectors indexing the Schubert polynomials can as well been seen as
    lehmer codes.

    INPUT:

    - ``R`` -- a ring
    - ``basis_name`` -- (default: canonical name) the name of the basis
     (used in repr)
    - ``basis_repr`` -- (default: ``Y``) the basis representation for elements
    - ``**keywords`` --  other keywords to send t the abstract polynomial ring

    OUTPUT:

    - The Multivariate polynomial algebra on x over ``R`` on the Schubert basis

    EXAMPLES::

        sage: from multipolynomial_bases import SchubertPolynomials
        sage: Schub = SchubertPolynomials(QQ)
        sage: Schub
        The Multivariate polynomial algebra on x over Rational Field on the Schubert basis of type A
        sage: Schub.an_element()
        2*Y[1, 0, 0] + Y[2, 2, 3] + Y[0, 0, 0] + 3*Y[0, 1, 0]
        sage: Schub[1,2,3] + Schub[2,2]
        Y[2, 2, 0] + Y[1, 2, 3]

     some operations::

        sage: Schub[2,2,3]^2
        Y[4, 4, 6] + Y[4, 5, 5]
        sage: Schub[2,3] * Schub[1,2]
        Y[4, 4] + Y[3, 5]
        sage: Schub[2,3].expand()
        xA[3, 2] + xA[2, 3]
        sage: Schub[2,3].to_expr()
        x1^3*x2^2 + x1^2*x2^3
        sage: pol = Schub[3,2,3] + Schub[3,1,1]
        sage: pol.divided_difference(1)
        Y[1, 2, 1] + Y[2, 2, 3]
        sage: pol.divided_difference(2)
        0
        sage: pol.isobaric_divided_difference(1)
        Y[2, 3, 3] + Y[1, 3, 1] + Y[2, 4, 2]
        sage: pol.isobaric_divided_difference(2)
        Y[3, 2, 3] + Y[3, 1, 1]


    some coercions::

        sage: x = Schub.monomial_basis_with_type()
        sage: pol = x[2,2,3] + x[3,1,2]; pol
        xA[2, 2, 3] + xA[3, 1, 2]
        sage: Schub(pol)
        Y[2, 2, 3] - Y[2, 3, 2] + Y[3, 1, 2] - Y[3, 2, 1]
        sage: from multipolynomial_bases import DemazureHatPolynomials
        sage: Dem = DemazureHatPolynomials(QQ)
        sage: Schub(Dem[1,0,2])
        -Y[3, 0, 0] + Y[1, 0, 2] + Y[2, 1, 0] - Y[1, 2, 0] - Y[2, 0, 1]

    """
    A = MultivariatePolynomialAlgebra(R, **keywords)
    return A.schubert_basis(basis_name=basis_name, basis_repr=basis_repr)


def DemazurePolynomials(R, group_type ="A", basis_name = None, basis_repr = "K", **keywords):
    r"""
        Creates the Demazure polynomials where demazure / key polynomials are indexed
        by vectors, as basis of multivariate polynomials on different bases,
         see :class:`MultivariatePolynomialAlgebra`

        Here is the definition we use for type A. For `v = (v_1, \cdots, v_n) \in \mathbb{N}^n`,
        we define

        `K_v = x_1^{v_1}x_2^{v_2}\cdots x_n^{v_n}` if `v` is a partition, i.e, if `v_1 \geq v_2 \geq \cdots \geq v_n`.

        Otherwise, we have for ` v_i > v_{i+1}`

        `K_{\cdots v_{i+1} v_i \cdots} = K_v \pi_i` where `\pi_i` is the
        ith isobar divided difference.

        The vectors indexing the key polynomials can as well been seen
        as lehmer codes.

        INPUT:

        - ``R`` -- a ring
        - ``group_type`` -- (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name`` -- (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``-- (default: ``K``) the basis representation for elements
        - ``**keywords`` --  other keywords to send t the abstract polynomial ring

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the Demazure basis
          of type ``group_type``

        EXAMPLES::

            sage: from multipolynomial_bases import DemazurePolynomials
            sage: Dem = DemazurePolynomials(QQ); Dem
            The Multivariate polynomial algebra on x over Rational Field on the Demazure basis of type A
            sage: Dem.an_element()
            2*K[1, 0, 0] + K[2, 2, 3] + K[0, 0, 0] + 3*K[0, 1, 0]
            sage: Dem[1,2,2] + Dem[2,3,1]
            K[1, 2, 2] + K[2, 3, 1]

        some operations::

            sage: Dem[2,2,3]^2
            K[4, 4, 6] + K[4, 5, 5]
            sage: Dem[2,3] * Dem[1,2]
            K[4, 4] + K[3, 5]
            sage: Dem[2,3].expand()
            xA[3, 2] + xA[2, 3]
            sage: Dem[2,3].to_expr()
            x1^3*x2^2 + x1^2*x2^3
            sage: pol = Dem[3,2,3] + Dem[3,1,1]
            sage: pol.divided_difference(1)
            K[1, 2, 1] + K[2, 2, 3] - K[2, 3, 2]
            sage: pol.divided_difference(2)
            0
            sage: pol.isobaric_divided_difference(1)
            K[2, 3, 3] + K[1, 3, 1]
            sage: pol.isobaric_divided_difference(2)
            K[3, 2, 3] + K[3, 1, 1]

        some coercions::

            sage: x = Dem.monomial_basis_with_type(); x
            The Multivariate polynomial algebra on x over Rational Field on the Ambient space basis of type A
            sage: pol = x[2,2,3] + x[3,1,2]; pol
            xA[2, 2, 3] + xA[3, 1, 2]
            sage: Dem(pol)
            K[2, 2, 3] - K[2, 3, 2] + K[3, 1, 2] - K[3, 2, 1]
            sage: from multipolynomial_bases import SchubertPolynomials
            sage: Schub = SchubertPolynomials(QQ)
            sage: Dem(Schub[1,0,2])
            K[1, 0, 2] + K[3, 0, 0]

    """
    A = MultivariatePolynomialAlgebra(R, **keywords)
    return A.demazure_basis(group_type=group_type, basis_name = basis_name, basis_repr=basis_repr)

def DemazureHatPolynomials(R, group_type ="A", basis_name = None, basis_repr = "^K", **keywords):
    r"""
        Creates the Demazure hat polynomials where demazure / key polynomials are indexed
        by vectors, as basis of multivariate polynomials on different bases,
         see :class:`MultivariatePolynomialAlgebra`

        Here is the definition we use for type A. For `v = (v_1, \cdots, v_n)
        \in \mathbb{N}^n`, we define

        `\hat{K}_v = x_1^{v_1}x_2^{v_2}\cdots x_n^{v_n}` if `v` is a partition, i.e,
         if `v_1 \geq v_2 \geq \cdots \geq v_n`.

        Otherwise, we have for ` v_i > v_{i+1}`

        `\hat{K}_{\cdots v_{i+1} v_i \cdots} = \hat{K}_v \hat{\pi}_i` where
         `\hat{\pi}_i` is the ith isobar hat divided difference.

        The vectors indexing the key polynomials can as well been seen
        as lehmer codes.

        INPUT:

        - ``R`` -- a ring
        - ``group_type`` -- (default: ``A``) the letter that represents the type of the weyl group
        - ``basis_name`` -- (default: canonical name) the name of the basis (used in repr)
        - ``basis_repr``-- (default: ``^K``) the basis representation for elements
        - ``**keywords`` --  other keywords to send t the abstract polynomial ring

        OUTPUT:

        - The Multivariate polynomial algebra on x over ``R`` on the Demazure basis
          of type ``group_type``

        EXAMPLES::

            sage: from multipolynomial_bases import DemazureHatPolynomials
            sage: HatDem = DemazureHatPolynomials(QQ); HatDem
            The Multivariate polynomial algebra on x over Rational Field on the Demazure hat basis of type A
            sage: HatDem.an_element()
            2*^K[1, 0, 0] + ^K[2, 2, 3] + ^K[0, 0, 0] + 3*^K[0, 1, 0]
            sage: HatDem[1,2,2] + HatDem[2,3,1]
            ^K[1, 2, 2] + ^K[2, 3, 1]

        some operations::

            sage: HatDem[2,2,3]^2
            ^K[4, 4, 6] - ^K[5, 4, 5] - ^K[4, 5, 5]
            sage: HatDem[2,3] * HatDem[1,2]
            -^K[4, 4] + ^K[3, 5]
            sage: HatDem[2,3].expand()
            xA[2, 3]
            sage: HatDem[2,3].to_expr()
            x1^2*x2^3
            sage: pol = HatDem[3,2,3] + HatDem[3,1,1]
            sage: pol.divided_difference(1)
            ^K[1, 2, 1] + ^K[2, 2, 3] + ^K[2, 1, 1]
            sage: pol.divided_difference(2)
            -^K[3, 2, 2]
            sage: pol.isobaric_divided_difference(1)
            ^K[3, 2, 3] + ^K[1, 3, 1] + ^K[2, 3, 3] + ^K[3, 1, 1]
            sage: pol.isobaric_divided_difference(2)
            ^K[3, 1, 1]
            sage: pol.hat_isobaric_divided_difference(1)
            ^K[2, 3, 3] + ^K[1, 3, 1]
            sage: pol.hat_isobaric_divided_difference(2)
            -^K[3, 2, 3]

        some coercions::

            sage: x = HatDem.monomial_basis_with_type()
            sage: pol = x[2,2,3] + x[3,1,2]; pol
            xA[2, 2, 3] + xA[3, 1, 2]
            sage: HatDem(pol)
            ^K[2, 2, 3] + ^K[3, 1, 2]
            sage: from multipolynomial_bases import SchubertPolynomials
            sage: Schub = SchubertPolynomials(QQ)
            sage: HatDem(Schub([1,0,2]))
            ^K[3, 0, 0] + ^K[1, 0, 2] + ^K[2, 1, 0] + ^K[1, 2, 0] + ^K[2, 0, 1]

    """
    A = MultivariatePolynomialAlgebra(R, **keywords)
    return A.demazure_hat_basis(group_type=group_type, basis_name = basis_name, basis_repr=basis_repr)

def GrothendieckPolynomials(R, basis_name = None, basis_repr= "G", **keywords):
    r"""
    The Grothendieck polynomial ring indexed by vectors as a basis of
    multivariate polynomials on different bases, see :class:`MultivariatePolynomialAlgebra`

    For double Grothendieck polynomials, see :class:`DoubleMultivariatePolynomialAlgebra`.

   Here is the definition we use. For `v = (v_1, \cdots, v_n) \in \mathbb{N}^n`, we define

    `G_v = x_1^{v_1}x_2^{v_2}\cdots x_n^{v_n}` if `v` is a partition, i.e, if `v_1 \geq v_2 \geq \cdots \geq v_n`.

    Otherwise, we have for ` v_i > v_{i+1}`

    `G_{\cdots v_{i+1} v_i-1 \cdots} = G_v \pi_i` where `\pi_i` is the ith isobar divided difference.

    The vectors indexing the Grothendieck polynomials can as well been seen as lehmer codes.


    The vectors indexing the Grothendieck polynomials can as well been seen as
    lehmer codes.

    INPUT:

    - ``R`` -- a ring
    - ``basis_name`` -- (default: canonical name) the name of the basis
     (used in repr)
    - ``basis_repr`` -- (default: ``G``) the basis representation for elements
    - ``**keywords`` --  other keywords to send to the abstract polynomial ring

    OUTPUT:

    - The Multivariate polynomial algebra on x over ``R`` on the Grothendieck basis

    EXAMPLES::

        sage: from multipolynomial_bases import GrothendieckPolynomials
        sage: Groth = GrothendieckPolynomials(QQ)
        sage: Groth
        The Multivariate polynomial algebra on x over Rational Field on the Grothendieck basis of type A, with positive exposants
        sage: Groth.an_element()
        2*G[1, 0, 0] + G[2, 2, 3] + G[0, 0, 0] + 3*G[0, 1, 0]
        sage: Groth[1,2,3] + Groth[2,2]
        G[2, 2, 0] + G[1, 2, 3]

     some operations::

        sage: Groth[2,2,3]^2
        G[4, 4, 6] - G[4, 5, 6] + G[4, 5, 5]
        sage: Groth[2,3] * Groth[1,2]
        -G[4, 5] + G[4, 4] + G[3, 5]
        sage: Groth[2,3].to_expr()
        -x1^3*x2^3 + x1^3*x2^2 + x1^2*x2^3
        sage: pol = Groth[3,2,3] + Groth[3,1,1]
        sage: pol.divided_difference(1)
        G[2, 2, 1] + G[1, 2, 1] + G[2, 2, 3]
        sage: pol.divided_difference(2)
        0
        sage: pol.isobaric_divided_difference(1)
        G[3, 3, 3] - G[3, 4, 3] - G[2, 4, 3] + G[2, 3, 3] + G[1, 3, 1] + G[2, 4, 2] + G[3, 3, 1] + G[3, 4, 2] + G[2, 3, 1]
        sage: pol.isobaric_divided_difference(2)
        G[3, 2, 3] + G[3, 1, 1]

    some coercions::

        sage: x = Groth.monomial_basis_with_type()
        sage: pol = x[2,2,3] + x[3,1,2]; pol
        xA[2, 2, 3] + xA[3, 1, 2]
        sage: Groth(pol)
        G[3, 3, 3] + G[3, 2, 2] + G[2, 3, 3] - G[2, 3, 2] - G[3, 3, 2] + G[2, 2, 3] - G[3, 2, 1] + G[3, 1, 2]
        sage: from multipolynomial_bases import SchubertPolynomials
        sage: Schub = SchubertPolynomials(QQ)
        sage: Groth(Schub[1,0,2])
        G[3, 0, 2] + G[2, 2, 2] + G[1, 2, 2] + G[2, 0, 2] + G[1, 1, 2] + G[1, 0, 2]


    """
    A = MultivariatePolynomialAlgebra(R, **keywords)
    return A.grothendieck_positive_basis(basis_name=basis_name, basis_repr=basis_repr)

