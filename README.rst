==================================
Bases of Multivariate polynomials
==================================

This package implements different bases of the algebra of multivariate polynomials
using `SageMath <http://www.sagemath.org>`_. This work has been started in 2010 and 
was distriuted in the Sage-combinat experimental package which does not exist anymore. 
That is why we know offer it as an extra Sage package. 

Installation
-------------

Requirement
~~~~~~~~~~~

sage 9.

Easy way
~~~~~~~~~

::

    $ sage -pip install multipolynomial_bases

From source
~~~~~~~~~~~

Download the source from the `github <https://github.com/VivianePons/multipolynomial-bases>`_
and run this command from the repo root::

    $ sage  -pip install --upgrade --no-index -v .

Equivalently, you can use the Makefile install command::

    $ make install

Usage
------

Once installed, you can use it in sage by importing the features::

    sage: from multipolynomial_bases import *
    sage: A.<x> = MultivariatePolynomialAlgebra(QQ)
    sage: A
    The Multivariate polynomial algebra on x over Rational Field

Documentation
-------------

The (partially written) documentation `is available here <http://openpyviv.com/multipolynomial-bases/>`_.

