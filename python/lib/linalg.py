#! /usr/bin/env python
#
"""
Basic linear algebra routines on top of Python's ``list`` type.

A vector is represented as a list of numbers; a matrix is a list of
vectors.

"""
# Copyright (C) 2012-2013 University of Zurich. All rights reserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
__version__ = '$Revision$'
__author__ = 'Riccardo Murri <riccardo.murri@gmail.com>'
__docformat__ = 'reStructuredText'


import math


def vector_sum(x, y):
    """Return the sum of two vectors `x` and `y`."""
    assert len(x) == len(y)
    return [ x[i]+y[i] for i in range(len(x)) ]


def dot_product(x, y):
    """Return the dot product of two vectors `x` and `y`."""
    assert len(x) == len(y)
    return sum((x[i]*y[i]) for i in range(len(x)))


def scalar_vector_product(alpha, x):
    """Given a scalar `alpha` and a vector `x`, return `alpha * x`."""
    return [ alpha*x_l for x_l in x ]


def matrix_vector_product(M, x):
    """
    Return the product of matrix `M` and vector `x`.
    """
    assert len(M) == 0 or len(M[0]) == len(x)
    return [ dot_product(m, x) for m in M ]


def norm(x):
    """
    Return the Euclidean norm of vector `x`.
    """
    return math.sqrt(sum(x[i]*x[i] for i in range(len(x))))


def identity_matrix(N, numtype=float):
    """
    Return the N by N identity matrix.
    """
    return [
        [ numtype(1 if (i==j) else 0) for j in range(N) ]
        for i in range(N)
    ]


def make_random_vector(n, numtype, sample_fn):
    """
    Return a random vector of length `n`.

    Entries are sampled from the given distribution function.

    :param int n:
      Length of the vector

    :param sample_fn:
      A callable that returns a single numeric value.  Should require no parameter.
    """
    assert n > 0
    return [ numtype(sample_fn()) for _ in range(n) ]


def make_random_matrix(m, n, numtype, sample_fn):
    """
    Return a matrix with `m` rows and `n` columns.
    """
    return [ make_random_vector(n, numtype, sample_fn) for _ in range(m) ]


def prettyprint_vector(x):
    return '[ ' + str.join(" ", (str(x_l) for x_l in x)) + ' ]'


def prettyprint_matrix(A, indent=4):
    return ('['
            + str.join('\n' + (" " * indent),
                       [prettyprint_vector(A_l) for A_l in A])
            + ']')
