#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
"""
Implementation of Raghavendra's linear systems solver
over finite fields.

See:

* http://www.eecs.berkeley.edu/~prasad/linsystems.pdf
* https://rjlipton.wordpress.com/2012/08/09/a-new-way-to-solve-linear-equations/

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
import os
import sys
import random  # used for non-repetitive random number

# local imports
from lib.modular import Modular
from lib.linalg import *
import lib.log

## test parameters
#
# vary these settings to test the algorithm over different base
# fields, matrix dimensions etc.

# be careful! already with q=3 we have N=7169 ...
q = 13

# vector length
n = 5

# a function returning a random number to be used as a vector
# coordinate; the default implementation chooses all coordinates
# uniformly from the integer interval [0,q-1].
random_elt = (lambda: random.randint(0, q-1))


## main program

log = lib.log.init()

# modular arithmetic modulo `q`
numtype = (lambda x: Modular(x, q))

def solve(A, b, N, sample_fn):
    """
    Return a list of solutions to the linear system `Ax = b`.

    :param A: a matrix with `m` rows and `n` columns, each being a `Modular(_,q)` instance.
    :param b: a vector of `n` elements, each element being a `Modular(_,q)` instance.
    :param sample_fn: a function returning a random number to be used as a vector coordinate
    """
    # Input:
    # * A is a matrix with m row vectors
    m = len(A)
    # * each row vector has n columns
    n = len(A[0])
    # b is a column vector with n rows
    assert m == len(b)
    # DEBUG
    log.debug("Given data:")
    log.debug("  m = %s", m)
    log.debug("  n = %s", n)
    log.debug("  A = %s", prettyprint_matrix(A, indent=8))
    log.debug("  b = %s", prettyprint_vector(b))
    log.debug("  N = %s", N)

    # Step 1. Sample N random vectors S_0 = {z_1, ... , z_N} where each z_i is i.i.d uniform.
    z = [ make_random_vector(n, numtype, sample_fn) for _ in range(N) ]
    #print "Initial choices of z's:"
    #for l, z_l in enumerate(z):
    #    print "  z_%d = %s" % (l, prettyprint_vector(z_l))

    # Step 2.
    for i in range(m): # for each constraint `A_i x = b_i` do:
        log.debug("Step 2, imposing %d-th constraint ...", i)
        # Step 2(a): selection
        t = [ z_l for z_l in z if (dot_product(A[i], z_l) == b[i]) ]
        # Step 2(b): recombination
        if len(t) == 0: # if T=\emptyset, then SYSTEM INFEASIBLE
            log.critical(("SYSTEM INFEASIBLE: No vector satisfies constraint A_%d", i))
            raise RuntimeError("SYSTEM INFEASIBLE: No vector satisfies constraint A_%d" % i)
        else:
            z = recombine(t, N)

    # final result
    return z


def recombine(t, N):
    """
    Recombination function.
    """
    s = [ ]
    for i in range(N):
        y = random.sample(t, q+1)
        # compute the sum of y's
        sum_y = [ 0 for _ in range(len(y[0])) ]
        for y_i in y:
            sum_y = vector_sum(sum_y, y_i)
        s.append(sum_y)
    return s


def _check_solution(A, b, z):
    for z_l in z:
        assert z_l == z[0]
    log.debug("OK: All final vectors are equal!")
    for i in range(len(b)):
        assert b[i] == dot_product(A[i], z[0])
    log.debug("OK: final vector(s) are a solution of the system Ax=b!")


def test_with_random_matrix(dim=5, N=None):
    A = make_random_matrix(dim, dim, numtype, random_elt)
    b = make_random_vector(dim, numtype, random_elt)

    n = len(b)
    if N is None:
        N = 1 + int(145*n*q*q*math.log(q))

    x = solve(A, b, N, random_elt)

    _check_solution(A, b, x)


def test_with_identity_matrix(dim, N=None):
    A = identity_matrix(dim, numtype)
    b = make_random_vector(dim, numtype, random_elt)

    n = len(b)
    if N is None:
        N = 1 + int(145*n*q*q*math.log(q))

    x = solve(A, b, N, random_elt)

    _check_solution(A, b, x)


if __name__ == '__main__':

    # Fix random numbers for debugging
    #random.seed(100)

    test_with_identity_matrix(n)
    test_with_random_matrix(n)
