#! /usr/bin/env python
#
"""
Implementation of Fliege's linear system solver, as described in:

  *A Randomized Parallel Algorithm with Run Time `O(n^2)`
  for Solving an `n` by `n` System of Linear Equations,*
  Joerg Fliege, arXiv:1209.3995v1.

The code in this program uses exact arithmetic (fractions of
arbitrary-precision integers) to avoid numerical stability issues.

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
__author__ = 'Benjamin Jonen <benjamin.jonen@bf.uzh.ch>, Riccardo Murri <riccardo.murri@gmail.com>'
__docformat__ = 'reStructuredText'


import math
import os
import sys
import random  # used for non-repetitive random number

import numpy as np
import numpy.linalg

from lib.linalg import *


## test parameters
#
# vary these settings to test the algorithm over different base
# fields, matrix dimensions etc.

# vector length
n = 5

# a function returning a random number to be used as a vector
# coordinate; the default implementation chooses all coordinates
# uniformly from the specified integer interval.
random_elt = (lambda: random.randint(0, 16))

# set `numtype` to the a function that takes a (small) integer number
# and returns it converted to a chosen exact/arbitrary-precision
# numeric type.  The constructors of the Python standard
# `fractions.Fraction` and `fractions.Decimal` types are both OK.
import fractions
numtype = fractions.Fraction


## main algorithm

def solve(A, b, sample_fn):
    """
    Return a list of (approximate) solutions to the linear system `Ax = b`.

    :param A: a NumPy `m` times `n` 2D-array, representing a matrix with `m` rows and `n` columns
    :param b: a Numpy 1D-array real-valued vector of `n` elements
    :param sample_fn: a function returning a random number to be used as a vector coordinate
    """
    # Input:
    # * A is a with m row vectors
    m = len(A)
    # * each row vector has n columns
    n = len(A[0])
    # b is a column vector with n rows
    assert m == len(b)
    # consistency checks
    assert n**2/n>4, "n is not large enough. Increase eq. system to fulfill n**2>4n"
    # DEBUG
    print "Given data:"
    print "  m = ", m
    print "  n = ", n
    print "  A = ", prettyprint_matrix(A, indent=8)
    print "  b = ", prettyprint_vector(b)

    # Step 2. Generate random sample of normals (fulfilling Assumption 1.)
    v = [ make_random_vector(n, numtype, sample_fn) for _ in range(n+1) ]
    print "Initial choices of v's:"
    for l, v_l in enumerate(v):
        print "  v_%d = %s" % (l, prettyprint_vector(v_l))

    # save all (i,j) pairs for later -- this is invariant in the Step 3 loop
    all_ij_pairs = [ (i,j) for i in range(n+1) for j in range(n+1) if i<j ]
    #print 'potential ij_pairs = ', all_ij_pairs
    assert len(all_ij_pairs) == ((n+1)*n / 2)

    # Step 3.
    for k in range(m): # loop over eqs
        print "Step 3, iteration %d starting ..." % k
        # Step 3(a): choose n+1 random pairs (i_l, j_l) with i_l < j_l
        # pick randomly from all_ij_pairs set to get n+1 random pairs
        ij_pairs = random.sample(all_ij_pairs, n+1)
        assert len(ij_pairs) == n+1
        print '  ij_pairs = ', ij_pairs
        # Step 3(b): recombine the v's; use x's as temporary storage
        try:
            x = [ rec(v[i], v[j], A[k], b[k])
                  for l, (i, j) in enumerate(ij_pairs)
              ]
        except ZeroDivisionError:
            print ('**** STOP in solve() ****')
            for l in range(len(x)):
                print ("x_%d = %s" % (l, prettyprint_vector(x[l])))
            raise
        # Step 3(d): rename x's -> v's
        for l in range(n+1):
            v[l] = x[l]
        # check progress
        _check_distance(A, b, v, k)

    # final result
    return v


def rec(u, v, a, beta):
    """
    Recombination function.
    """
    assert len(u) == len(v)
    assert len(v) == len(a)
    t0 = beta - dot_product(a, v)
    t1 = dot_product(a, [(u[i] - v[i]) for i in range(len(u))])
    if t1 == 0:
        print ("**** STOP in rec() ****")
        print ("u = %s" % prettyprint_vector(u))
        print ("v = %s" % prettyprint_vector(v))
        print ("a = %s" % prettyprint_vector(a))
        raise ZeroDivisionError()
    t = t0 / t1
    return  [ (t*u[i] + (1-t)*v[i]) for i in range(len(u)) ]


def _check_distance(A, b, vs, k=None):
    if k is None:
        k = len(b)
    if k == len(b):
        print ("Final values of v's:")
    else:
        print ("Values of v's after step %d" % k)
    for l, v_l in enumerate(vs):
      print "  v_%d = %s" % (l, prettyprint_vector(v_l))
    print "Distances of solutions computed by Fliege's algorithm:"
    for l, v_l in enumerate(vs):
        dist = norm([
            matrix_vector_product(A,v_l)[i] - b[i]
            for i in range(k)
        ])
        print ("  |Av_%d - b| = %g" % (l, dist))

    if k == len(b):
        print "Numpy's `linalg.solve` solution:"
        A_ = np.array([ [ float(A[i][j])
                          for j in range(len(A[0])) ]
                        for i in range(len(A)) ],
                      ndmin=2)
        b_ = np.array([ float(b[k]) for k in range(len(b)) ])
        v_prime = np.linalg.solve(A_, b_)
        print "  v' = %s" % v_prime
        print "Distance of Numpy's `linalg.solve` solution:"
        dist_prime = np.linalg.norm(np.dot(A_,v_prime) - b_)
        print ("  |Av' - b| = %g" % dist_prime)


def test_with_random_matrix(dim=5):
    A = np.random.randint(low=1, high=9,size=(dim, dim))
    b = np.random.randint(low=1, high=9,size=(dim,))

    rank = np.linalg.matrix_rank(A)
    assert rank == dim, 'Matrix A needs to have full rank. '

    x = solve(A, b, random_elt)

    _check_distance(A, b, x)


def test_with_identity_matrix(dim=5):
    A = identity_matrix(dim, numtype)
    b = make_random_vector(dim, numtype, random_elt)

    x = solve(A, b, random_elt)

    _check_distance(A, b, x)


if __name__ == '__main__':

    # Fix random numbers for debugging
    #np.random.seed(100)

    # Set numpy print options
    np.set_printoptions(linewidth=300, precision=5, suppress=True)

    test_with_identity_matrix()
    test_with_random_matrix()
