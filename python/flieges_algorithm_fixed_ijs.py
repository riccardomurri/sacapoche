#! /usr/bin/env python
#
"""
Implementation of Fliege's linear system solver, with the
following modifications:

- a fixed set of `(i.j)` pairs is used instead of drawing a new random
  set at each iteration;

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


import os
import sys
import random  # used for non-repetitive random number

import numpy as np
import numpy.linalg

import lib.log

## test parameters
#
# vary these settings to test the algorithm over different base
# fields, matrix dimensions etc.

# number of equations
dim = 5

# vector length
n = 5

# a function returning a random vector of length `n`;
# the default implementation chooses all coordinates
# uniformly from the flotaing-point range [0,1].
random_vector = np.random.rand


## main algorithm

log = lib.log.init()

# floating point constants (machine-dependent)
fpinfo = np.finfo(float)
EPSILON = fpinfo.eps
TINY = fpinfo.tiny
LARGE = fpinfo.max / 2.

def solve(A, b, sample_fn=random_vector):
    """
    Return a list of (approximate) solutions to the linear system `Ax = b`.

    :param A: a NumPy `m` times `n` 2D-array, representing a matrix with `m` rows and `n` columns
    :param b: a Numpy 1D-array real-valued vector of `n` elements
    :param sample_fn: a function returning a random 1D-vector; the default implementation chooses all coordinates uniformly from the interval [0,1].
    """
    # Input:
    # * A is a with m row vectors
    m = A.shape[0]
    # * each row vector has n columns
    n = A.shape[1]
    # b is a column vector with n rows
    assert m == b.shape[0]
    # consistency checks
    #assert n**2/n>4, "n is not large enough. Increase eq. system to fulfill n**2>4n"
    # DEBUG
    log.debug("Given data:")
    log.debug("  m = %s", m)
    log.debug("  n = %s", n)
    log.debug("  A = %s", str.join("\n       ", str(A).split('\n')))
    log.debug("  b = %s", b)

    # Step 2. Generate random sample of normals (fulfilling Assumption 1.)
    v = [ sample_fn(n) for _ in range(n+1) ] # n+1 random vectors of n elements each
    log.debug("Initial choices of v's:")
    for l, v_l in enumerate(v):
        log.debug("  v_%d = %s", l, v_l)

    ij_pairs = [ (j % (n+1), (j+1) % (n+1)) for j in range(0,n+1) ]
    log.debug('ij_pairs = %s', ij_pairs)
    assert len(ij_pairs) == len(v)

    # Step 3.
    for k in range(m): # loop over eqs
        log.debug("Step 3, iteration %d starting ...", k)
        # Step 3(a): CANCELED, we now use a fixed set of indices
        log.debug('  ij_pairs = %s', ij_pairs)
        # Step 3(b): recombine the v's; use x's as temporary storage
        x = [ rec(v[i], v[j], A[k,:], b[k])
              for l, (i, j) in enumerate(ij_pairs)
        ]
        # Step 3(d): rename x's -> v's
        for l in range(n+1):
            v[l] = x[l]
        # check progress
        _check_distance(A, b, v, k)

    # final result
    return v


def rec(u, v, a, beta, q=0):
    """
    Recombination function, as defined in Fliege (2012).

    If the optional parameter `q` is > 0, then ensure that the
    denominator in the `t` factor is larger than `q`, by substituting
    `v` for a random convex combination of `u` and `v`.

    """
    t0 = beta - np.dot(a, v)
    t1 = np.dot(a, (u - v))
    #assert np.abs(t1) > EPSILON
    while q > 0 and abs(t1) < q:
        # see remark 8 on page 7
        c = 1. + t0 * random.random()
        v = u + c*(v - u)
        t1 = np.dot(a, (u - v))
    if t1 == 0:
        raise ZeroDivisionError()
    elif t1 < TINY:
        log.warning("tiny denominator in rec(%s, %s, %s, %s)", u,v,a,beta)
    t = t0 / t1
    return  (t * u + (1. - t) * v)


def _check_distance(A, b, vs, k=None):
    if k is None:
        k = len(b)
    if k == len(b):
        log.debug("Final values of v's:")
    else:
        log.debug("Values of v's after iteration %d:", k)
    for l, v_l in enumerate(vs):
        log.debug("  v_%d = %s", l, v_l)
    if k == len(b):
        log.debug("Distances of solutions computed by Fliege's algorithm:")
    else:
        log.debug("Accuracy of solutions of the first %d equations:", (k+1))
    for l, v_l in enumerate(vs):
        dist = np.linalg.norm(np.dot(A[0:k,:],v_l) - b[0:k])
        log.debug("  |Av_%s - b| = %g", l, dist)

    if k == len(b):
        log.debug("Numpy's `linalg.solve` solution:")
        v_prime = np.linalg.solve(A,b)
        log.debug("  v' = %s", v_prime)
        log.debug("Distance of Numpy's `linalg.solve` solution:")
        dist_prime = np.linalg.norm(np.dot(A,v_prime) - b)
        log.debug("  |Av' - b| = %g", dist_prime)
    return dist


def test_with_random_matrix(dim=5):
    A = np.random.randint(low=1, high=9,size=(dim, dim))
    b = np.random.randint(low=1, high=9,size=(dim,))

    rank = np.linalg.matrix_rank(A)
    assert rank == dim, 'Matrix needs to have full rank. '

    sol = solve(A, b)

    _check_distance(A, b, sol)


def test_with_identity_matrix(dim=5):
    A = np.eye(dim)
    b = np.random.randint(low=1, high=9,size=(dim,))

    sol = solve(A, b)

    _check_distance(A, b, sol)

def test_with_fixed_v(dim=5,v=np.random.randn(dim + 1,n),N=1000):
    log.info('Starting test with %d trials' % N)
    log.debug('v = \n%s', v)
    def fixed_v(n):
        global ctr
        ctr += 1
        return v[ctr % (n +1)]
    for ix in range(N):
        global ctr
        ctr = -1

        A = np.random.randint(low=1, high=9,size=(dim, dim))
        b = np.random.randint(low=1, high=9,size=(dim,))

        sol = solve(A, b, sample_fn=fixed_v)

        dist = _check_distance(A, b, sol)
        if dist > 1.e-5:
            log.warning("Algorithm did not converge at trial %d; the 2-normed distance to solution is %15.10f", ix, dist)
            log.debug('A = %s', A)
            log.debug('b = %s', b)
    log.info('Completed test with %d trials', N)


if __name__ == '__main__':

    # Fix random numbers for debugging
    #np.random.seed(100)

    # Set numpy print options
    np.set_printoptions(linewidth=300, precision=5, suppress=True)

    #test_with_identity_matrix(n)
    #test_with_random_matrix(n)
    test_with_fixed_v()
