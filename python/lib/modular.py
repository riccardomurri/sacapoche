#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple-minded implementation of modular arithmetic in Python.
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


def egcd(a, b):
    """
    Extended GCD.

    Given integers `a` and `b`, return a triple `(g, u, v)` such that
    `g` is the GCD of `a` and `b`, and the Bézout identity
    `u*a + v*b = g` holds.

    See http://code.activestate.com/recipes/474129-extended-great-common-divisor-function/#c4
    """
    u, u1 = 1, 0
    v, v1 = 0, 1
    g, g1 = a, b
    while g1:
        q = g // g1
        u, u1 = u1, u - q * u1
        v, v1 = v1, v - q * v1
        g, g1 = g1, g - q * g1
    return g, u, v


class Modular(object):
    """A simple implementation of modular arithmetic in Python.

    Initialize a `Modular` instance with integer value and modulus::

      >>> x = Modular(3, 13)
      >>> y = Modular(9, 13)

    The modulus is available as attribute `.modulus` on each `Modular`
    instance::

      >>> x.modulus
      13

    `Modular` instances support the four arithmetic operations and the
    equality comparison::

      >>> x + y == Modular(-1, 13)
      True
      >>> x - y
      Modular(7, 13)

    """

    __slots__ = ('value', 'modulus')

    def __init__(self, value, modulus):
        self.modulus = modulus
        self.value = (value % modulus)

    def __int__(self):
        return self.value

    def __long__(self):
        return self.value

    def __str__(self):
        return ("%d(mod %d)" % (self.value, self.modulus))

    def __repr__(self):
        return ("Modular(%d, %d)" % (self.value, self.modulus))

    def inverse(self):
        g, u, v = egcd(self.value, self.modulus)
        assert g == 1
        # now `u*value + v*modulus = 1`, hence `u` is the inverse
        return Modular(u, self.modulus)

    # comparison operators; only == and != make sense for integers mod p

    def __eq__(self, other):
        assert self.modulus == other.modulus
        return self.value == other.value

    def __ne__(self, other):
        assert self.modulus == other.modulus
        return self.value == other.value

    # arithmetic operations

    def __add__(self, other):
        assert self.modulus == other.modulus
        return Modular(self.value + other.value, self.modulus)

    def __sub__(self, other):
        assert self.modulus == other.modulus
        return Modular(self.value - other.value, self.modulus)

    def __mul__(self, other):
        assert self.modulus == other.modulus
        return Modular(self.value * other.value, self.modulus)

    def __div__(self, other):
        assert self.modulus == other.modulus
        return Modular(self.value * other.inverse().value, self.modulus)

    # arithmetic operations, in-place modifiers

    def __iadd__(self, other):
        assert self.modulus == other.modulus
        self.value += other.value
        self.value %= self.modulus

    def __isub__(self, other):
        assert self.modulus == other.modulus
        self.value -= other.value
        self.value %= self.modulus

    def __imul__(self, other):
        assert self.modulus == other.modulus
        self.value *= other.value
        self.value %= self.modulus

    def __idiv__(self, other):
        assert self.modulus == other.modulus
        self.value *= other.inverse().value
        self.value %= self.modulus

    # arithmetic operations, reversed order variant
    #
    # Note: these are called when __op__(x,y) has failed, so we know
    # that `other` is *not* a Modular instance

    def __radd__(self, other):
        return Modular(other + self.value, self.modulus)

    def __rsub__(self, other):
        return Modular(other - self.value, self.modulus)

    def __rmul__(self, other):
        return Modular(other * self.value, self.modulus)

    def __rdiv__(self, other):
        other_ = Modular(other, self.modulus)
        return (other_ / self)


# if invoked as a standalone module, run tests

if __name__ == '__main__':
    import doctest
    doctest.testmod(name='modular',
                    optionflags=doctest.NORMALIZE_WHITESPACE)
