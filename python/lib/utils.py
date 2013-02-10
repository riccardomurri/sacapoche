#! /usr/bin/env python
#
"""
Miscellaneous utility functions.
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


import random


def randomized(seq):
    """
    Return the elements of `seq` one by one in a random order.
    """
    already = set()
    l = len(seq)
    while True:
        i = random.randint(0, l-1)
        if i in already:
            continue
        yield seq[i]
        already.add(i)
        if len(already) == l:
            break
