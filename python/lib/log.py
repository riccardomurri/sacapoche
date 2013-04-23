#! /usr/bin/env python
"""
Auxiliary functions to set up Python's stdlib logging.
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

import logging
import os
import sys


# copy level names from logging
from logging import CRITICAL, ERROR, WARNING, INFO, DEBUG


def init(level=logging.INFO, name=None,
         fmt='%(name)s %(levelname)-8s: %(message)s',
         longfmt='%(asctime)s[%(relativeCreated)d] %(levelname)-8s %(name)s[%(module)s.%(funcName)s(), line %(lineno)d]: %(message)s'):
    if name is None:
        name = os.path.basename(sys.argv[0])
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    logfile = logging.FileHandler(name + '.debug.log')
    logfile.setLevel(logging.DEBUG)
    formatter = logging.Formatter(longfmt)
    logfile.setFormatter(formatter)
    # create console handler with a higher log level
    console = logging.StreamHandler()
    console.setLevel(level)
    formatter = logging.Formatter(fmt)
    console.setFormatter(formatter)
    # add the handlers to logger
    logger.addHandler(console)
    logger.addHandler(logfile)
    # return the logging object to the caller
    return logger
