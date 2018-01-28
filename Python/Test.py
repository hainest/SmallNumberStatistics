#!/usr/bin/python

#########################################################################
# Test.py, July 21, 2012
# Copyright (c) 2012, Tim Haines
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see http://www.gnu.org/licenses/.
#
#########################################################################

from SmallNumberStatistics import binomialLimits, poissonLimits
import numpy
from sys import exit

#########################################################################
# NAME:
#       Test
#
# AUTHOR:
#       Tim Haines, thaines.astro@gmail.com
#
# PURPOSE:
#       This tests binomialLimits and poissonLimits of the
#       SmallNumberStatistics module for working functionality.
#
# CATEGORY:
#       Statistics and probability
#
# RETURNS:
#       The string 'OK' if all passed, else an error message.
#
# DEPENDENCIES:
#       SmallNumberStatistics.py
#       numpy
#
##########################################################################


nsuccess = 9
ntotal = 109

# make sure we are getting close to Gehrels
(u,l) = binomialLimits(nsuccess, ntotal)
if abs(u-0.118) > 0.1 or abs(l - 0.0561) > 0.1:
    exit("Binomial: Not close to Gehrels!\n")

# make sure default cl is working
(a, b) = binomialLimits(nsuccess, ntotal, 0.841344746068543)
if abs(u-a) > 1e-3 or abs(l-b) > 1e-3:
    exit("Binomial: The default cl is not working!\n");

# make sure the 1-sigma cl works
(a, b) = binomialLimits(nsuccess, ntotal, 1.0, sigma=1)
if abs(u-a) > 1e-3 or abs(l-b) > 1e-3:
    exit("Binomial: Not sigma is not working!\n");

# Test an extremum from Gehrels
nsuccess = 0;
ntotal = 100;
(u, l) = binomialLimits(nsuccess, ntotal, 3.0, sigma=1)
if abs(u-0.0639) > 1e-3 or l != 0:
    exit("Binomial: Not close to Gehrels - 2!\n")
    
# Test with list inputs
nsuccess = [5, 9];
ntotal = [10, 15];
(u, l) = binomialLimits(nsuccess, ntotal)
if not len(u) == 2 and len(l) == 2:
    exit("Binomial: list input produced insufficient output\n")

# Test with NumPy array inputs
(u, l) = binomialLimits(numpy.array(nsuccess), numpy.array(ntotal))
if not isinstance(u,numpy.ndarray) and len(u) == 2 and isinstance(l,numpy.ndarray) and len(l) == 2:
    exit("Binomial: NumPy array input produced insufficient output\n")

##########################################################################

k = 13

# make sure we are getting close to Gehrels
(u, l) = poissonLimits(k)
if abs(u-17.70) > 1e-2 or abs(l - 9.441) > 1e-2:
    exit("Poisson: Not close to Gehrels!\n")

# make sure default cl is working
(a, b) = poissonLimits(k, 0.841344746068543)
if abs(u-a) > 1e-3 or abs(l - b) > 1e-3:
    exit("Poisson: The default cl is not working!\n")

# make sure the 1-sigma cl works
(a, b) = poissonLimits(k, 1.0, sigma=1)
if abs(u-a) > 1e-3 or abs(l - b) > 1e-3:
    exit("Poisson: Not sigma is not working!\n")

# Test an extremum from Gehrels
k = 100
(u, l) = poissonLimits(k, 3.0, sigma=1)
if abs(u - 133.8) > 0.1 or abs(l - 72.65) > 0.1:
    exit("Poisson: Not close to Gehrels - 2!\n")

# Test with list inputs
k = [5,15]
(u, l) = poissonLimits(k)
if not len(u) == 2 and len(l) == 2:
    exit("Poisson: list input produced insufficient output\n")

# Test with NumPy array inputs
(u, l) = poissonLimits(numpy.array(k))
if not isinstance(u,numpy.ndarray) and len(u) == 2 and isinstance(l,numpy.ndarray) and len(l) == 2:
    exit("Poisson: NumPy array input produced insufficient output\n")

print('OK')
