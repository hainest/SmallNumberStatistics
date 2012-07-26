#!/usr/bin/perl

#########################################################################
# Test.pl, July 21, 2012
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

use strict;
use warnings;
use Carp qw(croak);
use SmallNumberStatistics qw(binomialLimits poissonLimits);

#########################################################################
# NAME:
#       Test
#
# AUTHOR:
#       Tim Haines, tdhq36@mail.umkc.edu
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
#       SmallNumberStatistics.pm
#
##########################################################################

my ($nsuccess, $ntotal, $u, $l, $cl, $opt);
$opt = {'sigma' => 1};

$nsuccess = 9;
$ntotal = 109;

# make sure we are getting close to Gehrels
($u, $l) = binomialLimits($nsuccess, $ntotal);
if (abs($u-0.118) > 0.1 || abs($l - 0.0561) > 0.1)
{
    croak("Binomial: Not close to Gehrels!\n");
}

# make sure default cl is working
($a, $b) = binomialLimits($nsuccess, $ntotal, 0.841344746068543);
if ( abs($u-$a) > 1e-3 || abs($l-$b) > 1e-3)
{
    croak("Binomial: The default cl is not working!\n");
}

# make sure the 1-sigma cl works
($a, $b) = binomialLimits($nsuccess, $ntotal, 1.0, {'sigma' => 1});
if (abs($u-$a) > 1e-3 || abs($l - $b) > 1e-3)
{
    croak("Binomial: Not sigma is not working!\n");
}

# Test an extremum from Gehrels
$nsuccess = 0;
$ntotal = 100;
($u, $l) = binomialLimits($nsuccess, $ntotal, 3.0, {'sigma' => 1});
if (abs($u - 0.0639) > 1e-3 || $l != 0)
{
    croak("Binomial: Not close to Gehrels - 2!\n");
}

# Test the first calling sequence for scalar input
$u = binomialLimits($nsuccess, $ntotal);
unless (@$u == 1 && @{$u->[0]} == 2)
{
	croak("Binomial: list ref output bad for scalar input\n");
}

# Test the list ref input
$nsuccess = [5, 9];
$ntotal = [10, 15];
($u, $l) = binomialLimits($nsuccess, $ntotal);
unless (@$u == 2 && @$l == 2)
{
	croak("Binomial: ref input produced insufficient output\n");
}

# Test the first calling sequence for list ref input
$u = binomialLimits($nsuccess, $ntotal);
unless (@$u == 2 && @{$u->[0]} == 2)
{
    croak("Binomial: list ref output bad for list ref input\n");
}

##########################################################################

my $k = 13;

# make sure we are getting close to Gehrels
($u, $l) = poissonLimits($k);
if (abs($u-17.70) > 1e-2 || abs($l - 9.441) > 1e-2)
{
    croak("Poisson: Not close to Gehrels!\n");
}

# make sure default cl is working
($a, $b) = poissonLimits($k, 0.841344746068543);
if (abs($u-$a) > 1e-3 || abs($l - $b) > 1e-3)
{
    croak("Poisson: The default cl is not working!\n");
}

# make sure the 1-sigma cl works
($a, $b) = poissonLimits($k, 1.0, {'sigma' => 1});
if (abs($u-$a) > 1e-3 || abs($l - $b) > 1e-3)
{
    croak("Poisson: Not sigma is not working!\n");
}

# Test an extremum from Gehrels
$k = 100;
($u, $l) = poissonLimits($k, 3.0, {'sigma' => 1});
if (abs($u - 133.8) > 0.1 || abs($l - 72.65) > 0.1)
{
    croak("Poisson: Not close to Gehrels - 2!\n");
}

# Test the first calling sequence for scalar input
$u = poissonLimits($k);
unless (@$u == 1 && @{$u->[0]} == 2)
{
    croak("Poisson: list ref output bad for scalar input\n");
}

# Test the list ref input
$k = [5, 15];
($u, $l) = poissonLimits($k);
unless (@$u == 2 && @$l == 2)
{
    croak("Poisson: ref input produced insufficient output\n");
}

# Test the first calling sequence for list ref input
$u = poissonLimits($k);
unless (@$u == 2 && @{$u->[0]} == 2)
{
    croak("Poisson: list ref output bad for list ref input\n");
}

##########################################################################

# Test PDL input/output
eval{use PDL;};
unless ($@)
{
    # Test the second calling sequence with piddles
    ($nsuccess, $ntotal) = (pdl(4,8), pdl(12,30));
    ($u, $l) = binomialLimits($nsuccess, $ntotal);
    unless ($u->nelem == $l->nelem)
    {
    	croak("Binomial: PDL input generates bad pdl output\n");
    }

    # Test the first calling sequence with piddles
    $u = binomialLimits($nsuccess, $ntotal);
    unless ($u->getdim(0) == 2 && $u->getdim(1) == 2)
    {
        croak("Binomial: PDL input generates bad 2D pdl output\n");
    }
    
    # Test the second calling sequence with piddles
    $k = pdl(4,8);
    ($u, $l) = poissonLimits($k);
    unless ($u->nelem == $l->nelem)
    {
        croak("Poisson: PDL input generates bad pdl output\n");
    }

    # Test the first calling sequence with piddles
    $u = poissonLimits($k);
    unless ($u->getdim(0) == 2 && $u->getdim(1) == 2)
    {
        croak("Poisson: PDL input generates bad 2D pdl output\n");
    }
}

print "OK\n";
