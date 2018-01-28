#########################################################################
# SmallNumberStatistics.pm, July 21, 2012
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
#########################################################################

package SmallNumberStatistics;
require Exporter;

@ISA       = qw(Exporter);
@EXPORT    = qw(poissonLimits binomialLimits);

use warnings;
use strict;
use Math::Cephes qw(ndtr pdtri bdtri);
use Carp qw(croak);

#########################################################################
# NAME:
#       poissonLimits
#
# AUTHOR:
#       Tim Haines, thaines.astro@gmail.com
#
# PURPOSE:
#       This function computes the single-sided upper and lower
#       confidence limits for the Poisson distribution.
#
# CATEGORY:
#       Statistics and probability
#
# CALLING SEQUENCE:
#       $p = poissonLimits($k, [$cl [, \%OPTIONS]])
#
#       ($upper, $lower) = poissonLimits($k, [, $cl [, \%OPTIONS]])
#
# INPUTS:
#       k:      A strictly nonnegative integer that specifies the
#               number of observed events. Can be a scalar, a list
#               reference, or a PDL object.
#        
# OPTIONAL INPUTS:
#       cl:     The confidence level in the interval [0, 1). The default
#               is 0.8413 (i.e., 1 sigma)
#
# OPTIONS:
#       sigma:  If this is set, then $cl is assumed to be a
#               multiple of sigma, and the actual confidence level
#               is computed from the standard normal distribution with
#               parameter $cl.
#
# RETURNS:
#       For the first calling sequence:
#
#           input type           output type
#       -----------------   -----------------------
#         scalar                list reference
#         PDL object            2D PDL object
#         list reference        list reference
#
#       Where the output is given as "list reference," the actual output
#       is a reference to a list of list references such that each sublist
#       is a two-element list with the first element containing the
#       upper confidence limit and the second element containing the
#       lower confidence limit.
#
#       For the second calling sequence:
#
#           input type           output type
#       -----------------   -----------------------
#         scalar                scalar
#         PDL object            1D PDL object
#         list reference        list reference
#
#       Where the output is listed as "list reference," each output is a list
#       reference to a simple list.
#
# REFERENCES:
#       N. Gehrels. Confidence limits for small numbers of events in astrophysical
#       data. The Astrophysical Journal, 303:336–346, April 1986.
#
# EXAMPLE:
#       Compute the confidence limits of seeing 20 events in 8
#       seconds at the 2.5 sigma using the first calling convention
#       with scalar input.
#
#           $p = poissonLimits(20, 2.5, {'sigma' => 1})
#               $p->[0]->[0] = 34.1875
#               $p->[0]->[1] = 10.5711
#
#       However, recall that the Poisson parameter is defined as the
#       average rate, so it is necessary to divide these values by
#       the time (or space) interval over which they were observed.
#
#       Since these are the confidence limits, the fraction would be
#       reported as
#
#           2.5 (+4.273, -1.321) observations per second
#
# DEPENDENCIES:
#       Math::Cephes (available from CPAN)
#
# NOTE:
#       Because pdtri ultimately solves for the inverse of the incomplete
#       complementary gamma function and the Cephes math library only allows
#       p < 0.5 for that distribution, then it issues a warning about possible
#       loss of precision when p >= 0.5. For this most part, this can be ignored.
#       However, this can have an effect for very small cl (hence, 1-cl ~= 1).
#
##########################################################################

sub poissonLimits
{
	my $opt;
	$opt = pop if ref($_[-1]) eq 'HASH';
	
	croak('Usage:
       $p = poissonLimits($k, [$cl [, \%OPTIONS]])

       ($upper, $lower) = poissonLimits($k, [, $cl [, \%OPTIONS]])
	') unless @_ >= 1;
   
    my ($k, $cl, @upper, @lower, $isPDL, $isScalar);
    ($k, $cl) = @_;

    # use a 1-sigma confidence level by default
    unless(defined $cl)
    {
        $cl = 1.0;
        $opt->{'sigma'} = 1;
    }
    
    if ($opt && $opt->{'sigma'})
    {
    	$cl = ndtr($cl);
    }
    
    # check for a piddle
    if(ref $k eq 'PDL')
    {
    	use PDL::Lite;
    	$isPDL = 1;
    	$k = [$k->list];
    }
    
    # Assume a scalar if not a reference
    if (ref $k eq '')
    {
    	$isScalar = 1;
    	$k = [$k];
    }
    
    for (@$k)
    {
        push @upper, pdtri($_, 1 - $cl);
        
        # See Gehrels (1986) for details
        if($_ == 0)
        {
        	push @lower, 0.0;
        }
        else
        {
            push @lower, pdtri($_ - 1, $cl);
        }
    }

    if(wantarray)
    {
    	if ($isScalar)
    	{
    		return ($upper[0], $lower[0]);
    	}
    	
    	if($isPDL)
	    {
	        return (pdl(\@upper), pdl(\@lower));
	    }
    	
    	return (\@upper, \@lower);
    }
    else
    {
    	if($isPDL)
	    {
	        return pdl(\@upper, \@lower);
	    }
	    
    	my @result;
    	
    	# collapse them into a 2D array
        for my $i (0..@upper-1)
        {
            push @result, [$upper[$i], $lower[$i]];
        }
        return \@result;
    }
}

#########################################################################
# NAME:
#       binomialLimits
#
# AUTHOR:
#       Tim Haines, thaines.astro@gmail.com
#
# PURPOSE:
#       This function computes the single-sided upper and lower
#       confidence limits for the binomial distribution.
#
# CATEGORY:
#       Statistics and probability
#
# CALLING SEQUENCE:
#       $p = binomialLimits($nsuccess, $ntotal, [, $cl [, \%OPTIONS]])
#
#       ($upper, $lower) = binomialLimits($nsuccess, $ntotal, [, $cl [, \%OPTIONS]])
#
# INPUTS:
#       nsuccess:   A strictly nonnegative integer that specifies the
#                   number of successes in $ntotal Bernoulli trials.
#       
#       ntotal:     An integer strictly greater than $nsuccess that
#                   specifies the number of Bernoulli trials.
#        
# OPTIONAL INPUTS:
#       cl:     The confidence level in the interval [0, 1]. The default
#               is 0.8413 (i.e., 1 sigma)
#
# OPTIONS:
#       sigma:  If this is set, then $cl is assumed to be a
#               multiple of sigma, and the actual confidence level
#               is computed from the standard normal distribution with
#               parameter $cl.
#
# RETURNS:
#       For the first calling sequence:
#
#           input type           output type
#       -----------------   -----------------------
#         scalar                list reference
#         PDL object            2D PDL object
#         list reference        list reference
#
#       Where the output is given as "list reference," the actual output
#       is a reference to a list of list references such that each sublist
#       is a two-element list with the first element containing the
#       upper confidence limit and the second element containing the
#       lower confidence limit.
#
#       For the second calling sequence:
#
#           input type           output type
#       -----------------   -----------------------
#         scalar                scalar
#         PDL object            1D PDL object
#         list reference        list reference
#
#       Where the output is listed as "list reference," each output is a list
#       reference to a simple list.
#
#       NOTE: For both calling sequences, both $nsuccess and $ntotal must be
#             of the same input type and the same length.
#
# REFERENCES:
#       N. Gehrels. Confidence limits for small numbers of events in astrophysical
#       data. The Astrophysical Journal, 303:336–346, April 1986.
#
# EXAMPLE:
#       I have a mass bin with 100 galaxies (20 reds and 80 blues)
#       and I am computing the fraction of reds to blues, then for this
#       bin NSUCCESS = 20 and NTOTAL = 100.
#
#       To compute the confidence limits at the 2.5 sigma level, using
#       the first calling convention with scalar inputs
#
#           $p = binomialLimits(20, 100, 2.5, {'sigma' => 1})
#               p->[0]->[0] = 0.31756
#               p->[0]->[1] = 0.11056
#
#       Since these are the confidence limits, the fraction would be
#       reported as
#
#           0.2 (+0.11756, -0.08944)
#
##########################################################################

sub binomialLimits
{
	my $opt;
    $opt = pop if ref($_[-1]) eq 'HASH';
    
    croak('Usage:
       $p = binomialLimits($nsuccess, $ntotal, [, $cl [, \%OPTIONS]])

       ($upper, $lower) = binomialLimits($nsuccess, $ntotal, [, $cl [, \%OPTIONS]])
    ') unless @_ >= 2;
   
    my ($nsuccess, $ntotal, $nfail, $cl, @upper, @lower, $isPDL, $isScalar);
    ($nsuccess, $ntotal, $cl) = @_;

    # use a 1-sigma confidence level by default
    unless(defined $cl)
    {
        $cl = 1.0;
        $opt->{'sigma'} = 1;
    }
    
    if ($opt && $opt->{'sigma'})
    {
        $cl = ndtr($cl);
    }
    
    # Make sure that both arguments have the same data type
    if (ref $nsuccess ne ref $ntotal)
    {
    	croak('$nsuccess and $ntotal must be the same type.');
    }
    
    # check if it's a piddle
    if(ref $nsuccess eq 'PDL')
    {
        use PDL::Lite;
        $isPDL = 1;
        $nsuccess = [$nsuccess->list];
        $ntotal = [$ntotal->list];
    }
    
    # Assume it's a scalar if it's not a reference
    if (ref $nsuccess eq '')
    {
        $isScalar = 1;
        $nsuccess = [$nsuccess];
        $ntotal = [$ntotal];
    }
    
    # Ensure the two inputs are the same length
    unless (@$nsuccess == @$ntotal)
    {
    	croak('$nsuccess and $ntotal must be the same length!');
    }

    for my $i (0..@$nsuccess-1)
    {
    	$nfail = $ntotal->[$i] - $nsuccess->[$i];
    	
    	# See Gehrels (1986) for details
    	if ($nfail == 0)
    	{
    		push @upper, 1.0;
    	}
    	else
    	{
            push @upper, bdtri($nsuccess->[$i], $ntotal->[$i], 1 - $cl);
    	}
        
        # See Gehrels (1986) for details
        if($nsuccess->[$i] == 0)
        {
            push @lower, 0.0;
        }
        else
        {
            push @lower, 1 - bdtri($nfail, $ntotal->[$i], 1 - $cl);
        }
    }

    if(wantarray)
    {
        if ($isScalar)
        {
            return ($upper[0], $lower[0]);
        }
        
        if($isPDL)
        {
            return (pdl(\@upper), pdl(\@lower));
        }
        
        return (\@upper, \@lower);
    }
    else
    {
        if($isPDL)
        {
            return pdl(\@upper, \@lower);
        }
        
        my @result;
        
        # collapse them into a 2D array
        for my $i (0..@upper-1)
        {
            push @result, [$upper[$i], $lower[$i]];
        }
        return \@result;
    }
}

1;
