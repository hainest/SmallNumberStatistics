; bdtri.pro, July 21, 2012
; Copyright (c) 2012, Tim Haines
;
;   This program is free software: you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation, either version 3 of the License, or
;   (at your option) any later version.
;
;   This program is distributed in the hope that it will be useful,
;   but WITHOUT ANY WARRANTY; without even the implied warranty of
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;   GNU General Public License for more details.
;
;   You should have received a copy of the GNU General Public License
;   along with this program.  If not, see http://www.gnu.org/licenses/.
;
; NAME:
;       BDTRI
;
; AUTHOR:
;       Tim Haines, University of Missouri - Kansas City
;       5110 Rockhill Road, Kansas City, MO 64110, USA
;       tdhq36@mail.umkc.edu
;
; PURPOSE:
;       This function computes the inverse binomial distribution.
;
; CATEGORY:
;       Statistics
;
; CALLING SEQUENCE:
;       p = BDTRI(K, N, Y)
;
; INPUTS:
;       K:      A strictly positive integer that specifies the number
;               of successes in N Bernoulli trials.
;       
;       N:      An integer strictly greater than K that specified the
;               number of Bernoulli trials.
;        
;       Y:      The cumulative probability of success.
;
; RETURNS:
;       p:      The event probability of success.
;
; EXAMPLE:
;       What is probability of a single event given one success in
;       5 trials and a cumulative probability of 0.4018776?
;         p_upper = bdtri(1,5,0.4018776)
;                 = 0.26646
;
;         p_lower = 1 - bdtri(4,5,0.4018776)
;                 = 0.16666
;
;       Note that this is the inverse problem to asking what the
;       cumulative probability of getting one six in five die rolls is.
;       Answer: 0.4018776. However, we can see that the inverse
;       problem does not produce a unique solution, but a range with
;       the "correct" answer being given by the lower bound
;       (1/6 = 0.16666...)
;
; REFERENCES:
;       This implementation is based on the method used in the Cephes
;       math library developed by Stephen L. Moshier.
;       http://www.netlib.org/cephes/
;
; ALGORITHM:
;       Find the event probability p such that the sum of the terms
;       0 through K of the Binomial probability density is equal to
;       the given cumulative probability Y using the inverse of the
;       incomplete beta integral (incbi) and the relation
;
;           1 - p = incbi( N-K, K+1, Y )
;
;       incbi(a,b,y) is found by solving for the probability x such that
;           incbet(a,b,x) - y = 0
;       by way of a simple bisection root-finding method.
;
; DEPENDENCIES:
;       IBETA       - IDL >4.0
;		BISECTION   - IDL >4.0
;                     Copyright (c) Erik Rosolowsky <eros@cosmic>
;                     https://people.ok.ubc.ca/erosolo/idl/lib/bisection.pro

function bdtri, k, n, y
    ; check the domain
    if k lt 0 then begin
        message, strjoin('bdtri domain error: k must be strictly positive')
    endif
    
    if n le k then begin
        message, strjoin('bdtri domain error: n must be strictly greater than k')
    endif
    
    if y lt 0.0 or y gt 1.0 then begin
        message, strjoin('bdtri domain error: y must be between 0 and 1, inclusive')
    endif
    
	common parms, a, b, p
	a = k
	b = n
	p = y

	return, 1 - bisection(0.5, 'incbi', itmax=200, radius=0.5, tol=1e-7, /double)
end

function incbi, x
	common parms
	return, ibeta(b-a, a+1, x[0]) - p
end
