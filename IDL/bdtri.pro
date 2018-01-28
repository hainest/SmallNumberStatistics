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
;       Tim Haines, thaines.astro@gmail.com
;
; PURPOSE:
;       This function computes the inverse binomial distribution.
;
; CATEGORY:
;       Statistics and probability
;
; CALLING SEQUENCE:
;       P = BDTRI(K, N, Y)
;
; INPUTS:
;       K:      A strictly nonnegative integer that specifies the number
;               of successes in N Bernoulli trials.
;       
;       N:      An integer strictly greater than K that specifies the
;               number of Bernoulli trials.
;        
;       Y:      The cumulative probability of success.
;
; RETURNS:
;       P:      The event probability of success.
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
;           1 - P = incbi( N-K, K+1, Y ).
;
;       incbi(a,b,y) is found by solving for the probability x such that
;
;           incbet(a,b,x) - y = 0
;
;       by way of a simple bisection root-finding method.
;
; DEPENDENCIES:
;       IBETA       - IDL >4.0
;       BISECTION   - IDL >4.0

function incbi, x
    compile_opt idl2, hidden
    common incbi_parms, a, b, p
    
	return, ibeta(b-a, a+1, x) - p
end

function bdtri, k, n, y
    compile_opt idl2
    
	; check the domain
    if k lt 0 then begin
        message, 'bdtri domain error: k must be strictly nonnegative'
    endif
    
    if floor(k) lt k then begin
		message, 'bdtri domain error: k must be an integer'
    endif
    
    if n le k then begin
        message, 'bdtri domain error: n must be strictly greater than k'
    endif
    
    if floor(n) lt n then begin
		message, 'bdtri domain error: n must be an integer'
	endif
    
    if y lt 0.0 or y gt 1.0 then begin
        message, 'bdtri domain error: y must be in the interval [0,1]'
    endif
    
	common incbi_parms
	a = k
	b = n
	p = y

    return, 1 - bisection(0.5, 'incbi', itmax=200, radius=0.5, tol=1e-7)
end
