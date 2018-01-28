; pdtri.pro, July 21, 2012
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
;       PDTRI
;
; AUTHOR:
;       Tim Haines, tdhq36@mail.umkc.edu
;
; PURPOSE:
;       This function computes the inverse Poisson distribution.
;
; CATEGORY:
;       Statistics and probability
;
; CALLING SEQUENCE:
;       M = PDTRI(K, Y)
;
; INPUTS:
;       K:      A strictly nonnegative integer that specifies the number
;               of events per time interval.
;        
;       Y:      The probability.
;
; RETURNS:
;       P:      The Poisson parameter (alias lambda)
;
; REFERENCES:
;       This implementation is based on the method used in the Cephes
;       math library developed by Stephen L. Moshier.
;       http://www.netlib.org/cephes/
;
; ALGORITHM:
;       Finds the Poisson variable M such that the integral from 0 to
;       M of the Poisson density is equal to the given probability Y
;       using the inverse of the complemented incomplete gamma integral
;       and the relation
;
;           M = igami( K+1, Y )
;
;       igami(a,p) is found by solving for the probability x such that
;
;           igamc( a, x ) = p
;
;       by way of a simple bisection root-finding method. igamc is the
;       complemented incomplete gamma integral given by
;
;           igamc(a, x) = 1 - igam(a, x)
;
;       The starting value, x0, for the bisection search is given by
;       the normal approximation to the Poisson distribution.
;
; DEPENDENCIES:
;       IGAMMA      - IDL >4.0
;       BISECTION   - IDL >4.0

function igami, x
    compile_opt idl2, hidden
    common igami_parms, a, p
    
	return, 1 - igamma(a, x, /double) - p
end

function pdtri, k, y
    ; check the domain
    if k lt 0 then begin
        message, 'pdtri domain error: k must be strictly nonnegative'
    endif

    if floor(k) lt k then begin
		message, 'pdtri domain error: k must be an integer'
    endif
    
    if y lt 0.0 or y ge 1.0 then begin
        message, 'pdtri domain error: y must be in the interval [0, 1)'
    endif
    
	common igami_parms
	a = k + 1
	p = y
    
    ; Set the starting point according to the normal approximation
    d = 1 / (9.0D * double(a))
    t = 1 - d - gauss_pdf(p) * sqrt(d)
    x0 = a * (t^3)

    return, bisection(x0, 'igami', itmax=200, radius=x0, tol=1e-7)
end
