; binomial_limits.pro, July 21, 2012
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
;       BINOMIAL_LIMITS
;
; AUTHOR:
;       Tim Haines, tdhq36@mail.umkc.edu
;
; PURPOSE:
;       This function computes the single-sided upper and lower
;		confidence limits for the binomial distribution.
;
; CATEGORY:
;       Statistics and probability
;
; CALLING SEQUENCE:
;       P = BINOMIAL_LIMITS(NSUCCESS, NTOTAL, [CL [,/SIGMA]])
;
; INPUTS:
;       NSUCCESS:	A strictly nonnegative integer that specifies the
;               	number of successes in NTOTAL Bernoulli trials. Can
;                   be an array.
;       
;       NTOTAL:		An integer strictly greater than NSUCCESS that
;               	specifies the number of Bernoulli trials. Can be
;                   an array.
;        
; OPTIONAL INPUTS:
;       CL:     The confidence level in the interval [0, 1]. The default
;               is 0.8413 (i.e., 1 sigma)
;
;       SIGMA:  If this keyword is set, then CL is assumed to be a
;               multiple of sigma, and the actual confidence level
;               is computed from the standard normal distribution with
;               parameter CL.
;
; RETURNS:
;       A two-dimensional array such that the array P[*,0] contains the upper
;       single-sided confidence limits, and the array P[*,1] contains the lower
;       single-sided confidence limits.
;
; REFERENCES:
;       N. Gehrels. Confidence limits for small numbers of events in astrophysical
;       data. The Astrophysical Journal, 303:336–346, April 1986.
;
; EXAMPLE:
;       I have a mass bin with 100 galaxies (20 reds and 80 blues)
;       and I am computing the fraction of reds to blues, then for this
;       bin NSUCCESS = 20 and NTOTAL = 100.
;
;       To compute the confidence limits at the 2.5 sigma level, use
;
;           p = binomial_limits(20, 100, 2.5, /sigma)
;               p = [ [0.31756], [0.11056] ]
;
;       Since these are the confidence limits, the fraction would be
;       reported as
;
;           0.2 (+0.11756, -0.08944)
;
; DEPENDENCIES:
;       BDTRI   - IDL >5.3

function binomial_limits, nsuccess, ntotal, cl, sigma=sigma
	; Resolve the 'incbi' routine in bdtri.pro so that 'bisection'
	; can find it with call_function
	resolve_routine, 'bdtri', /is_function
	
	; Set CL = 1 sigma by default
    if ~keyword_set(cl) then begin
        cl = 1.0D
        sigma = 1
    endif
    
    if keyword_set(sigma) then begin
        cl = gaussint(cl)
    endif
    
    if n_elements(nsuccess) ne n_elements(ntotal) then begin
        message, 'nsuccess and total must have same size'
    endif
    
    p = make_array(n_elements(nsuccess), 2, /double)
    
    for i=0, n_elements(nsuccess)-1 do begin
        nfail = ntotal[i] - nsuccess[i];
        
        ; See Gehrels (1986) for details
        if nfail eq 0 then begin
            p[i,0] = 1.0D
        endif else begin
            p[i,0] = bdtri(nsuccess[i], ntotal[i], 1 - cl)
        endelse
        
        ; See Gehrels (1986) for details        
        if nsuccess[i] eq 0 then begin
            p[i,1] = 0.0D
        endif else begin
            p[i,1] = 1 - bdtri(nfail, ntotal[i], 1 - cl)
        endelse
    endfor

	return, p
end
