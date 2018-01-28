; poisson_limits.pro, July 21, 2012
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
;       POISSON_LIMITS
;
; AUTHOR:
;       Tim Haines, thaines.astro@gmail.com
;
; PURPOSE:
;       This function computes the single-sided upper and lower
;       confidence limits for the Poisson distribution.
;
; CATEGORY:
;       Statistics and probability
;
; CALLING SEQUENCE:
;       P = POISSON_LIMITS(K, [CL [,/SIGMA]])
;
; INPUTS:
;       K:      A strictly nonnegative integer that specifies the
;               number of observed events. Can be an array.
;        
; OPTIONAL INPUTS:
;       CL:     The confidence level in the interval [0, 1). The default
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
;       To compute the confidence limits of seeing 20 events in 8
;       seconds at the 2.5 sigma, use
;
;           p = poisson_limits(20, 2.5, /sigma)
;               p = [ [34.1875], [10.5711] ]
;
;       However, recall that the Poisson parameter is defined as the
;       average rate, so it is necessary to divide these values by
;       the time (or space) interval over which they were observed.
;
;       Since these are the confidence limits, the fraction would be
;       reported as
;
;           2.5 (+4.273, -1.321) observations per second
;
; DEPENDENCIES:
;       PDTRI   - IDL >5.3

function poisson_limits, k, cl, sigma=sigma

    ; Set CL = 1 sigma by default
    if ~keyword_set(cl) then begin
        cl = 1.0
        sigma = 1
    endif
    
    if keyword_set(sigma) then begin
        cl = gaussint(cl)
    endif
   
    p = make_array(n_elements(k), 2, /double)
    
    for i=0, n_elements(k)-1 do begin
        p[i,0] = pdtri(k[i], 1 - cl)
        
        ; See Gehrels (1986) for details
        if k[i] eq 0 then begin
            p[i,1] = 0.0D
        endif else begin
            p[i,1] = pdtri(k[i] - 1, cl)
        endelse
    endfor

    return, p
end
