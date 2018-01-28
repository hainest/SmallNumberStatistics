; bisection.pro, July 21, 2012
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
;       BISECTION
;
; AUTHOR:
;       Tim Haines, thaines.astro@gmail.com
;
; PURPOSE:
;       This function finds a root of a function bisecting a region
;       around the given starting value.
;
;       NOTE: Only positive roots are considered.
;
; CATEGORY:
;       Numerical analysis
;
; CALLING SEQUENCE:
;       X = BISECTION(X0, FNAME [, ITMAX=ITMAX, TOL=TOL, RADIUS=RADIUS])
;
; INPUTS:
;       X0:         The starting point to look for a root.
;
;       FNAME:      A string literal containing the name of the function
;                   for which the root is being found.
;
; OPTIONAL INPUTS:
;       ITMAX:      Maximum number of iterations. Default is 50.
;
;       TOL:        Relative error between candidate root values.
;                   Default is 1e-4.
;
;       RADIUS:     The search radius around X0. Default is 0.5.
;
; RETURNS:
;       X:      The root of the function. An error is thrown if no
;               root is found in ITMAX iterations.
;
; REFERENCES:
;       This is a heavily modified version of the bisection routine
;       written by Erik Rosolowsky (<eros@cosmic>). The original may
;       be found at https://people.ok.ubc.ca/erosolo/idl/lib/bisection.pro
;
; ALGORITHM:
;       Attempt to find a candidate root in the interval
;       [0, (2^itmax)*radius] defined by a change in sign of FUNC on
;       one such interval.
;
;       If a candidate root is found, refine the search by halving the
;       interval around the candidate root found above until the
;       error tolerance is reached.

function bisection, xinit, fname, itmax = itmax, tol = tol, radius = radius

    if not keyword_set(itmax) then itmax = 50
    if not keyword_set(tol) then tol = 1e-4
    if not keyword_set(radius) then radius = 0.5

    ; make sure the given values are of the proper data type
    xinit = double(xinit)
    radius = double(radius)
    tol = double(tol)
    itmax = fix(itmax)
    
    x = xinit
    iter = 0
    
    repeat begin
        x1 = x+radius
        x0 = x-radius
        
        ; clip lower x value to zero since negative roots are unallowed
        if x0 < 0.0 then x0 = 0
        
        y1 = call_function(fname, x1)
        y0 = call_function(fname, x0)
        
        iter = iter+1
        
        if iter gt itmax then begin
            msg = STRJOIN([ $
                        'Unable to locate a root inside a radius of', $
                        strtrim(string(iter*x0),2), $
                        ' of the starting point x0 = ', $
                        strtrim(string(x0),2) $
                    ], /single)
            message, msg
        endif
        
        radius = radius*2
    endrep until (y1*y0 le 0)
    
    toler = abs(x0-x1)
    
    ; Refine the search by halving the interval around the candidate
    ; root found above until the error tolerance is reached
    while toler gt tol do begin
        xold = x
        x = (x0+x1)*0.5
        y = call_function(fname, x)
        
        x0 = x0+(x-x0)*(y*y0 gt 0) 
        y0 = y0+(call_function(fname, x)-y0)*(y*y0 gt 0) 
        
        x1 = x1+(x-x1)*(y*y1 gt 0) 
        y1 = y1+(call_function(fname, x)-y1)*(y*y1 gt 0) 
        
        toler = abs(x0-x1)
        iter = iter+1
        
        if iter gt itmax then begin
            message, 'Exceeded maximum number of iterations.'
        endif
    endwhile

    return, x
end
