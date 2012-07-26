; Test.pro, July 21, 2012
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
;       TEST
;
; AUTHOR:
;       Tim Haines, tdhq36@mail.umkc.edu
;
; PURPOSE:
;       This function tests binomial_limits and poisson_limits for
;       working functionality.
;
; CATEGORY:
;       Statistics and probability
;
; CALLING SEQUENCE:
;       TEST
;
; RETURNS:
;       The string 'OK' if all passed, else an error message.
;
; DEPENDENCIES:
;       BINOMIAL_LIMITS
;       POISSON_LIMITS

pro test

    nsuccess = 9
    ntotal = 109
    
    ; make sure we are getting close to Gehrels
    p = binomial_limits(nsuccess, ntotal)
    if abs(p[0,0]-0.118) gt 1e-2 or abs(p[0,1] - 0.0561) gt 1e-2 then begin
            message, 'Binomial: Not close to Gehrels!'
    endif
    
    ; make sure default cl is working
    p2 = binomial_limits(nsuccess, ntotal, 0.841344746068543)
    if abs(p2[0,0]-p[0,0]) gt 1e-3 or abs(p2[0,1] - p[0,1]) gt 1e-3 then begin
        message, 'Binomial: The default cl is not working!'
    endif
    
    ; make sure the 1-sigma cl works
    p2 = binomial_limits(nsuccess, ntotal, 1.0, /sigma)
    if abs(p2[0,0]-p[0,0]) gt 1e-3 or abs(p2[0,1]-p[0,1]) gt 1e-3 then begin
        message, 'Binomial: Not sigma is not working'
    endif
    
    ; Test an extremum from Gehrels
    nsuccess = 0
    ntotal = 100
    p = binomial_limits(nsuccess, ntotal, 3.0, /sigma)
    if abs(p[0,0]-0.0639) gt 1e-3 or abs(p[0,1]) ne 0 then begin
        message, 'Binomial: Not close to Gehrels - 2!'
    endif
    
    k = 13
    
    ; make sure we are getting close to Gehrels
    p = poisson_limits(k)
    if abs(p[0,0]-17.70) gt 1e-2 or abs(p[0,1]-9.441) gt 1e-2 then begin
        message, 'Poisson: Not close to Gehrels!'
    endif
    
    ; make sure default cl is working
    p2 = poisson_limits(k, 0.841344746068543)
    if abs(p2[0,0]-p[0,0]) gt 1e-3 or abs(p2[0,1]-p[0,1]) gt 1e-3 then begin
        message, 'Poisson: The default cl is not working!'
    endif
    
    ; make sure the 1-sigma cl works
    p2 = poisson_limits(k, 1.0, /sigma)
    if abs(p2[0,0]-p[0,0]) gt 1e-3 or abs(p2[0,1]-p[0,1]) gt 1e-3 then begin
        message, 'Poisson: Not sigma is not working!'
    endif
    
    ; Test an extremum from Gehrels
    k = 100
    p = poisson_limits(k, 3.0, /sigma)
    if abs(p[0,0]-133.8) gt 0.1 or abs(p[0,1]-72.65) gt 0.1 then begin
        message, 'Poisson: Not close to Gehrels - 2!'
    endif
   
   print, 'OK' 
end