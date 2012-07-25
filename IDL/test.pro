function test, a, b
    if ~keyword_set(a) then print, 'a is not set!'
    a = 1
    if ~keyword_set(a) then print, 'a is not set!'
    
    if ~keyword_set(b) then print, 'b is not set!'
    b = 1
    if ~keyword_set(b) then print, 'b is not set!'
end
