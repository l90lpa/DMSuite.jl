
# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams
"""
    fourint(fk, x)

Evaluates the trigonometric interpolant of the data (xk[j], fk[j]), where
xk are equidistant nodes, at the points x. Requires two or more data points

# Arguments
- fk:  vector of y-coordinates of data, at equidistant points xk[j] = (j-1)*2*pi/N,  j = 1...N.
- x:   vector of x-values where interpolant is to be evaluated.

# Outputs
- t:    vector of interpolated values.

# Details
The code implements the barycentric formula; see page 46 in
P. Henrici; Applied & Computational Complex Analysis III; Wiley; 1986.
(Note that if some fk .> 1/eps, with eps the machine epsilon
the value of eps in the code may have to be reduced.)
"""
function fourint(fk, x)
    
    fk = fk[:]; x = x[:];           # Make sure data are column vectors.
        
    N = length(fk) 
    M = length(x)
        
    xk = (2*pi/N) * (0:N-1)         # Compute equidistant points
     
    w = (-1.0) .^ (0:N-1)           # Weights for trig interpolation
        
    x2 = x/2; xk2 = xk/2
     
    D = repeat(x2,1,N) - repeat(xk2,1,M)'   # Compute quantities x-x[k]
    
    if iseven(N)                  
        D = 1.0 ./ tan.(D + eps()*(D .== 0))       #  Formula for N even
    else
        D = 1.0 ./ sin.(D + eps()*(D .== 0))       #  Formula for N odd
    end
    
    t = D*(w .* fk) ./ (D*w)            # Evaluate interpolant as matrix-vector products.

    return t
end