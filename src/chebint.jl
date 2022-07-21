
# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams
"""
    chebint(fk, x)

Evaluates the polynomial interpolant of the data points (xk[j], fk[j]), where 
xk[j] are the Chebyshev nodes, at the points x. 
Requires two or more data points.

# Arguments 
- fk: vector of y-coordinates of data; at Chebyshev points xk[j] = cos((j-1)*pi/(N-1)), j = 1...N.
- x: vector of x-values where polynomial interpolant is to be evaluated.

# Outputs
- p: vector of interpolated values.

# Details
The code implements the barycentric formula; see page 252 in
P. Henrici; Essentials of Numerical Analysis; Wiley; 1982.
"""
function chebint(fk, x)
    # Note: that if some fk > 1/eps, with eps the machine epsilon 
    # the value of eps in the code may have to be reduced.
    
    fk = fk[:]; x = x[:];                       # Make sure data are column vectors.
    
    N = length(fk); 
    M = length(x)
         
    xk = cheb1Extrema(N-1)
    
    w = Vector((-1.0) .^ (0:N-1));              # w = weights for Chebyshev formula
    w[1] = w[1]/2; w[N] = w[N]/2
     
    D = repeat(x,1,N) - repeat(xk,1,M)';        # Compute quantities x-x[k] & their reciprocals.
    D = 1 ./ (D + eps() * (D .== 0));                  
      
    p = D*(w .* fk) ./ (D*w);                   # Evaluate interpolant as matrix-vector products.
    
    return p
end