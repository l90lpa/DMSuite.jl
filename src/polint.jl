


"""
    polint(xk, fk, x, alpxk, alpx)
    
Computes the polynomial interpolant of the data [xk, fk] under the provided weights, alpxk and alpx.

# Arguments
- xk:    vector of x-coordinates of data [assumed distinct].
- fk:    vector of y-coordinates of data.
- x:     vector of x-values where polynomial interpolant is to be evaluated.
- alpxk: vector of weight values sampled at the points xk.
- alpx:  vector of weight values sampled at the points x.

# Outputs
- p:    vector of interpolated values.

# Details
The code implements the barycentric formula; see page 252 in
P. Henrici; Essentials of Numerical Analysis; Wiley; 1982.

Except for certain nice node distributions
polynomial interpolation of high-degree is an ill-conditioned
problem.  This code does not test for conditioning so use with
care.
"""
function polint(xk, fk, x, alpxk, alpx)
    # Note: that if some fk > 1/eps, with eps the machine epsilon
    # the value of eps in the code may have to be reduced
    
    fk = fk./alpxk
    
    x = x[:]                              # Make sure the data are column vectors
    xk = xk[:];  fk = fk[:]
    alpxk = alpxk[:]; alpx = alpx[:]

    N = length(xk) 
    M = length(x)

    D = repeat(xk,1,N) - repeat(xk,1,N)'  # Compute the weights w[k]
    D[diagind(D)] .= 1 
    w = 1 ./ prod(D,dims=1)'                       
 
    D = repeat(x,1,N) - repeat(xk,1,M)'   # Compute quantities x-x[k] & their reciprocals. 
    D = 1 ./ (D .+ (eps() * (D==0)))                 
  
    p = alpx .* (D*(w .* fk) ./ (D*w))    # Evaluate interpolant as matrix-vector products.
                                            

    p
end

"""
    polint(xk, fk, x)

Evaluates the polynomial interpolant of the data [xk, fk], at the points x. Requires two or more data points.

# Arguments
- xk: vector of x-coordinates of data [assumed distinct].
- fk: vector of y-coordinates of data.
- x:  vector of x-values where polynomial interpolant is to be evaluated.

# Outputs
- p:    Vector of interpolated values.

# Details
Except for certain nice node distributions polynomial interpolation of high-degree is an ill-conditioned
problem.  This code does not test for conditioning so use with care.
""" 
function polint(xk, fk, x)

    
    alpxk = ones(Float64, length(xk))
    alpx = ones(Float64, length(x))
    
    polint(xk,fk,x,alpxk,alpx)    
end