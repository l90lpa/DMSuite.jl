
# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams
"""
    lagdif(N, M, b)
    
Computes the differentiation matrices D1, D2, ..., DM on Laguerre points.

# Arguments
- N: number of points, i.e., order of differentiation matrices [integer].
- M: number of derivatives required [integer]. Note that M must satisfy, 0 < M < N-1.
- b: scaling parameter [real, positive].

# Outputs
- x: vector of nodes (zeros of Laguerre polynomial of degree N-1 plus x = 0), all scaled by the parameter b.
- DM: DM[1:N,1:N,ell] contains ell-th derivative matrix, ell=1,...,M.
"""
function lagdif(N, M, b)
     
    x = lagroots(N-1)                    # Compute Laguerre roots()
    x = [0; x]                           # Add a node at x=0 to facilitate the implementation of BCs.
    
    alpha = exp.(-x ./ 2)                    # Compute weights.
    
    beta = Matrix{Float64}(undef, M, N)
    for ell = 1:M                               # Set up beta matrix s.t.
        beta[ell,:] .= (-0.5)^ell               # beta(ell,j) is (l^th derivative
    end                                         # alpha(x))/alpha(x)
                                                # evaluated at x = x[j].
    
    DM = poldif(x, alpha, beta)             # Compute differentiation matrix [b=1].
    
    x = x ./ b                              # Scale nodes by the factor b.
    
    for ell = 1:M                           # Adjust for b not equal to 1.
        DM[:,:,ell] .*= (b^ell)
    end
    
    return x, DM
end