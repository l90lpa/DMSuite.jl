
# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams
"""
    herdif(N, M, b)

Computes the differentiation matrices D1, D2, ..., DM on Hermite points.

# Arguments
- N: number of points, i.e., order of differentiation matrices [integer].
- M: number of derivatives required [integer]. Note that M must satisfy, 0 < M < N-1.
- b: scaling parameter [real, positive].

# Outputs
- x:    vector of nodes (zeros of Hermite polynomial of degree N scaled by the parameter b.)
- DM:   DM[1:N,1:N,ell] contains ell-th derivative matrix, ell=1,..,M.
"""
function herdif(N, M, b)

    
    x = herroots(N)                      # Compute Hermite roots.
    
    alpha = exp.(-x .^ 2/2)              # Compute weights.
    
    beta = Matrix{Float64}(undef, M+1, N)
    beta[1,:] .= 1                        # Set up beta matrix s.t. beta(l,j) =
    beta[2,:] = -x                        # (l-th derivative of alpha(x))/alpha(x)
                                          # evaluated at x = x[j].
    for ell = 3:M+1                         
        beta[ell,:] = -x' .* beta[ell-1,:]' - (ell-2)*beta[ell-2,:]'
    end
    
    beta = beta[2:end,:]                  # Remove initializing row from beta()
    
    DM = poldif(x, alpha, beta)           # Compute differentiation matrix [b=1].
    
    x = x ./ b                            # Scale nodes by the factor b.
    
    for ell = 1:M                         # Adjust for b not equal to 1.
        DM[:,:,ell] .*= (b^ell)
    end
    
    return x, DM
end