using LinearAlgebra

# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams

"""
    poldif(x, alpha, B)

Computes the differentiation matrices D1, D2, ..., DM on arbitrary nodes, x, with associated weights, alpha.

# Arguments
- x: vector of N distinct nodes.
- alpha: vector of weight values alpha(x), evaluated at x = x[k].
- B: matrix of size M x N;  where M is the highest derivative required.  It should contain the quantities 
B[ell,j] = beta(ell,j) = (ell-th derivative of alpha(x[j]))/alpha(x[j]). Note that M must satisfy, 0 < M < N-1.

# Outputs
- DM: DM[1:N,1:N,ell] contains ell-th derivative matrix, ell=1,...,M.
"""
function poldif(x, alpha, B)


    N = length(x)                      
    x = x[:]                     # Make sure x is a column vector 
                              
    alpha = alpha[:]             # Make sure alpha is a column vector.
    M = size(B,1)                # First dimension of B is the number of derivative matrices to be computed.
        
    XX = repeat(x,1,N)
    DX = XX - XX'                  # DX contains entries x[k]-x[j].
        
    DX[diagind(DX)] .= 1           # Put 1's one the main diagonal.
        
    c = alpha .* prod(DX,dims=2)   # Quantities c[j].
        
    C = repeat(c,1,N)
    C = C ./ C'                    # Matrix with entries c[k]/c[j].
        
    Z = 1 ./ DX                    # Z contains entries 1/(x[k]-x[j]) with zeros on the diagonal.
    Z[diagind(Z)] = zeros(N,1)              
        
    X = Z'                         # X is same as Z', but with diagonal entries removed.
    X = reshape([X[l] for l in CartesianIndices(X) if l[1] ≠ l[2]],N-1,N)
        
    Y = ones(N-1,N)                     # Y is matrix of cumulative sums
    D = Matrix{Float64}(I, (N,N))       # D differentiation matrices.

    DM = Array{Float64}(undef,N,N,M)
    for ell = 1:M
        Y = cumsum([B[ell,:]'; ell*Y[1:N-1,:] .* X], dims=1)     # Diagonals
        D = ell*Z .* (C .* repeat(diag(D),1,N) - D)              # Off-diagonals
        D[diagind(D)] = Y[N,:]                                   # Correct the diagonal
        DM[:,:,ell] = D                                          # Store the current D
    end

    return DM
end

"""
    poldif(x, M)
    
Computes the differentiation matrices D1, D2, ..., DM on arbitrary nodes, x, with unit weights.

# Arguments
- x: vector of N distinct nodes.
- M: the number of derivatives required [integer]. Note that M must satisfy, 0 < M < N-1.

# Outputs
- DM: DM[1:N,1:N,ell] contains ell-th derivative matrix, ell=1,...,M.
"""
function poldif(x, M)

    N = length(x)                      

    alpha = ones(N)              
    B = zeros(Float64,M,N)

    return poldif(x, alpha, B)
end