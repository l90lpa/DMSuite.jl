using LinearAlgebra
using ToeplitzMatrices

# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams

"""
    sincdif(N, M, h)
    
Computes sinc the differentiation matrices D1, D2, ..., DM on equidistant points.

# Arguments
- N: number of points; i.e.; order of differentiation matrix.
- M: number of derivatives required [integer]. Note that M must satisfy, 0 < M < N-1.
- h: step-size (real, positive).

# Outputs
- x:  vector of nodes.
- DM: DM[1:N,1:N,ell] contains ell-th derivative matrix, ell=1,...,M.
"""
function sincdif(N, M, h)
    
    k = 1:N-1
    t = k*pi
    x = Vector(h*(-(N-1)/2:(N-1)/2))
    
    sigma = zeros(N-1)
    
    DM = Array{Float64}(undef,N,N,M)
    for ell = 1:M
        sigma = (-ell*sigma .+ imag.(exp.(im*t)* im^ell)) ./t
        col = (pi/h)^ell*[imag(im^(ell+1))/(ell+1); sigma]
        row = (-1.0)^ell*col
        row[1] = col[1]
        DM[:,:,ell] = Matrix(Toeplitz(col,row))
    end
    
    return x, DM
end