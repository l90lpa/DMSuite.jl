using LinearAlgebra
using ToeplitzMatrices

function sincdif(N, M, h)

#  The function [x, DM] = sincdif(N, M, h) computes sinc the
#  differentiation matrices D1; D2; ...; DM on equidistant points.
#
#  Input:
#  N:    Number of points; i.e.; order of differentiation matrix.
#  M:    Number of derivatives required [integer].
#  h:    Step-size (real, positive).
#
#  Note:  0 < M .< N-1.
#
#  Output:
#  x:    Vector of nodes.
#  DM:   DM[1:N,1:N,l] contains l-th derivative matrix, l=1..M.
#
#  J.A.C. Weideman; S.C. Reddy 1998.  Help lines corrected 
#  by JACW; March/April 2003.
    
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