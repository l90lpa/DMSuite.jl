using LinearAlgebra
using ToeplitzMatrices
using FFTW

function fourdif(N,m)
#
# The function [x, DM] = fourdif(N,m) computes the m'th derivative Fourier 
# spectral differentiation matrix on grid with N equispaced points in [0,2pi)
# 
#  Input:
#  N:        Size of differentiation matrix.
#  M:        Derivative required [non-negative integer]
#
#  Output:
#  x:        Equispaced points 0, 2pi/N, 4pi/N, ... , (N-1)2pi/N
#  DM:       m'th order differentiation matrix
#
# 
#  Explicit formulas are used to compute the matrices for m=1 & 2. 
#  A discrete Fouier approach is employed for m>2. The program 
#  computes the first column and first row & then uses the 
#  toeplitz command to create the matrix.

#  For m=1 & 2 the code implements a "flipping trick" to
#  improve accuracy suggested by W. Don & A. Solomonoff in 
#  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 [1994].
#  The flipping trick is necesary since sin t can be computed to high
#  relative precision when t is small whereas sin (pi-t) cannot.
#
#  S.C. Reddy; J.A.C. Weideman 1998.  Corrected for MATLAB R13 
#  by JACW; April 2003.
     
    
    x = Vector(2 * pi * (0:N-1) / N)          # gridpoints
    h = 2*pi/N                                # grid spacing
    zi = im
    kk = (1:N-1)
    n1 = Int64(floor((N-1)/2)); n2 = Int64(ceil((N-1)/2))

    if m == 0                                 # compute first column of zeroth derivative matrix; which is identity
        col1 = [1; zeros(N-1)]                
        row1 = col1                            
        
    elseif m == 1                             # compute first column of 1st derivative matrix
        if rem(N,2) == 0                       
            topc = cot.((1:n2) * h/2)
            col1 = [0; 0.5*((-1.0) .^ kk) .* [topc; -reverse(topc[1:n1], dims = 1)]]
        else
            topc = csc.((1:n2) * h/2)
            col1 = [0; 0.5*((-1.0) .^ kk) .* [topc; reverse(topc[1:n1], dims = 1)]]
        end
        row1 = -col1      # first row
        
    elseif m == 2                              # compute first column of 2nd derivative matrix
        if rem(N,2) == 0                         
            topc = csc.((1:n2) * h/2) .^ 2
            col1 = [-pi^2/3/h^2-1/6; -0.5*((-1.0) .^ kk) .* [topc; reverse(topc[1:n1], dims = 1)]]
        else
            topc = csc.((1:n2) * h/2) .* cot.((1:n2) * h/2)
            col1 = [-pi^2/3/h^2+1/12; -0.5*((-1.0).^kk) .* [topc; -reverse(topc[1:n1], dims = 1)]]
        end
        row1 = col1        # first row 
        
    else                                       # employ FFT to compute 1st column of matrix for m>2
        N1 = Int64(floor((N-1)/2))                    
        mwave = iseven(N) ? zi*[(0:N1); (-N/2)*rem(m+1,2); (-N1:-1)] : zi*[(0:N1); (-N1:-1)]
        col1 = vec(real(ifft((mwave .^ m) .* fft([1; zeros(N-1)]))))
        if rem(m,2) == 0
            row1 = col1;   # first row even derivative
        else
            col1 = [0; col1[2:N]]; 
            row1 = -col1;  # first row odd derivative
        end
    end
    DM = Matrix(Toeplitz(col1,row1));                   

    return x, DM
end