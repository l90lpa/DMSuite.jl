using LinearAlgebra
using ToeplitzMatrices

function chebdif(N, M)

#  The function [x, DM] =  chebdif(N,M) computes the differentiation 
#  matrices D1; D2; ...; DM on Chebyshev nodes. 
# 
#  Input:
#  N:        Size of differentiation matrix.        
#  M:        Number of derivatives required [integer].
#  Note:     0 .< M <= N-1.
#
#  Output:
#  x:        The Chebyshev nodes
#  DM:       DM[1:N,1:N,ell] contains ell-th derivative matrix, ell=1..M.
#
#  The code implements two strategies for enhanced 
#  accuracy suggested by W. Don & S. Solomonoff in 
#  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 [1994].
#  The two strategies are [a] the use of trigonometric 
#  identities to avoid the computation of differences 
#  x[k]-x[j] & [b] the use of the "flipping trick"
#  which is necessary since sin t can be computed to high
#  relative precision when t is small whereas sin (pi-t) cannot.
#  Note added May 2003:  It may; in fact; be slightly better not to
#  implement the strategies [a] & [b].   Please consult the following
#  paper for details:   "Spectral Differencing with a Twist"; by
#  R. Baltensperger & M.R. Trummer; to appear in SIAM J. Sci. Comp. 

#  J.A.C. Weideman; S.C. Reddy 1998.  Help notes modified by 
#  JACW; May 2003.
    
    n1 = Int64(floor(N/2)); n2  = Int64(ceil(N/2));  # Indices used for flipping trick.
    
    k = (0:N-1)                                      # Compute theta vector.
    th = k*pi/(N-1)
    
    x = (sin.(pi * (N-1:-2:1-N)' / (2*(N-1))))';        # Compute Chebyshev points.
    
    T = repeat(th/2,1,N);                
    DX = 2 * sin.(T'+T) .* sin.(T'-T);               # Trigonometric identity. 
    DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:], dims = 2), dims = 1)];   # Flipping trick. 
    DX[diagind(DX)] .= 1.0                           # Put 1's on the main diagonal of DX.
    
    
    C = Matrix(SymmetricToeplitz((-1.0).^k));        # C is the matrix with 
    C[1,:] = C[1,:]*2; C[N,:] = C[N,:]*2;            # entries c[k]/c[j]
    C[:,1] = C[:,1]/2; C[:,N] = C[:,N]/2
    
    Z = 1 ./ DX;                                     # Z contains entries 1/(x[k]-x[j])  
    Z[diagind(Z)] .= 0.0;                            # with zeros on the diagonal.

    D = Matrix{Float64}(I, (N,N));                   # D contains diff. matrices.
    
    DM = Array{Float64}(undef,N,N,M)
    for ell = 1:M
        D = ell*Z.*(C.*repeat(diag(D),1,N) - D);     # Off-diagonals
        D[diagind(D)] = -sum(D, dims=2)              # Correct main diagonal of D
        DM[:,:,ell] = D;                             # Store current D in DM
    end
    
    return x, DM
end