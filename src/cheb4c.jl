

"""
    cheb4c(N)

Computes the fourth derivative matrix on Chebyshev interior points, incorporating 
the clamped boundary conditions:

u[1]=u"[1]=u[-1]=u"[-1]=0.

# Arguments
- N:     N-2 = order of differentiation matrix. (The interpolant has degree N+1.)

# Outputs
- x:      interior Chebyshev points [vector of length N-2]
- D4:     fourth derivative matrix  [size (N-2)x(N-2)]

# Details
The code implements two strategies for enhanced 
accuracy suggested by W. Don & S. Solomonoff in 
SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
The two strategies are (a) the use of trigonometric 
identities to avoid the computation of differences 
x[k]-x[j] & (b) the use of the "flipping trick"
which is necessary since sin t can be computed to high
relative precision when t is small whereas sin(pi-t) cannot.
"""
function cheb4c(N)

    n1 = Int64(floor(N/2-1))               # n1, n2 are indices used 
    n2 = Int64(ceil(N/2-1))                # for the flipping trick.

    k = 1:N-2;                             # Compute theta vector.
    th = k*pi/(N-1);                 

    x = cheb1Extrema(N-1)[2:end-1] # interior Chebyshev points.

    s = [sin.(th[1:n1]); reverse(sin.(th[1:n2]))];   # s = sin(theta, dims = 1)
                               
    alpha = s .^ 4;                           # Compute weight function
    beta1 = -4*s .^ 2 .* x ./ alpha;          # & its derivatives.
    beta2 =  4*(3*x .^ 2 .- 1) ./ alpha;   
    beta3 = 24*x ./ alpha
    beta4 = 24  ./  alpha
    B = [beta1'; beta2'; beta3'; beta4']

    T = repeat(th/2,1,N-2);                
    DX = 2*sin.(T'+T) .* sin.(T'-T);     # Trigonometric identity 
    DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:], dims = 2), dims = 1)];   # Flipping trick. 
    DX[diagind(DX)] .= 1;                # Put 1's on the main diagonal of DX.

    ss = s .^ 2 .* (-1) .^ k;            # Compute the matrix with entries
    S = repeat(ss, 1, N-2)               # c[k]/c[j]
    C = S ./ S';                      

    Z = 1 ./ DX;                         # Z contains entries 1/(x[k]-x[j]).
    Z[diagind(Z)] .= 0;                  # with zeros on the diagonal.

    X = Z';                              # X is same as Z'; but with diagonal entries removed.
    X = reshape([X[l] for l in CartesianIndices(X) if l[1] ≠ l[2]],N-3,N-2)

    Y = ones(N-3,N-2);                   # Initialize Y & D vectors.
    D = Matrix{Float64}(I, (N-2,N-2));   # Y contains matrix of cumulative sums
                                         # D scaled differentiation matrices.
    DM = Array{Float64}(undef,N-2,N-2,4);

    for ell = 1:4
        Y = cumsum([B[ell,:]'; ell*Y[1:N-3,:] .* X], dims=1); # Recursion for diagonals
        D = ell*Z .* (C .* repeat(diag(D),1,N-2) - D);        # Off-diagonal
        D[diagind(D)] = Y[N-2,:];                             # Correct the diagonal
        DM[:,:,ell] = D;                                      # Store D in DM
    end

    D4 = DM[:,:,4];                  # Extract fourth derivative matrix

   return x, D4

end