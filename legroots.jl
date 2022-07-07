function legroots(N);

#  The function r = legroots(N) computes the roots of the 
#  Legendre polynomial of degree N.
#
#  J.A.C. Weideman, S.C. Reddy 1998.
    
    n = 1:N-1                       #  Indices
    d = n ./ sqrt.(4*n .^ 2 .- 1)    #  Create subdiagonals
    J = SymTridiagonal(zeros(N), d) # Jacobi matrix
    r = eigvals(J)                  #  Compute eigenvalues

    return r
end