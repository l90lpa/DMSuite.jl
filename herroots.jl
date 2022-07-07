using LinearAlgebra

function herroots(N)

#  The function r = herroots(N) computes the roots of the 
#  Hermite polynomial of degree N.
#
#  J.A.C. Weideman; S.C. Reddy 1998.
 
    J = SymTridiagonal(zeros(N), sqrt.(1:N-1)) # Jacobi matrix
    r = eigvals(J) ./ sqrt(2)                  # Compute eigenvalues

    return r
end