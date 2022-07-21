using LinearAlgebra

# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams
"""
    herroots(N)

Computes the roots of the Hermite polynomial of degree N.
"""
function herroots(N)

    J = SymTridiagonal(zeros(N), sqrt.(1:N-1)) # Jacobi matrix
    r = eigvals(J) ./ sqrt(2)                  # Compute eigenvalues

    return r
end