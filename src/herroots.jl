

"""
    herroots(N)

Computes the roots of the Hermite polynomial of degree N.
"""
function herroots(N)

    J = SymTridiagonal(zeros(N), sqrt.(1:N-1)) # Jacobi matrix
    r = eigvals(J) ./ sqrt(2)                  # Compute eigenvalues

    return r
end