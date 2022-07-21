

"""
    legroots(N)

Computes the roots of the Legendre polynomial of degree N.
"""
function legroots(N);
    
    n = 1:N-1                       #  Indices
    d = n ./ sqrt.(4*n .^ 2 .- 1)    #  Create subdiagonals
    J = SymTridiagonal(zeros(N), d) # Jacobi matrix
    r = eigvals(J)                  #  Compute eigenvalues

    r
end