using LinearAlgebra

# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams
"""
    lagroots(N)
    
Computes the roots of the Laguerre polynomial of degree N.
"""
function lagroots(N);

    
    J = SymTridiagonal(Vector(1:2:2*N-1), -Vector(1:N-1)) # Jacobi matrix
    r = eigvals(J)                                        # Compute eigenvalues
    
    return r
end 
    