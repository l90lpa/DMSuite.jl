using LinearAlgebra

function lagroots(N);

#  The function r = lagroots(N) computes the roots of the 
#  Laguerre polynomial of degree N.
#
#  J.A.C. Weideman, S.C. Reddy 1998.
    
    J = SymTridiagonal(Vector(1:2:2*N-1), -Vector(1:N-1)) # Jacobi matrix
    r = eigvals(J)                                        # Compute eigenvalues
    
    return r
end 
    