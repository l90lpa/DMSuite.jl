
"""
    cheb1Roots(N)

Computes the roots of the Chebyshev polynomial of the first kind of degree N. AKA, Chebyshev-Gauss (CG) Points, or 
Chebyshev points of the first kind.
"""
function cheb1Roots(N)
    
    x = cos.(pi .* (2 .* (1:N) .- 1) ./ (2*N))

    if iseven(N)
        x[div(N,2)+1:end] = -reverse(x[1:div(N,2)])
    else
        x[div(N-1,2)+2:end] = -reverse(x[1:div(N-1,2)])
        x[div(N-1,2)+1] = 0.0
    end

    x
end



"""
    cheb1Extrema(N)

Computes the extrama of the Chebyshev polynomial of the first kind of degree N. AKA, Chebyshev-Gauss-Lobatto (CGL) Points, or 
Chebyshev points of the second kind.
"""
function cheb1Extrema(N)
    N += 1
    sin.(pi * (N-1:-2:1-N) / (2*(N-1)))
end