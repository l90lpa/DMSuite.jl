using LinearAlgebra
using AbstractFFTs

# Originally implemented in Matlab by S.C. Reddy & J.A.C. Weideman, implemented in Julia by L.P. Adams
"""
    fourdifft(f,M)

Computes the M-th derivative of the function f[x] using the Fourier differentiation process. 
The function is assumed to be 2pi-periodic & the input data values f should 
correspond to samples of the function at N equispaced points on [0, 2pi).
The Fast Fourier Transform is used.

# Arguments
- f: vector of samples of f[x] at x = 0, 2pi/N, 4pi/N, ... , (N-1)2pi/N
- M: derivative required [non-negative integer]

# Outputs
- Dmf:     M-th derivative of f
"""
function fourdifft(f,M)
    
    
    f = f[:]                       # Make sure f is a column vector
    N = length(f)
        
    N1 = Int64(floor((N-1)/2))           # Set up wavenumbers           
    wave = iseven(N) ? im*[(0:N1);  (-N/2)*rem(M+1,2); (-N1:-1)] : im*[(0:N1); (-N1:-1)]
    
    Dmf = ifft(((wave) .^ M) .* fft(f))   # Transform to Fourier space, take deriv & return to physical space.
                                         
    if maximum(abs.(imag.(f))) == 0           # Real data in real derivative out
        Dmf = real(Dmf)                
    end 
    return Dmf
end