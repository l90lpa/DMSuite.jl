

"""
    chebdifft(f,M)

Computes the M'th approximate Chebyshev derivatives of the data vector f, using the FFT.

# Arguments
- f: vector of length N containing function values at the Chebyshev points x[k] = cos((k-1)*pi/(N-1)), k = 1...N.
- M: derivative required [positive integer]
    
# Outputs
- Dmf: vector containing approximate M'th derivative

# Details
A Fast Fourier Transform is used compute the Chebyshev cofficients
of the data vector. A recursion formula is used to compute the
Chebyshev coefficients for each derivative. A FFT is then used again
to compute the derivatives in physical space.
"""
function chebdifft(f,M)

    f=f[:];                                      # Make sure f is a vector     
    N=length(f);      
    a0=fft([f; reverse(f[2:N-1], dims = 1)]);    # Extend & compute fft()
    a0=a0[1:N] .* [0.5; ones(N-2,1); 0.5] ./ (N-1);    # a0 contains Chebyshev coefficients of f
    
    a=[a0 zeros(N,M)];                             # Recursion formula
    for ell=1:M                                    # for computing coefficients
        a[N-ell,ell+1]=2*(N-ell)*a[N-ell+1,ell];   # of ell'th derivative 
        for k=N-ell-2:-1:1
            a[k+1,ell+1]=a[k+3,ell+1]+2*(k+1)*a[k+2,ell]
        end
        a[1,ell+1]=a[2,ell]+a[3,ell+1]/2
    end
    
    back=[2*a[1,M+1]; a[2:N-1,M+1]; 2*a[N,M+1]; reverse(a[2:N-1,M+1], dims = 1)]
    Dmf=0.5*fft(back);                        # Transform back to
    Dmf=Dmf[1:N];                             # physical space
    
    if maximum(abs.(imag.(f))) == 0           # Real data in real derivative out
        Dmf = real(Dmf)                
    end    
    
    Dmf
end