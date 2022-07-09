using AbstractFFTs

function sincdifft(f, M, h)

#  Dmf = sincdifft(f, M, h) computes the m-th derivative of the function
#  f[x] using the sinc differentiation process.   The function is() 
#  assumed to be defined on the entire real line & the input values
#  correspond to samples of the function at N equispaced points
#  symmetric with respect to the origin.
#
#  Input:
#  f:    Vector of samples of f[x] at h*[-(N-1)/2:(N-1)/2]
#  M:    Number of derivatives required [integer].
#  h:    Step-size (real, positive).
#
#  Note:  0 < M .< N-1.
#
#  Output:
#  Dmf:   m-th derivative of f
#
#  J.A.C. Weideman; S.C. Reddy 2000.   Corrected for complex data
#  by JACW; April 2003. 
    
    f = f[:]'                 # Ensure f is a row vector
    N = length(f)     
    t = pi*(1:N-1)'           
    
    sigma = zeros(1,N-1)
                                # Compute first column & row of diff matrix
    for l = 1:M
        sigma = (-l*sigma .+ imag.(exp.(im*t) .* im^l)) ./ t
    end
    col = (pi/h)^M*[imag(im^(M+1))/(M+1) sigma]
    row = (-1)^M*col 
    row[1] = col[1]
    
    # Imbed first row of Toeplitz matrix into bigger circulant matrix:
                                
    rowbig = [row zeros(1,2^nextpow(2,2*N)-2*N+1)  reverse(col[2:N]', dims = 2)]
    
    # Multiply circulant matrix with data vector by using FFT:
    
    NN = length(rowbig)
    e = NN*ifft(rowbig)           # Eigenvalues of circulant matrix.
    fhat = fft([f zeros(1,NN-N)]) # Take FFT of padded data vector
    Dmf = ifft(e.*fhat)           # multiply the result by e-values &
    Dmf = Dmf[1:N]               # take inverse FFT. 
    
    if maximum(abs.(imag.(f))) == 0           # Real data in real derivative out
        Dmf = real(Dmf)                
    end 
    return Dmf
end