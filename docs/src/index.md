# DMSuite.jl

Spectral differentiation and interpolation methods under various bases (Chebyshev, Hermite, Laguerre, Legendre, Fourier, Sinc) in Julia. This suite of functionality was originally proposed and implemented in Matlab by, J. A. Weideman and S. C. Reddy (see References [1]).

## Theory of Spectral Methods

For information on spectral methods one could consult, __Spectral Methods in Matlab__ (by L.N. Trefethen) for a more practical introduction, and/or __Spectral Methods: Fundamentals in Single Domains__ (by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang) for a more theoretical treatment.

## Usage

DMSuite does not depend on a specific FFT implementation, instead it is built against the abstract interface [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl). Therefore to use a method such as `chebdifft`, one must install and load an FFT implementation such as [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) or [FastTransforms.jl](https://github.com/JuliaApproximation/FastTransforms.jl).

## Example

Solve the BVP: $u'' = \exp(4x)$, subject to the BCs,  $u(-1)=u(1)=0$
```julia
using FastTransforms
using DMSuite
using Plots

N = 16
M = 2
x,DM = chebdif(N+1,M)
D2 = DM[2:N,2:N,2]       # impose boundary conditions
f = exp.(4 .* x[2:N])    # impose boundary conditions

u = D2 \ f     # solve
u = [0;u;0]    # set boundary values

xx = -1:0.01:1
uu = chebint(u,xx)
plot(xx,uu)
```

## Alternative Libraries

An alternative to DMSuite, at a higher level of abstraction and with a wider feature set, is [ApproxFun.jl](https://juliaapproximation.github.io/ApproxFun.jl/latest/)