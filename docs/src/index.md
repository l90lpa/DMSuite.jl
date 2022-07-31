# DMSuite.jl

Spectral differentiation and interpolation methods on polynomial (Chebyshev, Hermite, Laguerre, Legendre), and trigonometric (Fourier, Sinc) basis functions in Julia. This suite of functionality was originally proposed and implemented in Matlab by, J. A. Weideman and S. C. Reddy (see References [1]).

## Theory of Spectral Methods

For information on spectral methods I recommend consulting __Spectral Methods in Matlab__ (by L.N. Trefethen) and/or __Spectral Methods: Fundamentals in Single Domains__ (by C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang)

## Example

Solve the BVP: $u'' = \exp(4x)$ subject to $u(-1)=u(1)=0$
```julia
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