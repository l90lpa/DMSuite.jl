# DMSuite.jl

Spectral differentiation and interpolation methods under various bases (Chebyshev, Hermite, Laguerre, Legendre, Fourier, Sinc) in Julia. This suite of functionality was originally proposed and implemented in Matlab by, J. A. Weideman and S. C. Reddy [1].

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

## References

1. J. A. Weideman and S. C. Reddy. 2000. __A MATLAB differentiation matrix suite__. ACM Trans. Math. Softw. 26, 4 (Dec. 2000), 465–519. [https://doi.org/10.1145/365723.365727](https://doi.org/10.1145/365723.365727)
1. W. S. Don and A. Solomonoff. 1995. __Accuracy and Speed in Computing the Chebyshev Collocation Derivative__. SIAM J. Sci. Comput. 16, 6 (Nov. 1995), 1253–1268. [https://doi.org/10.1137/0916073](https://doi.org/10.1137/0916073)
1. R. Baltensperger and M. R. Trummer. 2003. __Spectral Differencing with a Twist__. SIAM J. Sci. Comput. 24, 5 (Jan. 2003), 1465–1487. [https://doi.org/10.1137/S1064827501388182](https://doi.org/10.1137/S1064827501388182)
1. P. Henrici. __Essentials of Numerical Analysis__. 1982. Wiley. (barycentric formula see page 252)
1. P. Henrici. __Applied & Computational Complex Analysis III__. 1986. Wiley. (barycentric formula see page 46)
