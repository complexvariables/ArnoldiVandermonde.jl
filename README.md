# ArnoldiVandermonde

[![Build Status](https://github.com/complexvariables/ArnoldiVandermonde.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/complexvariables/ArnoldiVandermonde.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/complexvariables/ArnoldiVandermonde.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/complexvariables/ArnoldiVandermonde.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

`ArnoldiVandermonde.jl` is a Julia package for constructing and working with Arnoldi-based Vandermonde matrices. These matrices give a basis for polynomials orthogonalized over a discrete point set, which can enable robust function approximation.

See [Brubeck, Nakatsukasa, & Trefethen (2021)](https://doi.org/10.1137/19M130100X) for an introduction to the Arnoldi–Vandermonde algorithm.

## Installation

To install the package, use the Julia package manager:

```julia
using Pkg
Pkg.add("ArnoldiVandermonde")
```

## Usage

Here is a simple example of how to use `ArnoldiVandermonde.jl`:

```julia
using ArnoldiVandermonde, LinearAlgebra

θ = (0:199) / 100
z = [3cospi(t) + im*sinpi(t) for t in θ]     # points on an ellipse
B = ArnoldiBasis(z, 10)                      # create a basis of degree 10
Q = vectors(B)                               # 200 × 11 matrix of basis vectors
cond(Q)                                      # ≈ 1
f(z) = (z^8 - z^5 + 2) / 3^8
y = f.(z)
norm(Q * (Q \ y) - y)                        # ≈ machine roundoff
```

To evaluate the resulting approximation at off-node points, you can create an `ArnoldiPolynomial`. 

```julia
p = ArnoldiPolynomial(Q \ y, B)
# Small error inside the ellipse:
maximum(abs(f(x) - p(x)) for x in range(-3, 3, 10000))
p = B \ y                                    # does the same thing
p = B \ f                                    # does the same thing
```

The `project` function does a simple iteration on the node set and degree to find a good polynomial approximation to a function over a real interval.

```julia
f(x) = sin(exp(x))
p = project(f, 0, 1)                         # stops at degree 17
p = project(f, BigFloat(0), 1)               # stops at degree 77
x = range(BigFloat(0), 1, 5000)
maximum(abs, f.(x) - p(x))                   # ≈ 3e-72
```

## Notes

- Evaluation of an `ArnoldiPolynomial` directly on an array can be much faster than broadcasting scalar evaluation.
- The same algorithms will work for any vectors that support linear combination and inner products. But this package is written with fixed-dimensional `AbstractVectors` in mind for simplicity. In fact, the resulting basis is a view of a standard `Matrix`.
- If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request on [GitHub](https://github.com/complexvariables/ArnoldiVandermonde.jl).

## Similar work

- The `ArnoldiFit` variant of `fit` provided in [`Polynomials`](https://github.com/JuliaMath/Polynomials.jl/tree/f1ee58659f4ca0000e929945524b3d613bcd8391) has similar functionality. But the present package provides easy access to the basis vectors, automatic iteration of degree and grid in `project`, and (at this writing) is faster in basis construction.
- See also [this thread](https://discourse.julialang.org/t/julia-analog-to-rs-poly-function-orthogonal-polynomials-by-qr-decomposition/128640) about using [Krylov.jl](https://juliaregistries.github.io/General/packages/redirect_to_repo/Krylov) to accomplish this.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
