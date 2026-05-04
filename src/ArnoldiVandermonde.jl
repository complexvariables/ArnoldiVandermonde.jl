module ArnoldiVandermonde

using LinearAlgebra
export ArnoldiBasis, ArnoldiPolynomial, nodes, degree, vectors,
    evaluate, evaluate!, increment!, project

const RCFloat{T} = Union{T, Complex{T}} where {T<:AbstractFloat}

"""
    ArnoldiBasis (type)

Well-conditioned representation for polynomials on a discrete point set.

# Fields
- `nodes`: evaluation points
- `Q`: orthogonal basis
- `H`: orthogonalization coefficients
"""
mutable struct ArnoldiBasis{T}
    nodes::Vector{T}
    Q::Matrix{T}    # orthonormal basis vectors
    H::UpperHessenberg{T, Matrix{T}}    # upper Hessenberg matrix of orthogonalization coefficients
    degree::Int
    function ArnoldiBasis{T}(nodes::AbstractVector{T}, Q::AbstractMatrix{T}, H::AbstractMatrix{T}, degree::Int) where {T<:RCFloat}
        size(Q, 1) == length(nodes) || throw(ArgumentError("Number of rows in Q must match number of nodes"))
        return new{T}(nodes, Q, UpperHessenberg(H), degree)
    end
    function ArnoldiBasis{T}(nodes::AbstractVector, Q::AbstractMatrix, H::AbstractMatrix, degree::Int) where {T<:RCFloat}
        nodes = convert.(T, nodes)
        Q = convert.(T, Q)
        H = convert.(T, H)
        return ArnoldiBasis{T}(nodes, Q, UpperHessenberg(H), degree)
    end
    function ArnoldiBasis{T}(B::ArnoldiBasis{S}) where {T<:RCFloat,S}
        return ArnoldiBasis{T}(B.nodes, B.Q, B.H, B.degree)
    end
end
Base.eltype(B::ArnoldiBasis) = eltype(B.nodes)
degree(b::ArnoldiBasis) = b.degree
nodes(b::ArnoldiBasis) = b.nodes
vectors(b::ArnoldiBasis) = view(b.Q, :, 1:b.degree+1)

# Main basis constructor
"""
    ArnoldiBasis(z, m)

Construct an ArnoldiBasis of degree `m` on the nodes `z`. The nodes must be distinct and the number of nodes must be at least one more than the degree.
"""
function ArnoldiBasis(z::AbstractVector=ComplexF64[], m::Integer=0; max_degree=max(m, 60))
    n = length(z)
    if n < m + 1
        throw(ArgumentError("Number of nodes must be at least one more than the degree"))
    end
    T = eltype(float(z))
    v = Vector{T}(undef, n)
    Q = similar(v, n, max_degree + 1)
    H = similar(v, max_degree + 1, max_degree)
    Q[:, 1] .= 1
    for k in 1:m
        _increment!(Q, H, z, k, n)
    end
    return ArnoldiBasis{T}(z, Q, H, m)
end

"""
    increment!(B)

Increment the degree of the ArnoldiBasis `B` by one.
"""
function increment!(B::ArnoldiBasis{T}) where {T}
    z, Q, H = nodes(B), B.Q, B.H
    n = length(z)
    if degree(B) == n
        throw(ArgumentError("Cannot grow basis beyond number of nodes"))
    end
    if degree(B) == size(Q, 2) - 1
        throw(ArgumentError("Maximum degree reached; cannot grow basis"))
    end
    B.degree += 1
    _increment!(Q, H, z, degree(B), n)
    return B
end

# Perform one pass of Arnoldi orthogonalization
function _increment!(Q, H, z, m, n)
    v = view(Q, :, m+1)
    v .= z .* view(Q, :, m)
    for k in 1:m
        Qk = view(Q, :, k)
        H[k, m] = dot(Qk, v) / n
        axpy!(-H[k, m], Qk, v)
    end
    H[m+1, m] = norm(v) / sqrt(n)
    lmul!(1/H[m+1, m], v)
    return nothing
end

# COV_EXCL_START
function Base.show(io::IO, ::MIME"text/plain", B::ArnoldiBasis)
    ioc = IOContext(io, :compact => get(io, :compact, true))
    n = length(B.nodes)
    deg = degree(B)
    print(ioc, "Arnoldi basis of degree $deg on $n nodes")
end

function Base.show(io::IO, B::ArnoldiBasis)
    n = length(B.nodes)
    m = size(B.Q, 2)
    print(IOContext(io, :compact => true), "$n×$m Arnoldi basis")
end
# COV_EXCL_STOP

"""
    ArnoldiPolynomial (type)

Polynomial representation using an ArnoldiBasis.

# Fields
- `coeff`: vector of coefficients
- `basis`: ArnoldiBasis for the polynomial
"""
struct ArnoldiPolynomial{T} <: Function
    coeff::Vector{T}
    basis::ArnoldiBasis{T}
    tmp::Vector{T}   # temporary workspace for evaluation
    function ArnoldiPolynomial{T}(coeff::AbstractVector{T}, basis::ArnoldiBasis{T}) where {T<:RCFloat}
        if length(coeff) != degree(basis) + 1
            throw(ArgumentError("Incompatible coefficient and basis sizes"))
        end
        tmp = similar(coeff)
        return new{T}(coeff, basis, tmp)
    end
    function ArnoldiPolynomial{T}(p::ArnoldiPolynomial{S}) where {T<:RCFloat,S}
        return ArnoldiPolynomial{T}(convert.(T, p.coeff), ArnoldiBasis{T}(p.basis))
    end
end

"""
    ArnoldiPolynomial(coeff, basis)

Construct an ArnoldiPolynomial with coefficients `coeff` and basis `basis`.
"""
function ArnoldiPolynomial(
    coeff::AbstractVector=[0im],
    basis::ArnoldiBasis{S}=ArnoldiBasis()
    ) where {S}
    T = promote_type(eltype(coeff), S)
    coeff = convert.(T, coeff)
    return ArnoldiPolynomial{T}(coeff, ArnoldiBasis{T}(basis))
end

Base.eltype(p::ArnoldiPolynomial) = eltype(p.coeff)
Base.length(p::ArnoldiPolynomial) = length(p.coeff)
degree(p::ArnoldiPolynomial) = length(p.coeff) - 1
nodes(p::ArnoldiPolynomial) = nodes(p.basis)

# COV_EXCL_START
function Base.show(io::IO, ::MIME"text/plain", p::ArnoldiPolynomial)
    ioc = IOContext(io,:compact => get(io, :compact, true))
    len = length(p)
    n = length(nodes(p))
    if len==0
        print(ioc, "Empty $(typeof(p)) Arnoldi polynomial")
    else
        print(ioc, "Arnoldi polynomial of degree $(len-1) on $n nodes")
    end
end

function Base.show(io::IO, p::ArnoldiPolynomial)
    deg = length(p) - 1
    print(IOContext(io, :compact => true), "Arnoldi polynomial of degree $deg")
end
# COV_EXCL_STOP


"""
    p(z)
    evaluate(p, z)

Evaluate the ArnoldiPolynomial `p` at `z`.
"""
(p::ArnoldiPolynomial)(z) = evaluate(p, z)

function evaluate(p::ArnoldiPolynomial{T}, z::Number) where {T}
    g = p.coeff[1]
    m = degree(p.basis)
    q = p.tmp
    H = p.basis.H
    q[1] = one(T)
        for k in 1:m
        v = q[k+1]
        v = z * q[k]
        for j in 1:k
            v -= H[j, k] * q[j]
        end
        q[k+1] = v / H[k+1, k]
        g += p.coeff[k+1] * q[k+1]
    end
    return g
end

function evaluate(p::ArnoldiPolynomial{T}, z::AbstractArray) where {T}
    g = fill(p.coeff[1], length(z))
    evaluate!(g, p, z)
    return g
end

"""
    evaluate!(g::AbstractArray, p::ArnoldiPolynomial, z::AbstractArray)

Evaluate the ArnoldiPolynomial `p` at the points `z`, using `g` to calculate the result in-place.
"""
function evaluate!(g, p::ArnoldiPolynomial{T}, z::AbstractArray) where {T}
    g .= p.coeff[1]
    m = degree(p.basis)
    Q = Matrix{T}(undef, length(z), m+1)
    Q[:, 1] .= 1
    H = p.basis.H
    for k in 1:m
        v = view(Q, :, k+1)
        v .= z .* view(Q, :, k)
        for j in 1:k
            axpy!(-H[j, k], view(Q, :, j), v)
        end
        lmul!(1/H[k+1, k], v)
        axpy!(p.coeff[k+1], v, g)
    end
    return nothing
end

"""
    \\(B::ArnoldiBasis, data)

Find the least-squares projection of the data onto the basis `B`, returning an ArnoldiPolynomial.
The data can be either a vector of values at the nodes or a function that will be evaluated at the nodes.

# Example
````julia-repl
julia> z = cispi.(range(0, 2, 200));   # discretize unit circle

julia> B = ArnoldiBasis(z, 10);

julia> p = B \\ cos
Arnoldi polynomial of degree 10 on 200 nodes

julia> maximum(abs, p.(z) - cos.(z))
2.1199597752184118e-9
````
"""
Base.:\(B::ArnoldiBasis, f::Function) = B \ f.(nodes(B))

function Base.:\(B::ArnoldiBasis, y::AbstractVector)
    length(y) == length(nodes(B)) || throw(ArgumentError("Length of data must match number of nodes"))
    c = vectors(B) \ y
    return ArnoldiPolynomial(c, B)
end

"""
    project(f, a, b; tol=100ε, max_degree=100, error_norm=x -> norm(x, Inf))

Project a function `f` in the least-squares sense onto the interval `[a, b]`, returning an ArnoldiPolynomial.
The projection iteratively grows the basis until the infinity norm of the error is less than `tol` or the degree reaches `max_degree`. The `error_norm` function can be used to specify a different norm for measuring the error.
"""
function project(f::Function, a::Real, b::Real;
    tol = 100eps(promote_type(typeof(float(a)), typeof(float(b)))),
    max_degree = 100,
    error_norm = x -> norm(x, Inf)
    )
    deg = 0
    n = 32
    x = range(a, b, n+1)
    y = f.(x)
    px = similar(y)
    B = ArnoldiBasis(x, deg; max_degree)
    p = []
    while true
        p = B \ y
        evaluate!(px, p, x)
        if error_norm(y - px) < tol
            break
        end
        if deg >= max_degree
            @warn "Maximum degree reached without convergence"
            break
        end
        deg += 1
        increment!(B)
        if length(x) < 3deg
            n *= 2
            x = range(a, b, n+1)
            y = f.(x)
            B = ArnoldiBasis(x, deg; max_degree)
            px = similar(y)
        end
    end
    return p
end

end # module