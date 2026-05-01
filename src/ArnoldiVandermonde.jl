module ArnoldiVandermonde

using LinearAlgebra
export ArnoldiBasis, ArnoldiPolynomial, nodes, degree, evaluate, vectors, increment!

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
function ArnoldiBasis(z::AbstractVector=ComplexF64[], m::Integer=0; max_degree=60)
    n = length(z)
    if n < m + 1
        throw(ArgumentError("Number of nodes must be at least one more than the degree"))
    end
    T = eltype(float(z))
    v = Vector{T}(undef, n)
    Q = similar(v, n, max_degree + 1)
    H = similar(v, max_degree + 1, max_degree)
    Q[:, 1] .= 1
    for m in 1:m
        v = z .* Q[:, m]
        for k in 1:m
            Qk = view(Q, :, k)
            H[k, m] = dot(Qk, v) / n
            v .-= H[k, m] * Qk
        end
        H[m+1, m] = norm(v) / sqrt(n)
        @. Q[:, m+1] = v / H[m+1, m]
    end
    return ArnoldiBasis{T}(z, Q, H, m)
end

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
    m = degree(B)
    v = z .* Q[:, m]
    for k in 1:m
        Qk = view(Q, :, k)
        H[k, m] = dot(Qk, v) / n
        v -= H[k, m] * Qk
    end
    H[m+1, m] = norm(v) / sqrt(n)
    @. Q[:, m+1] = v / H[m+1, m]
    return B
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
    function ArnoldiPolynomial{T}(coeff::AbstractVector{T}, basis::ArnoldiBasis{T}) where {T<:RCFloat}
        if length(coeff) != degree(basis) + 1
            throw(ArgumentError("Incompatible coefficient and basis sizes"))
        end
        return new{T}(coeff, basis)
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
    H = p.basis.H
    m = degree(p.basis)
    Q = fill(one(T), m+1)
    for k in 1:m
        v = z .* Q[k]
        for j in 1:k
            v -= H[j, k] * Q[j]
        end
        Q[k+1] = v / H[k+1, k]
        g += p.coeff[k+1] * Q[k+1]
    end
    return g
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


function project(f::Function, a::Real, b::Real; tol=100eps(), max_degree=60)
    x = range(a, b, 100)
    y = f.(x)
    deg = 0
    p = []
    while true
        B = ArnoldiBasis(x, deg)
        p = B \ f
        resid = f.(x) - p.(x)
        if maximum(abs, resid) < tol
            break
        end
        if deg >= max_degree
            @warn "Maximum degree reached without convergence"
            break
        end
        deg += 1
        increment!(B)
    end
    return p
end

end # module