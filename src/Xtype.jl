"""
    AbstractXIndexMapping

Abstract type for index mappings between linear index `i` and multi-index `(n, Rij, is, it, iu)`.

A concrete implementation must provide:
- `linear_index(mapping, n, Rij, is, it, iu)`: Convert multi-index to linear index
- `multi_index(mapping, i)`: Convert linear index to multi-index (n, Rij, is, it, iu)
- `total_length(mapping)`: Total number of valid elements

# Constraints
- Elements are only stored when (is-1) + (it-1) + (iu-1) is odd
- The mapping must be bijective on the valid domain
"""
abstract type AbstractXIndexMapping end

"""
    linear_index(mapping::AbstractXIndexMapping, n::Int, Rij::Int, is::Int, it::Int, iu::Int) -> Int

Convert multi-index (n, Rij, is, it, iu) to linear index i.

Must satisfy: 1 <= i <= total_length(mapping)
"""
function linear_index end

"""
    multi_index(mapping::AbstractXIndexMapping, i::Int) -> (n, Rij, is, it, iu)

Convert linear index i to multi-index (n, Rij, is, it, iu).

Must satisfy: multi_index(mapping, linear_index(mapping, n, Rij, is, it, iu)) == (n, Rij, is, it, iu)
"""
function multi_index end

"""
    total_length(mapping::AbstractXIndexMapping) -> Int

Return the total number of valid elements in the mapping.
"""
function total_length end

"""
    is_valid_multi_index(is::Int, it::Int, iu::Int) -> Bool

Check if a multi-index corresponds to a valid element.
Elements with (is-1) + (it-1) + (iu-1) even are invalid.
"""
@inline function is_valid_multi_index(is::Int, it::Int, iu::Int)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    return (ns + nt + nu) % 2 == 1
end

"""
    XVector{T, M <: AbstractXIndexMapping}

A vector storing X values with customizable index mapping.

# Fields
- `data::Vector{T}`: Linear storage of values
- `mapping::M`: Index mapping strategy

# Usage
```julia
# Create with a mapping
mapping = SomeXIndexMapping(n_flavors, n_pairs, N)
xvec = XVector{Float64}(mapping)

# Access by multi-index
val = xvec[n, Rij, is, it, iu]
xvec[n, Rij, is, it, iu] = new_val

# Access underlying data directly
val = xvec.data[i]
```
"""
struct XVector{T,M<:AbstractXIndexMapping}
    data::Vector{T}
    mapping::M

    function XVector{T}(mapping::M) where {T,M<:AbstractXIndexMapping}
        len = total_length(mapping)
        data = zeros(T, len)
        new{T,M}(data, mapping)
    end
end

# Multi-index access
Base.@propagate_inbounds function Base.getindex(
    xv::XVector,
    n::Int,
    Rij::Int,
    is::Int,
    it::Int,
    iu::Int,
)
    i = linear_index(xv.mapping, n, Rij, is, it, iu)
    return xv.data[i]
end

Base.@propagate_inbounds function Base.setindex!(
    xv::XVector,
    val,
    n::Int,
    Rij::Int,
    is::Int,
    it::Int,
    iu::Int,
)
    i = linear_index(xv.mapping, n, Rij, is, it, iu)
    xv.data[i] = val
    return val
end

# Standard array interface
Base.size(xv::XVector) = size(xv.data)
Base.length(xv::XVector) = length(xv.data)
Base.eltype(::Type{XVector{T,M}}) where {T,M} = T
Base.similar(xv::XVector{T,M}) where {T,M} = XVector{T}(xv.mapping)

# Fill operations
Base.fill!(xv::XVector, val) = (fill!(xv.data, val); xv)
Base.zero(xv::XVector{T,M}) where {T,M} = XVector{T}(xv.mapping)

"""
    DefaultXIndexMapping

Default index mapping for X array, optimized for N even.
Skips elements where (is-1) + (it-1) + (iu-1) is even.

Layout: n varies fastest, then Rij, then is, then it, then iu.

For N even, exactly half of the (is, it, iu) combinations are valid,
and the count can be computed using simple arithmetic.
"""
struct DefaultXIndexMapping <: AbstractXIndexMapping
    n_flavors::Int
    n_pairs::Int
    N::Int

    function DefaultXIndexMapping(n_flavors::Int, n_pairs::Int, N::Int)
        @assert iseven(N) "N must be even for DefaultXIndexMapping"
        new(n_flavors, n_pairs, N)
    end
end

"""
    count_valid_before(N::Int, is::Int, it::Int, iu::Int) -> Int

Count valid (is, it, iu) combinations before (is, it, iu).
For N even, this uses arithmetic based on parity arguments.
"""
@inline function count_valid_before(N::Int, is::Int, it::Int, iu::Int)
    # Complete iu layers before this one
    complete_iu_layers = iu - 1
    # Each complete iu layer has N*N elements, exactly half valid
    count = complete_iu_layers * (N * N ÷ 2)

    # Within current iu layer, complete it rows before this one
    complete_it_rows = it - 1
    # Each complete it row has N elements, exactly half valid
    count += complete_it_rows * (N ÷ 2)

    # Within current it row, count valid is values before this one
    # is is valid if (is-1) + (it-1) + (iu-1) is odd
    # Regardless of which parity is required, counting elements in an
    # alternating sequence {1,3,5,...} or {2,4,6,...} before position is
    # always gives (is-1)÷2
    count += (is - 1) ÷ 2

    return count
end

@inline function linear_index(
    m::DefaultXIndexMapping,
    n::Int,
    Rij::Int,
    is::Int,
    it::Int,
    iu::Int,
)
    @boundscheck begin
        @assert 1 <= n <= m.n_flavors "n must be in 1:$(m.n_flavors)"
        @assert 1 <= Rij <= m.n_pairs "Rij must be in 1:$(m.n_pairs)"
        @assert 1 <= is <= m.N "is must be in 1:$(m.N)"
        @assert 1 <= it <= m.N "it must be in 1:$(m.N)"
        @assert 1 <= iu <= m.N "iu must be in 1:$(m.N)"
        @assert is_valid_multi_index(is, it, iu) "Invalid: (is-1)+(it-1)+(iu-1) must be odd"
    end

    # Number of valid (is, it, iu) blocks before this one
    num_blocks_before = count_valid_before(m.N, is, it, iu)

    # Base offset (1-indexed)
    base_offset = num_blocks_before * m.n_flavors * m.n_pairs

    # Offset within current (is, it, iu) block: n fastest, then Rij
    block_offset = (Rij - 1) * m.n_flavors + n

    return base_offset + block_offset
end

@inline function multi_index(m::DefaultXIndexMapping, i::Int)
    @boundscheck @assert 1 <= i <= total_length(m) "Linear index out of bounds"

    N = m.N
    half_N = N ÷ 2

    # Size of each (is, it, iu) block
    block_size = m.n_flavors * m.n_pairs

    # Which valid block (0-indexed)
    block_idx = (i - 1) ÷ block_size
    offset_in_block = (i - 1) % block_size

    # Each iu layer has N*N/2 valid combinations
    elements_per_iu = N * half_N
    iu_idx = block_idx ÷ elements_per_iu
    remaining = block_idx % elements_per_iu

    # Within iu layer, each it row has N/2 valid combinations
    it_idx = remaining ÷ half_N
    is_offset = remaining % half_N

    iu = iu_idx + 1
    it = it_idx + 1

    # Find is based on parity requirement
    # We want (is-1) + (it-1) + (iu-1) to be odd
    # So (is-1) must have parity: odd - (it-1) - (iu-1) mod 2
    target_parity = (1 - (it - 1) - (iu - 1)) % 2

    if target_parity == 0
        # ns = is-1 must be even: ns = 0, 2, 4, ... → is = 1, 3, 5, ...
        is = 2 * is_offset + 1
    else
        # ns = is-1 must be odd: ns = 1, 3, 5, ... → is = 2, 4, 6, ...
        is = 2 * is_offset + 2
    end

    # Extract n and Rij from offset_in_block
    Rij = offset_in_block ÷ m.n_flavors + 1
    n = offset_in_block % m.n_flavors + 1

    return (n, Rij, is, it, iu)
end

@inline function total_length(m::DefaultXIndexMapping)
    # Total valid (is, it, iu) combinations: N^3 / 2 for even N
    num_valid_combinations = (m.N * m.N * m.N) ÷ 2
    return num_valid_combinations * m.n_flavors * m.n_pairs
end

# Export public interface
export AbstractXIndexMapping, DefaultXIndexMapping
export linear_index, multi_index, total_length, is_valid_multi_index
export XVector
