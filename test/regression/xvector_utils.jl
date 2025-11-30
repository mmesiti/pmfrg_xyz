# Utility functions for converting between Array{T,5} and XVector formats
# Used for testing XVector with existing regression data

using PMFRG_xyz

"""
    array_to_xvector(X_array::Array{T,5}, Par) -> XVector{T}

Convert X from Array{T,5} format to XVector format.
Only copies valid elements where (is-1) + (it-1) + (iu-1) is odd.
"""
function array_to_xvector(X_array::Array{T,5}, Par) where {T}
    # Extract dimensions
    n_flavors, n_pairs, N, _, _ = size(X_array)
    @assert n_flavors == 42 "Expected 42 flavors"
    @assert n_pairs == Par.System.Npairs
    @assert N == Par.NumericalParams.N

    # Create XVector mapping
    mapping = PMFRG_xyz.DefaultXIndexMapping(n_flavors, n_pairs, N)
    X_xvec = PMFRG_xyz.XVector{T}(mapping)

    # Copy valid elements from Array to XVector
    for iu in 1:N, it in 1:N, is in 1:N
        if !PMFRG_xyz.is_valid_multi_index(1, 1, is, it, iu)
            continue
        end
        for Rij in 1:n_pairs, n in 1:n_flavors
            X_xvec[n, Rij, is, it, iu] = X_array[n, Rij, is, it, iu]
        end
    end

    return X_xvec
end

"""
    xvector_to_array(X_xvec::XVector{T}, Par) -> Array{T,5}

Convert X from XVector format back to Array{T,5} format.
Only copies valid elements; invalid elements are set to zero.
"""
function xvector_to_array(X_xvec::PMFRG_xyz.XVector{T,MType}, Par) where {T,MType}
    # Extract parameters
    n_flavors = 42
    n_pairs = Par.System.Npairs
    N = Par.NumericalParams.N

    # Create Array with zeros
    X_array = zeros(T, n_flavors, n_pairs, N, N, N)

    # Copy valid elements from XVector to Array
    for iu in 1:N, it in 1:N, is in 1:N
        if !PMFRG_xyz.is_valid_multi_index(1, 1, is, it, iu)
            continue
        end
        for Rij in 1:n_pairs, n in 1:n_flavors
            X_array[n, Rij, is, it, iu] = X_xvec[n, Rij, is, it, iu]
        end
    end

    return X_array
end

export array_to_xvector, xvector_to_array
