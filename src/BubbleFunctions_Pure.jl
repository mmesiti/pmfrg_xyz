"""
Stacked version of SigmaType for efficient flavor indexing.
Allows accessing flavor i (1=x, 2=y, 3=z) via indexing: sigma[i]
"""
struct StackedSigmaType{T}
    data::Vector{Array{T,2}}  # [x, y, z] arrays
end

# Conversion from SigmaType to StackedSigmaType
StackedSigmaType(sigma::SigmaType{T}) where {T} =
    StackedSigmaType{T}([sigma.x, sigma.y, sigma.z])

# Allow indexing: sigma[1] gives x, sigma[2] gives y, sigma[3] gives z
Base.getindex(sigma::StackedSigmaType, i::Int) = sigma.data[i]

"""
Context for computing propagators in pure bubble functions.
Contains iSigma and DiSigma (derivative of iSigma) for all three flavors (x, y, z),
used to compute Katanin-corrected propagators.
Uses StackedSigmaType for efficient flavor indexing.
"""
struct PropagatorContext{T}
    iSigma::StackedSigmaType{T}   # From State, indexed by flavor
    DiSigma::StackedSigmaType{T}  # From Deriv, indexed by flavor
    T::T                           # Temperature
end

# Conversion constructor from SigmaType
PropagatorContext(iSigma::SigmaType{T}, DiSigma::SigmaType{T}, temp::T) where {T} =
    PropagatorContext{T}(StackedSigmaType(iSigma), StackedSigmaType(DiSigma), temp)

"""
Compute X-type propagator: -iSKat[i](x, nw) * iG[j](x, nw + ns).
Used in X-bubble functions (particle-particle-like diagrams).

Arguments:
- ctx: PropagatorContext with stacked iSigma and DiSigma
- i, j: Flavor indices (1=x, 2=y, 3=z)
- x: Spatial site index
- nw: Matsubara frequency index
- ns: Frequency shift (from is-1)
"""
@inline function propX(
    ctx::PropagatorContext{T},
    i::Int,
    j::Int,
    x::Int,
    nw::Int,
    ns::Int,
) where {T}
    iSKat_i = iSKat_(ctx.iSigma[i], ctx.DiSigma[i], x, nw, ctx.T)
    iG_j = iG_(ctx.iSigma[j], x, nw + ns, ctx.T)
    return -iSKat_i * iG_j
end

"""
Compute Y-type propagator: -iSKat[i](xi, nw) * iG[j](xj, nw - nt).
Used in Y-bubble functions (particle-hole-like diagrams).

Arguments:
- ctx: PropagatorContext with stacked iSigma and DiSigma
- i, j: Flavor indices (1=x, 2=y, 3=z)
- xi, xj: Spatial site indices (for non-local pairs)
- nw: Matsubara frequency index
- nt: Frequency shift (from it-1)
"""
@inline function propY(
    ctx::PropagatorContext{T},
    i::Int,
    j::Int,
    xi::Int,
    xj::Int,
    nw::Int,
    nt::Int,
) where {T}
    iSKat_i = iSKat_(ctx.iSigma[i], ctx.DiSigma[i], xi, nw, ctx.T)
    iG_j = iG_(ctx.iSigma[j], xj, nw - nt, ctx.T)
    return -iSKat_i * iG_j
end

####################################################
######### PURE VERTEX ACCESS FOR BUBBLE COMPUTATION
####################################################

"""
Pure function to get a single vertex value with flavor transformation.
Returns Gamma[n_transformed, R, s+1, t+1, u+1] after applying flavor transformation.

Arguments:
- Gamma: 5D vertex array [21, Npairs, N, N, N]
- n: Flavor index (1-21, from fd module)
- s, t, u: Frequency indices (0-based, will be converted via ConvertFreqArgs)
- flavTransf: Tuple{Bool,Bool,Bool} for block transformations
- R: Pair index
- N: Frequency mesh size

Returns: Single scalar vertex value
"""
@inline function getVertex(
    Gamma::AbstractArray{T,5},
    n::Int,
    s::Int,
    t::Int,
    u::Int,
    flavTransf::Tuple{Bool,Bool,Bool},
    R::Int,
    N::Int,
) where {T}
    @inbounds begin
        # Apply flavor transformation based on block
        block = div(n + 2, 6)  # 0=diagonal, 1=block1, 2=block2, 3=block3

        if block != 0 && flavTransf[block]
            # Cyclic permutation: (a,b,c,d,e,f) → (d,e,f,a,b,c)
            block_start = 4 + (block - 1) * 6
            offset = n - block_start
            new_offset = (offset + 3) % 6  # Right shift by 3
            n_transf = block_start + new_offset
        else
            n_transf = n
        end

        # Convert frequency arguments (apply abs and cutoffs)
        s_abs, t_abs, u_abs = ConvertFreqArgs(s, t, u, N)

        return Gamma[n_transf, R, s_abs+1, t_abs+1, u_abs+1]
    end
end


"""
Pure bubble functions for PMFRG_xyz.

This file contains 42 pure functions (21 X-type + 21 Y-type) that compute
individual flavor contributions to the bubble diagrams in the FRG flow equations.

Each pure function:
- Takes scalar inputs (Rij, is, it, iu) plus arrays (Gamma, ctx, System, lenIntw)
- Includes the Matsubara frequency (nw) sum internally
- Returns a scalar contribution for one flavor

This refactoring enables clearer understanding of data reuse opportunities.
"""

####################################################
######### X-TYPE BUBBLE FUNCTIONS (21 functions)
####################################################
# X-bubble functions compute particle-particle-like contributions
# Stored in X[1:21, Rij, is, it, iu]

"""
Compute X contribution for flavor yy (fd.yy = 2).
Formula: -V12[yy]*V34[yy]*P[2,2] - V12[yz1]*V34[zy1]*P[3,3] - V12[yx1]*V34[xy1]*P[1,1]

This is a reference implementation for the X-type functions.
"""
function Xyy(
    Rij::Int,
    is::Int,
    it::Int,
    iu::Int,
    Gamma::AbstractArray{T,5},
    ctx::PropagatorContext{T},
    System,
    lenIntw::Int,
) where {T}

    (; Nsum, siteSum, invpairs) = System
    N = size(Gamma, 3)
    ns, nt, nu = is - 1, it - 1, iu - 1

    result = zero(T)

    for nw = -lenIntw:lenIntw-1
        # Compute mixed frequencies
        wpw1, wpw2, _, _, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nw)

        # Frequency arguments for vertices
        s1, t1, u1 = ConvertFreqArgs(ns, wpw1, -wpw2, N)
        s2, t2, u2 = ConvertFreqArgs(ns, -wmw3, -wmw4, N)

        # Sum over lattice sites for this pair
        for k_spl = 1:Nsum[Rij]
            ki = siteSum.ki[k_spl, Rij]
            kj = siteSum.kj[k_spl, Rij]
            m = siteSum.m[k_spl, Rij]
            xk = siteSum.xk[k_spl, Rij]

            # Determine R indices based on flavor transformation
            # flavTransf12[1] = wpw1 * wpw2 > 0
            # flavTransf34[1] = wmw3 * wmw4 < 0
            R12 = (wpw1 * wpw2 > 0) ? invpairs[ki] : ki
            R34 = (wmw3 * wmw4 < 0) ? invpairs[kj] : kj

            # Accumulate contribution (note the minus signs)
            # propX(ctx, i, i, xk, nw, ns) = -(iSKat_(ctx.iSigma[i], ctx.DiSigma[i], xk, nw, ctx.T) * iG_(ctx.iSigma[i], xk, nw + ns, ctx.T))
            result += (
                -Gamma[fd.yy, R12, s1+1, t1+1, u1+1] *
                Gamma[fd.yy, R34, s2+1, t2+1, u2+1] *
                (
                    -m *
                    iSKat_(ctx.iSigma[2], ctx.DiSigma[2], xk, nw, ctx.T) *
                    iG_(ctx.iSigma[2], xk, nw + ns, ctx.T)
                ) -
                Gamma[ifelse(wpw1 * wpw2 > 0, fd.zy1, fd.yz1), R12, s1+1, t1+1, u1+1] *
                Gamma[ifelse(wmw3 * wmw4 < 0, fd.yz1, fd.zy1), R34, s2+1, t2+1, u2+1] *
                (
                    -m *
                    iSKat_(ctx.iSigma[3], ctx.DiSigma[3], xk, nw, ctx.T) *
                    iG_(ctx.iSigma[3], xk, nw + ns, ctx.T)
                ) -
                Gamma[ifelse(wpw1 * wpw2 > 0, fd.xy1, fd.yx1), R12, s1+1, t1+1, u1+1] *
                Gamma[ifelse(wmw3 * wmw4 < 0, fd.yx1, fd.xy1), R34, s2+1, t2+1, u2+1] *
                (
                    -m *
                    iSKat_(ctx.iSigma[1], ctx.DiSigma[1], xk, nw, ctx.T) *
                    iG_(ctx.iSigma[1], xk, nw + ns, ctx.T)
                )
            )
        end
    end

    return result
end

# TODO: Implement remaining 20 X-functions:
# Xxx, Xzz, Xxy1, Xxz1, Xyz1, Xyx1, Xzx1, Xzy1,
# Xxy2, Xxz2, Xyz2, Xyx2, Xzx2, Xzy2,
# Xxy3, Xxz3, Xyz3, Xyx3, Xzx3, Xzy3

####################################################
######### Y-TYPE BUBBLE FUNCTIONS (21 functions)
####################################################
# Y-bubble functions compute particle-hole-like contributions
# Stored in X[22:42, Rij, is, it, iu]
# Only computed for non-local pairs (Rij not in OnsitePairs)

"""
Compute Y contribution for flavor xx (fd.xx = 1).
Formula: (V13[xx]*V24[xx]*P[1,1] + V13[xy2]*V24[xy2]*P[2,2] + V13[xz2]*V24[xz2]*P[3,3])
       + (V31[xx]*V42[xx]*PT[1,1] + V31[xy2]*V42[xy2]*PT[2,2] + V31[xz2]*V42[xz2]*PT[3,3])

This is a reference implementation for the Y-type functions.
"""
function Yxx(
    Rij::Int,
    is::Int,
    it::Int,
    iu::Int,
    Gamma::AbstractArray{T,5},
    ctx::PropagatorContext{T},
    System,
    lenIntw::Int,
) where {T}

    (; Npairs, invpairs, PairTypes, OnsitePairs) = System
    N = size(Gamma, 3)

    # Y-bubble only for non-local pairs
    if Rij in OnsitePairs
        return zero(T)
    end

    (; xi, xj) = PairTypes[Rij]
    ns, nt, nu = is - 1, it - 1, iu - 1

    result = zero(T)

    for nw = -lenIntw:lenIntw-1
        # Compute mixed frequencies (different from X-type)
        _, wpw2, wpw4, wmw1, wmw3, _ = mixedFrequencies(ns, nt, nu, nw)

        # Flavor transformations for 4 vertices
        flavTransf13 = (nt * wmw3 < 0, wmw1 * wmw3 > 0, wmw1 * nt > 0)
        flavTransf24 = (nt * wpw4 < 0, wpw2 * wpw4 > 0, wpw2 * nt > 0)
        flavTransf31 = (nt * wmw1 > 0, wmw3 * wmw1 > 0, wmw3 * nt < 0)
        flavTransf42 = (nt * wpw2 > 0, wpw4 * wpw2 > 0, wpw4 * nt < 0)

        # Frequency arguments
        s13, t13, u13 = ConvertFreqArgs(wmw1, nt, wmw3, N)
        s24, t24, u24 = ConvertFreqArgs(wpw2, -nt, -wpw4, N)
        s31, t31, u31 = ConvertFreqArgs(wmw3, nt, -wmw1, N)
        s42, t42, u42 = ConvertFreqArgs(-wpw4, -nt, wpw2, N)

        # R indices based on flavor transformations
        R13 = flavTransf13[1] ? invpairs[Rij] : Rij
        R24 = flavTransf24[1] ? invpairs[Rij] : Rij
        R31 = flavTransf31[1] ? invpairs[Rij] : Rij
        R42 = flavTransf42[1] ? invpairs[Rij] : Rij

        # Get all needed vertex values (4 vertices, 3 flavors each = 12 values)
        V13_xx = getVertex(Gamma, fd.xx, s13, t13, u13, flavTransf13, R13, N)
        V13_xy2 = getVertex(Gamma, fd.xy2, s13, t13, u13, flavTransf13, R13, N)
        V13_xz2 = getVertex(Gamma, fd.xz2, s13, t13, u13, flavTransf13, R13, N)

        V24_xx = getVertex(Gamma, fd.xx, s24, t24, u24, flavTransf24, R24, N)
        V24_xy2 = getVertex(Gamma, fd.xy2, s24, t24, u24, flavTransf24, R24, N)
        V24_xz2 = getVertex(Gamma, fd.xz2, s24, t24, u24, flavTransf24, R24, N)

        V31_xx = getVertex(Gamma, fd.xx, s31, t31, u31, flavTransf31, R31, N)
        V31_xy2 = getVertex(Gamma, fd.xy2, s31, t31, u31, flavTransf31, R31, N)
        V31_xz2 = getVertex(Gamma, fd.xz2, s31, t31, u31, flavTransf31, R31, N)

        V42_xx = getVertex(Gamma, fd.xx, s42, t42, u42, flavTransf42, R42, N)
        V42_xy2 = getVertex(Gamma, fd.xy2, s42, t42, u42, flavTransf42, R42, N)
        V42_xz2 = getVertex(Gamma, fd.xz2, s42, t42, u42, flavTransf42, R42, N)

        # Compute propagators: P[i,j] and PT[i,j] (transposed)
        P11 = propY(ctx, 1, 1, xi, xj, nw, nt)
        P22 = propY(ctx, 2, 2, xi, xj, nw, nt)
        P33 = propY(ctx, 3, 3, xi, xj, nw, nt)

        PT11 = propY(ctx, 1, 1, xj, xi, nw, nt)
        PT22 = propY(ctx, 2, 2, xj, xi, nw, nt)
        PT33 = propY(ctx, 3, 3, xj, xi, nw, nt)

        # Accumulate contribution (P term + PT term, note positive signs)
        result += (
            (V13_xx * V24_xx * P11 + V13_xy2 * V24_xy2 * P22 + V13_xz2 * V24_xz2 * P33) + (
                V31_xx * V42_xx * PT11 +
                V31_xy2 * V42_xy2 * PT22 +
                V31_xz2 * V42_xz2 * PT33
            )
        )
    end

    return result
end

# TODO: Implement remaining 20 Y-functions:
# Yyy, Yzz, Yxy1, Yxz1, Yyz1, Yyx1, Yzx1, Yzy1,
# Yxy2, Yxz2, Yyz2, Yyx2, Yzx2, Yzy2,
# Yxy3, Yxz3, Yyz3, Yyx3, Yzx3, Yzy3
