module PMFRG_xyz

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

using RecursiveArrayTools
using SpinFRGLattices,
    OrdinaryDiffEqLowOrderRK, DiffEqCallbacks, RecursiveArrayTools, StructArrays
using SpinFRGLattices.StaticArrays
using Unroll
using MuladdMacro
using FastBroadcast

setZero!(a::AbstractArray{T,N}) where {T,N} = fill!(a, zero(T))

function setZero!(PartArr::ArrayPartition)
    for arr in PartArr.x
        fill!(arr, 0.0)
    end
end

"""Recursively sets structure to zero"""
function setZero!(a::T) where {T}
    for f in fieldnames(T)
        setZero!(getfield(a, f))
    end
    return a
end

# In my convention instead of γ I use iΣ
# Since iΣ now carries three flavors I create a struct for it
struct SigmaType{T}
    x::Array{T,2}
    y::Array{T,2}
    z::Array{T,2}
end

# Previously there existed a Γ::VertexType with VertexType containing
# the three Γ-flavors. Since here I have 21 I opted to add one
# array-dimension as opposed to enlarge the struct.
struct StateType{T}
    f_int::Vector{T}
    iSigma::SigmaType{T}
    Gamma::Array{T,5}
end

# the XYZ model may give different Χ_x, Χ_y, X_z
struct Observables{T}
    Chi_x::Vector{T}
    Chi_y::Vector{T}
    Chi_z::Vector{T}
end

# np_vec is removed because
# ns = np_vec[is] is the same
# as simply ns = is - 1
struct NumericalParams{T<:Real}
    N::Int

    accuracy::T
    temp_min::T
    temp_max::T

    lenIntw::Int
    lenIntw_acc::Int
end

# Remnant from the old code that I havent implemented yet
struct OptionParams
    use_symmetry::Bool
    minimal_output::Bool
end

# The code doesnt work for some reason if I name this struct
# OneLoopParams
struct OneLoopParams{T,SType}
    System::SType
    NumericalParams::NumericalParams{T}
    Options::OptionParams
end

# Similar to Gamma I give X an extra dimension as opposed to create
# A BubbleType struct for it
struct OneLoopWorkspace{T,ParType}
    State::StateType{T}
    Deriv::StateType{T}
    X::Array{T,5}
    Par::ParType
end

# For a general Vertex there can be 3^4 = 81 flavor combinations
# In the XYZ model the SO(3) symmetry breaks down to a residual Klein-4 Symmetry
# This means that the Vertex function can only depend on two distinct flavors at most
# This gives 21 different Vertex functions.
# In my convention I dont use X and ̃X but just one big array called X.
# If I need to acces the ̃X part (which in my convention I name Y) I just go X[21 + flavor]
getVDims(Par) = (
    21,
    Par.System.Npairs,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
)
getBubbleVDims(Par) = (
    42,
    Par.System.Npairs,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
    Par.NumericalParams.N,
)
_getFloatType(Par) = typeof(Par.NumericalParams.accuracy)

function SigmaType(NUnique::Int, N::Int, type = Float64)
    return SigmaType(
        zeros(type, NUnique, N),
        zeros(type, NUnique, N),
        zeros(type, NUnique, N),
    )
end
SigmaType(Par) = SigmaType(Par.System.Npairs, Par.NumericalParams.N)

function StateType(NUnique::Int, N::Int, VDims::Tuple, type = Float64)
    return StateType(zeros(type, NUnique), SigmaType(tpye, NUnique, N), zeros(type, VDims))
end
StateType(Par) =
    StateType(Par.System.NUnique, Par.NumericalParams.N, getVDims(Par), _getFloatType(Par))
StateType(f_int, iSigma_x, iSigma_y, iSigma_z, Gamma) =
    StateType(f_int, SigmaType(iSigma_x, iSigma_y, iSigma_z), Gamma)
RecursiveArrayTools.ArrayPartition(x) =
    ArrayPartition(x.f_int, x.iSigma.x, x.iSigma.y, x.iSigma.z, x.Gamma)
StateType(Arr::ArrayPartition) = StateType(Arr.x...)

# The constructor of this is just blind-copied. To this day I dont really understand
# the purpose of lenIntw and lenIntw_acc
function NumericalParams(;
    N::Integer = 24,
    accuracy = 1e-6,
    temp_min = exp(-10.0),
    temp_max = exp(10.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = 2 * maximum((N, lenIntw)),
)

    return NumericalParams(N, accuracy, temp_min, temp_max, lenIntw, lenIntw_acc)
end

function OneLoopWorkspace(State, Deriv, X, Par)
    setZero!(Deriv)
    setZero!(X)

    return OneLoopWorkspace(StateType(State.x...), StateType(Deriv.x...), X, Par)
end

OptionParams(; use_symmetry::Bool = true, MinimalOutput::Bool = false, kwargs...) =
    OptionParams(use_symmetry, MinimalOutput)
Params(System; kwargs...) =
    OneLoopParams(System, NumericalParams(; kwargs...), OptionParams(; kwargs...))

#############################################################
######### PROPAGATORS ## PROPAGATORS ## PROPAGATORS #########
#############################################################

### Propagators will depend on an additional flavor
### Instead of modifying the propagators, I will simply use them
### as V_, by doing iG_(iSigma.x, ...)

function get_w(nw, T)
    return pi * (2 * nw + 1)
end

function get_sign_iw(nw::Integer, N::Integer)
    # s = sign(nw)
    nw_bounds = min(nw, N - 1)  ### used to be min(abs(nw),...), but nw is set positive in iSigma_
    return nw_bounds + 1        ### used to be s * ...
end

### Sigma inputted as State.iSigma, which is Array{T, 2}
function iSigma_(iSigma::AbstractArray, x::Integer, nw::Integer)
    N = size(iSigma, 2)
    s = 1
    if nw < 0
        nw = -nw - 1
        s = -1
    end
    iw = get_sign_iw(nw, N)
    return s * iSigma[x, iw]
end

function iG_(iSigma::AbstractArray, x::Integer, nw::Integer, T::Real)
    w = get_w(nw, T)
    return 1.0 / (w * sqrt(T) + iSigma_(iSigma, x, nw))
end

### by differentiating the above inverse by T
function iS_(iSigma::AbstractArray, x::Integer, nw::Integer, T::Real)
    w = get_w(nw, T)
    return -iG_(iSigma, x, nw, T)^2 * w / (2.0 * sqrt(T))
end

### Katanin requires (d/dΛ)iΣ, which in the original code is iSigma_(DSigma, x, nw)
### might be wrong here though.
function iSKat_(
    iSigma::AbstractArray,
    DSigma::AbstractArray,
    x::Integer,
    nw::Integer,
    T::Real,
)
    w = get_w(nw, T)
    return -iG_(iSigma, x, nw, T)^2 * (w / (2.0 * sqrt(T)) + iSigma_(DSigma, x, nw))
end

####################################################
######### VERTICES ## VERTICES ## VERTICES #########
####################################################

# In the Heisenberg case these are the Vertex' Symmetries 
#     s <--> -s
#     t <--> -t, i <--> j
#     u <--> -u, i <--> j
# In the XYZ model a change of frequency sign also means a change
# of flavor type. I separate the Vertex flavors into four blocks.
# Transformations of flavors only transform within those blocks.
function ConvertFreqArgs(ns, nt, nu, Nw)
    ns, nt, nu = abs.((ns, nt, nu))

    ns = min(ns, Nw - 1 - (ns + Nw - 1) % 2) ### weird cutoff, idk why
    nt = min(nt, Nw - 1 - (nt + Nw - 1) % 2)
    nu = min(nu, Nw - 1 - (nu + Nw - 1) % 2)

    return ns, nt, nu
end

using LinearAlgebra
using SparseArrays


function V_(
    Vertex::AbstractArray,
    n::Int,
    ns::Int,
    nt::Int,
    nu::Int,
    isFlavorTransform::Tuple{Bool,Bool,Bool},
    Rij::Integer,
    Rji::Integer,
    N::Integer,
)

    @inbounds begin
        # isFlavorTransform = (nt * nu < 0, ns * nu < 0, ns * nt < 0)
        block = div(n + 2, 6)

        if (block != 0 && isFlavorTransform[block])
            # This transforms a block (a,b,c,d,e,f) into (d,e,f,a,b,c)
            # The layout of the flidx-module is hence *very* important

            block_start = 4 + (block - 1) * 6
            offset = n - block_start
            new_offset = (offset + 3) % 6 # cyclic permutation, right shift by 3

            n_transf = block_start + new_offset
        else
            n_transf = n
        end

        ns, nt, nu = ConvertFreqArgs(ns, nt, nu, N)
        Rij = ifelse(isFlavorTransform[1], Rji, Rij)
        return Vertex[n_transf, Rij, ns+1, nt+1, nu+1]
    end
end

"Optimized, in-place version of V_ to be used in addX! and addY!"
@inline function Vert!(V, Gamma, s, t, u, flavTransf, R)
    G = @view Gamma[:, R, s+1, t+1, u+1]

    @unroll for i = 1:3
        V[i] = G[i]
    end

    if flavTransf[1]
        @unroll for i = 4:6
            V[i] = G[i+3]
        end
        @unroll for i = 7:9
            V[i] = G[i-3]
        end
    else
        @unroll for i = 4:9
            V[i] = G[i]
        end
    end

    if flavTransf[2]
        @unroll for i = 10:12
            V[i] = G[i+3]
        end
        @unroll for i = 13:15
            V[i] = G[i-3]
        end
    else
        @unroll for i = 10:15
            V[i] = G[i]
        end
    end

    if flavTransf[3]
        @unroll for i = 16:18
            V[i] = G[i+3]
        end
        @unroll for i = 19:21
            V[i] = G[i-3]
        end
    else
        @unroll for i = 16:21
            V[i] = G[i]
        end
    end
end


function mixedFrequencies(ns, nt, nu, nwpr)
    nw1 = Int((ns + nt + nu - 1) / 2)
    nw2 = Int((ns - nt - nu - 1) / 2)
    nw3 = Int((-ns + nt - nu - 1) / 2)
    nw4 = Int((-ns - nt + nu - 1) / 2)

    wpw1 = nwpr + nw1 + 1
    wpw2 = nwpr + nw2 + 1
    wpw3 = nwpr + nw3 + 1
    wpw4 = nwpr + nw4 + 1
    wmw1 = nwpr - nw1
    wmw2 = nwpr - nw2
    wmw3 = nwpr - nw3
    wmw4 = nwpr - nw4

    return wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4
end

# This defines the Vertex flavors. module was the fastest option
module flidx
const xx = 1
const yy = 2
const zz = 3
const xy1 = 4
const xz1 = 5
const yz1 = 6
const yx1 = 7
const zx1 = 8
const zy1 = 9
const xy2 = 10
const xz2 = 11
const yz2 = 12
const yx2 = 13
const zx2 = 14
const zy2 = 15
const xy3 = 16
const xz3 = 17
const yz3 = 18
const yx3 = 19
const zx3 = 20
const zy3 = 21
end

struct ThreadLocalBuffersT{T}
    V12_addX::Array{T,3}
    V34_addX::Array{T,3}
    X_sum_addX::Array{T,3}
    X_sum_addY::Array{T,3}
    spropX::Array{T,3}
    spropY::Array{T,4}
    V13_addY::Vector{T}
    V24_addY::Vector{T}
    V31_addY::Vector{T}
    V42_addY::Vector{T}

end

function get_ThreadLocalBuffers(
    N,
    System,
    iuh_blocksize::Int,
    ComputeType::Type{<:AbstractFloat},
    nbuffers = Threads.nthreads(),
)::Vector{ThreadLocalBuffersT{ComputeType}}
    (; Npairs, NUnique) = System

    [
        ThreadLocalBuffersT(
            zeros(ComputeType, iuh_blocksize, 21, Npairs),
            zeros(ComputeType, iuh_blocksize, 21, Npairs),
            zeros(ComputeType, iuh_blocksize, 21, Npairs),
            zeros(ComputeType, iuh_blocksize, 21, Npairs),
            zeros(ComputeType, 3, 3, NUnique),
            zeros(ComputeType, 3, 3, NUnique, NUnique),
            zeros(ComputeType, 21),
            zeros(ComputeType, 21),
            zeros(ComputeType, 21),
            zeros(ComputeType, 21),
        ) for _ = 1:nbuffers
    ]
end

# The main bottleneck seems to me to be located in the creation of large
# arrays of size 42 and 21 and the continued calling fo the V_ function.
function addX!(
    X_sum_addX::Array{T,3},
    Gamma::Array{T,5},
    System::Geometry,
    N::Int64,
    is::Integer,
    it::Integer,
    nwpr::Integer,
    Props::Array{T,3},
    Buffers::ThreadLocalBuffersT{T},
    sitesum_split,
    iuh_start::Integer,
    block_length::Integer,
) where {T}
    (; Npairs, invpairs) = System
    ns = is - 1
    nt = it - 1
    (; V12_addX, V34_addX) = Buffers

    for iuh_local = 1:block_length
        iuh_global = iuh_start + iuh_local - 1
        iu_parity = (is + it) % 2
        iu = (iuh_global - 1) * 2 + 1 + (1 - iu_parity)

        nu = iu - 1

        wpw1, wpw2, _, _, _, _, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

        flavTransf12 = (wpw1 * wpw2 > 0, ns * wpw2 > 0, ns * wpw1 < 0)
        flavTransf34 = (wmw3 * wmw4 < 0, ns * wmw4 > 0, ns * wmw3 > 0)

        # get fields of siteSum struct as Matrices for better use of LoopVectorization
        s1, t1, u1 = ConvertFreqArgs(ns, wpw1, -wpw2, N)
        s2, t2, u2 = ConvertFreqArgs(ns, -wmw3, -wmw4, N)

        swap_R12 = flavTransf12[1]
        swap_R34 = flavTransf34[1]


        @inbounds for ki = 1:Npairs

            R12 = swap_R12 ? invpairs[ki] : ki
            R34 = swap_R34 ? invpairs[ki] : ki
            Vert!((@view V12_addX[iuh_local, :, ki]), Gamma, s1, t1, u1, flavTransf12, R12)
            Vert!((@view V34_addX[iuh_local, :, ki]), Gamma, s2, t2, u2, flavTransf34, R34)
        end
    end

    fill!((@view X_sum_addX[1:block_length, :, :]), zero(T))

    function loop_sitesum_split(f, sitesum_split)
        for (; S, Nsum_split, Rijmin, Rijmax, kimin, kimax) in sitesum_split,
            Rij = Rijmin:Rijmax,
            ki = kimin:kimax

            Rij_iblock = Rij - Rijmin + 1
            ki_iblock = ki - kimin + 1

            for k_spl = 1:Nsum_split[ki_iblock, Rij_iblock]
                (; kj, m, xk) = S[k_spl, ki_iblock, Rij_iblock]

                f(Rij, ki, kj, m, xk)
            end
        end
    end

    function loop_sitesum_split_v2_blocked(f, sitesum_split)
        unrolled_kj_m_xk, cumsum_ns = sitesum_split
        @inline @inbounds start(Rij, ki) =
            let i = ki + (Rij - 1) * Npairs
                if i > 1
                    cumsum_ns[i-1] + 1
                else
                    1
                end
            end
        @inline @inbounds stop(Rij, ki) = cumsum_ns[ki+(Rij-1)*Npairs]

        nblocks = Int(ceil(Npairs / 15)) # max_blocksize=15
        blocksize = Npairs / nblocks

        @inline blockstart(iblock) = round(Int8, blocksize * (iblock - 1)) + 1
        @inline blockstop(iblock) = round(Int8, min(blocksize * iblock, Npairs))

        for Rijblock = 1:nblocks,
            kiblock = 1:nblocks,
            Rij = blockstart(Rijblock):blockstop(Rijblock),
            ki = blockstart(kiblock):blockstop(kiblock)

            for k_spl = start(Rij, ki):stop(Rij, ki)
                (; kj, m, xk) = unrolled_kj_m_xk[k_spl]

                f(Rij, ki, kj, m, xk)
            end
        end
    end

    function loop_sitesum_split_v2_NOblocked(f, sitesum_split)
        unrolled_kj_m_xk, cumsum_ns = sitesum_split
        @inline @inbounds start(Rij, ki) =
            let i = ki + (Rij - 1) * Npairs
                if i > 1
                    cumsum_ns[i-1] + 1
                else
                    1
                end
            end
        @inline @inbounds stop(Rij, ki) = cumsum_ns[ki+(Rij-1)*Npairs]

        nblocks = Int(ceil(Npairs / 15)) # max_blocksize=15
        blocksize = Npairs / nblocks

        @inline blockstart(iblock) = round(Int8, blocksize * (iblock - 1)) + 1
        @inline blockstop(iblock) = round(Int8, min(blocksize * iblock, Npairs))


        for Rij = 1:Npairs, ki = 1:Npairs
            for k_spl = start(Rij, ki):stop(Rij, ki)
                (; kj, m, xk) = unrolled_kj_m_xk[k_spl]

                f(Rij, ki, kj, m, xk)

            end
        end

    end


    @inbounds @muladd begin
        ###
        ###

        #for Rij in 1:Npairs, ki in 1:Npairs
        #    for k_spl = start(Rij, ki):stop(Rij, ki)
        #        (; kj, m, xk) = unrolled_kj_m_xk[k_spl]
        ####

        #loop_sitesum_split_v2_NOblocked(sitesum_split) do Rij, ki, kj, m, xk
        loop_sitesum_split(sitesum_split) do Rij, ki, kj, m, xk
            Ptm = @SMatrix [m * Props[i, j, xk] for i = 1:3, j = 1:3]

            @.. @inbounds @fastmath begin
                X_sum_addX[1:block_length, flidx.yy, Rij] +=
                    -V12_addX[1:block_length, flidx.yy, ki] *
                    V34_addX[1:block_length, flidx.yy, kj] *
                    Ptm[2, 2] -
                    V12_addX[1:block_length, flidx.yz1, ki] *
                    V34_addX[1:block_length, flidx.zy1, kj] *
                    Ptm[3, 3] -
                    V12_addX[1:block_length, flidx.yx1, ki] *
                    V34_addX[1:block_length, flidx.xy1, kj] *
                    Ptm[1, 1]
                X_sum_addX[1:block_length, flidx.zz, Rij] +=
                    -V12_addX[1:block_length, flidx.zz, ki] *
                    V34_addX[1:block_length, flidx.zz, kj] *
                    Ptm[3, 3] -
                    V12_addX[1:block_length, flidx.zx1, ki] *
                    V34_addX[1:block_length, flidx.xz1, kj] *
                    Ptm[1, 1] -
                    V12_addX[1:block_length, flidx.zy1, ki] *
                    V34_addX[1:block_length, flidx.yz1, kj] *
                    Ptm[2, 2]
                X_sum_addX[1:block_length, flidx.xx, Rij] +=
                    -V12_addX[1:block_length, flidx.xx, ki] *
                    V34_addX[1:block_length, flidx.xx, kj] *
                    Ptm[1, 1] -
                    V12_addX[1:block_length, flidx.xy1, ki] *
                    V34_addX[1:block_length, flidx.yx1, kj] *
                    Ptm[2, 2] -
                    V12_addX[1:block_length, flidx.xz1, ki] *
                    V34_addX[1:block_length, flidx.zx1, kj] *
                    Ptm[3, 3]

                ### Xab1 += -Vaa Vab1 - Vab1 Vbb - Vac1 Vcb1
                X_sum_addX[1:block_length, flidx.xy1, Rij] +=
                    -V12_addX[1:block_length, flidx.xx, ki] *
                    V34_addX[1:block_length, flidx.xy1, kj] *
                    Ptm[1, 1] -
                    V12_addX[1:block_length, flidx.xy1, ki] *
                    V34_addX[1:block_length, flidx.yy, kj] *
                    Ptm[2, 2] -
                    V12_addX[1:block_length, flidx.xz1, ki] *
                    V34_addX[1:block_length, flidx.zy1, kj] *
                    Ptm[3, 3]
                X_sum_addX[1:block_length, flidx.xz1, Rij] +=
                    -V12_addX[1:block_length, flidx.xx, ki] *
                    V34_addX[1:block_length, flidx.xz1, kj] *
                    Ptm[1, 1] -
                    V12_addX[1:block_length, flidx.xz1, ki] *
                    V34_addX[1:block_length, flidx.zz, kj] *
                    Ptm[3, 3] -
                    V12_addX[1:block_length, flidx.xy1, ki] *
                    V34_addX[1:block_length, flidx.yz1, kj] *
                    Ptm[2, 2]
                X_sum_addX[1:block_length, flidx.yx1, Rij] +=
                    -V12_addX[1:block_length, flidx.yy, ki] *
                    V34_addX[1:block_length, flidx.yx1, kj] *
                    Ptm[2, 2] -
                    V12_addX[1:block_length, flidx.yx1, ki] *
                    V34_addX[1:block_length, flidx.xx, kj] *
                    Ptm[1, 1] -
                    V12_addX[1:block_length, flidx.yz1, ki] *
                    V34_addX[1:block_length, flidx.zx1, kj] *
                    Ptm[3, 3]
                X_sum_addX[1:block_length, flidx.yz1, Rij] +=
                    -V12_addX[1:block_length, flidx.yy, ki] *
                    V34_addX[1:block_length, flidx.yz1, kj] *
                    Ptm[2, 2] -
                    V12_addX[1:block_length, flidx.yz1, ki] *
                    V34_addX[1:block_length, flidx.zz, kj] *
                    Ptm[3, 3] -
                    V12_addX[1:block_length, flidx.yx1, ki] *
                    V34_addX[1:block_length, flidx.xz1, kj] *
                    Ptm[1, 1]
                X_sum_addX[1:block_length, flidx.zx1, Rij] +=
                    -V12_addX[1:block_length, flidx.zz, ki] *
                    V34_addX[1:block_length, flidx.zx1, kj] *
                    Ptm[3, 3] -
                    V12_addX[1:block_length, flidx.zx1, ki] *
                    V34_addX[1:block_length, flidx.xx, kj] *
                    Ptm[1, 1] -
                    V12_addX[1:block_length, flidx.zy1, ki] *
                    V34_addX[1:block_length, flidx.yx1, kj] *
                    Ptm[2, 2]
                X_sum_addX[1:block_length, flidx.zy1, Rij] +=
                    -V12_addX[1:block_length, flidx.zz, ki] *
                    V34_addX[1:block_length, flidx.zy1, kj] *
                    Ptm[3, 3] -
                    V12_addX[1:block_length, flidx.zy1, ki] *
                    V34_addX[1:block_length, flidx.yy, kj] *
                    Ptm[2, 2] -
                    V12_addX[1:block_length, flidx.zx1, ki] *
                    V34_addX[1:block_length, flidx.xy1, kj] *
                    Ptm[1, 1]

                ### Xab2 += -Vab2 Vab2 - Vab3 Vba3
                X_sum_addX[1:block_length, flidx.xy2, Rij] +=
                    -V12_addX[1:block_length, flidx.xy2, ki] *
                    V34_addX[1:block_length, flidx.xy2, kj] *
                    Ptm[1, 2] -
                    V12_addX[1:block_length, flidx.xy3, ki] *
                    V34_addX[1:block_length, flidx.yx3, kj] *
                    Ptm[2, 1]
                X_sum_addX[1:block_length, flidx.xz2, Rij] +=
                    -V12_addX[1:block_length, flidx.xz2, ki] *
                    V34_addX[1:block_length, flidx.xz2, kj] *
                    Ptm[1, 3] -
                    V12_addX[1:block_length, flidx.xz3, ki] *
                    V34_addX[1:block_length, flidx.zx3, kj] *
                    Ptm[3, 1]
                X_sum_addX[1:block_length, flidx.yx2, Rij] +=
                    -V12_addX[1:block_length, flidx.yx2, ki] *
                    V34_addX[1:block_length, flidx.yx2, kj] *
                    Ptm[2, 1] -
                    V12_addX[1:block_length, flidx.yx3, ki] *
                    V34_addX[1:block_length, flidx.xy3, kj] *
                    Ptm[1, 2]
                X_sum_addX[1:block_length, flidx.yz2, Rij] +=
                    -V12_addX[1:block_length, flidx.yz2, ki] *
                    V34_addX[1:block_length, flidx.yz2, kj] *
                    Ptm[2, 3] -
                    V12_addX[1:block_length, flidx.yz3, ki] *
                    V34_addX[1:block_length, flidx.zy3, kj] *
                    Ptm[3, 2]
                X_sum_addX[1:block_length, flidx.zx2, Rij] +=
                    -V12_addX[1:block_length, flidx.zx2, ki] *
                    V34_addX[1:block_length, flidx.zx2, kj] *
                    Ptm[3, 1] -
                    V12_addX[1:block_length, flidx.zx3, ki] *
                    V34_addX[1:block_length, flidx.xz3, kj] *
                    Ptm[1, 3]
                X_sum_addX[1:block_length, flidx.zy2, Rij] +=
                    -V12_addX[1:block_length, flidx.zy2, ki] *
                    V34_addX[1:block_length, flidx.zy2, kj] *
                    Ptm[3, 2] -
                    V12_addX[1:block_length, flidx.zy3, ki] *
                    V34_addX[1:block_length, flidx.yz3, kj] *
                    Ptm[2, 3]

                ### Xab3 += -Vab2 Vab3 - Vab3 Vba2
                X_sum_addX[1:block_length, flidx.xy3, Rij] +=
                    -V12_addX[1:block_length, flidx.xy2, ki] *
                    V34_addX[1:block_length, flidx.xy3, kj] *
                    Ptm[1, 2] -
                    V12_addX[1:block_length, flidx.xy3, ki] *
                    V34_addX[1:block_length, flidx.yx2, kj] *
                    Ptm[2, 1]
                X_sum_addX[1:block_length, flidx.xz3, Rij] +=
                    -V12_addX[1:block_length, flidx.xz2, ki] *
                    V34_addX[1:block_length, flidx.xz3, kj] *
                    Ptm[1, 3] -
                    V12_addX[1:block_length, flidx.xz3, ki] *
                    V34_addX[1:block_length, flidx.zx2, kj] *
                    Ptm[3, 1]
                X_sum_addX[1:block_length, flidx.yx3, Rij] +=
                    -V12_addX[1:block_length, flidx.yx2, ki] *
                    V34_addX[1:block_length, flidx.yx3, kj] *
                    Ptm[2, 1] -
                    V12_addX[1:block_length, flidx.yx3, ki] *
                    V34_addX[1:block_length, flidx.xy2, kj] *
                    Ptm[1, 2]
                X_sum_addX[1:block_length, flidx.yz3, Rij] +=
                    -V12_addX[1:block_length, flidx.yz2, ki] *
                    V34_addX[1:block_length, flidx.yz3, kj] *
                    Ptm[2, 3] -
                    V12_addX[1:block_length, flidx.yz3, ki] *
                    V34_addX[1:block_length, flidx.zy2, kj] *
                    Ptm[3, 2]
                X_sum_addX[1:block_length, flidx.zx3, Rij] +=
                    -V12_addX[1:block_length, flidx.zx2, ki] *
                    V34_addX[1:block_length, flidx.zx3, kj] *
                    Ptm[3, 1] -
                    V12_addX[1:block_length, flidx.zx3, ki] *
                    V34_addX[1:block_length, flidx.xz2, kj] *
                    Ptm[1, 3]
                X_sum_addX[1:block_length, flidx.zy3, Rij] +=
                    -V12_addX[1:block_length, flidx.zy2, ki] *
                    V34_addX[1:block_length, flidx.zy3, kj] *
                    Ptm[3, 2] -
                    V12_addX[1:block_length, flidx.zy3, ki] *
                    V34_addX[1:block_length, flidx.yz2, kj] *
                    Ptm[2, 3]
            end
        end
    end

    return
end

function addY!(
    X_sum_addY::Array{T,3},
    Gamma::Array{T,5},
    System::Geometry,
    N::Int64,
    is::Integer,
    it::Integer,
    nwpr::Integer,
    Props::Array{T,4},
    Buffers::ThreadLocalBuffersT{T},
    iuh_start::Integer,
    block_length::Integer,
) where {T}
    (; Npairs, invpairs, PairTypes, OnsitePairs) = System
    ns = is - 1
    nt = it - 1

    (; V13_addY, V24_addY, V31_addY, V42_addY) = Buffers
    fill!(X_sum_addY, zero(T))
    for iuh_local = 1:block_length
        iuh_global = iuh_start + iuh_local - 1
        iu_parity = (is + it) % 2
        iu = (iuh_global - 1) * 2 + 1 + (1 - iu_parity)

        nu = iu - 1
        _, wpw2, _, wpw4, wmw1, _, wmw3, _ = mixedFrequencies(ns, nt, nu, nwpr)
        flavTransf13 = (nt * wmw3 < 0, wmw1 * wmw3 > 0, wmw1 * nt > 0)
        flavTransf24 = (nt * wpw4 < 0, wpw2 * wpw4 > 0, wpw2 * nt > 0)
        flavTransf31 = (nt * wmw1 > 0, wmw3 * wmw1 > 0, wmw3 * nt < 0)
        flavTransf42 = (nt * wpw2 > 0, wpw4 * wpw2 > 0, wpw4 * nt < 0)


        V13 = V13_addY
        V24 = V24_addY
        V31 = V31_addY
        V42 = V42_addY


        s13, t13, u13 = ConvertFreqArgs(wmw1, nt, wmw3, N)
        s24, t24, u24 = ConvertFreqArgs(wpw2, -nt, -wpw4, N)
        s31, t31, u31 = ConvertFreqArgs(wmw3, nt, -wmw1, N)
        s42, t42, u42 = ConvertFreqArgs(-wpw4, -nt, wpw2, N)

        swap13 = flavTransf13[1]
        swap24 = flavTransf24[1]
        swap31 = flavTransf31[1]
        swap42 = flavTransf42[1]

        # Xtilde only defined for nonlocal pairs Rij != Rii
        @inbounds @muladd for Rij = 1:Npairs
            Rij in OnsitePairs && continue
            # loop over all left hand side inequivalent pairs Rij
            # loop over all left hand side inequivalent pairs Rij
            #Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij)
            (; xi, xj) = PairTypes[Rij]

            R13 = swap13 ? invpairs[Rij] : Rij
            R24 = swap24 ? invpairs[Rij] : Rij
            R31 = swap31 ? invpairs[Rij] : Rij
            R42 = swap42 ? invpairs[Rij] : Rij

            Vert!(V13, Gamma, s13, t13, u13, flavTransf13, R13)
            Vert!(V24, Gamma, s24, t24, u24, flavTransf24, R24)
            Vert!(V31, Gamma, s31, t31, u31, flavTransf31, R31)
            Vert!(V42, Gamma, s42, t42, u42, flavTransf42, R42)


            P = @SMatrix [Props[i, j, xi, xj] for i = 1:3, j = 1:3]
            PT = @SMatrix [Props[j, i, xj, xi] for i = 1:3, j = 1:3]

            ### Yaa = Vaa Vaa + Vab2 Vab2 + Vac2 Vac2 + (w -- -w + t)

            X_sum_addY[iuh_local, flidx.xx, Rij] =
                X_sum_addY[iuh_local, flidx.xx, Rij] + (
                    (
                        V13[flidx.xx] * V24[flidx.xx] * P[1, 1] +
                        V13[flidx.xy2] * V24[flidx.xy2] * P[2, 2] +
                        V13[flidx.xz2] * V24[flidx.xz2] * P[3, 3]
                    ) + (
                        V31[flidx.xx] * V42[flidx.xx] * PT[1, 1] +
                        V31[flidx.xy2] * V42[flidx.xy2] * PT[2, 2] +
                        V31[flidx.xz2] * V42[flidx.xz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.yy, Rij] =
                X_sum_addY[iuh_local, flidx.yy, Rij] + (
                    (
                        V13[flidx.yy] * V24[flidx.yy] * P[2, 2] +
                        V13[flidx.yx2] * V24[flidx.yx2] * P[1, 1] +
                        V13[flidx.yz2] * V24[flidx.yz2] * P[3, 3]
                    ) + (
                        V31[flidx.yy] * V42[flidx.yy] * PT[2, 2] +
                        V31[flidx.yx2] * V42[flidx.yx2] * PT[1, 1] +
                        V31[flidx.yz2] * V42[flidx.yz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.zz, Rij] =
                X_sum_addY[iuh_local, flidx.zz, Rij] + (
                    (
                        V13[flidx.zz] * V24[flidx.zz] * P[3, 3] +
                        V13[flidx.zx2] * V24[flidx.zx2] * P[1, 1] +
                        V13[flidx.zy2] * V24[flidx.zy2] * P[2, 2]
                    ) + (
                        V31[flidx.zz] * V42[flidx.zz] * PT[3, 3] +
                        V31[flidx.zx2] * V42[flidx.zx2] * PT[1, 1] +
                        V31[flidx.zy2] * V42[flidx.zy2] * PT[2, 2]
                    )
                )

            ### Yab1 = Vab3 Vab3 + Vab1 Vab1 + (w -- -w + t)

            X_sum_addY[iuh_local, flidx.xy1, Rij] =
                X_sum_addY[iuh_local, flidx.xy1, Rij] + (
                    (
                        V13[flidx.xy3] * V24[flidx.xy3] * P[2, 1] +
                        V13[flidx.xy1] * V24[flidx.xy1] * P[1, 2]
                    ) + (
                        V31[flidx.xy3] * V42[flidx.xy3] * PT[2, 1] +
                        V31[flidx.xy1] * V42[flidx.xy1] * PT[1, 2]
                    )
                )

            X_sum_addY[iuh_local, flidx.xz1, Rij] =
                X_sum_addY[iuh_local, flidx.xz1, Rij] + (
                    (
                        V13[flidx.xz3] * V24[flidx.xz3] * P[3, 1] +
                        V13[flidx.xz1] * V24[flidx.xz1] * P[1, 3]
                    ) + (
                        V31[flidx.xz3] * V42[flidx.xz3] * PT[3, 1] +
                        V31[flidx.xz1] * V42[flidx.xz1] * PT[1, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.yx1, Rij] =
                X_sum_addY[iuh_local, flidx.yx1, Rij] + (
                    (
                        V13[flidx.yx3] * V24[flidx.yx3] * P[1, 2] +
                        V13[flidx.yx1] * V24[flidx.yx1] * P[2, 1]
                    ) + (
                        V31[flidx.yx3] * V42[flidx.yx3] * PT[1, 2] +
                        V31[flidx.yx1] * V42[flidx.yx1] * PT[2, 1]
                    )
                )

            X_sum_addY[iuh_local, flidx.yz1, Rij] =
                X_sum_addY[iuh_local, flidx.yz1, Rij] + (
                    (
                        V13[flidx.yz3] * V24[flidx.yz3] * P[3, 2] +
                        V13[flidx.yz1] * V24[flidx.yz1] * P[2, 3]
                    ) + (
                        V31[flidx.yz3] * V42[flidx.yz3] * PT[3, 2] +
                        V31[flidx.yz1] * V42[flidx.yz1] * PT[2, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.zx1, Rij] =
                X_sum_addY[iuh_local, flidx.zx1, Rij] + (
                    (
                        V13[flidx.zx3] * V24[flidx.zx3] * P[1, 3] +
                        V13[flidx.zx1] * V24[flidx.zx1] * P[3, 1]
                    ) + (
                        V31[flidx.zx3] * V42[flidx.zx3] * PT[1, 3] +
                        V31[flidx.zx1] * V42[flidx.zx1] * PT[3, 1]
                    )
                )

            X_sum_addY[iuh_local, flidx.zy1, Rij] =
                X_sum_addY[iuh_local, flidx.zy1, Rij] + (
                    (
                        V13[flidx.zy3] * V24[flidx.zy3] * P[2, 3] +
                        V13[flidx.zy1] * V24[flidx.zy1] * P[3, 2]
                    ) + (
                        V31[flidx.zy3] * V42[flidx.zy3] * PT[2, 3] +
                        V31[flidx.zy1] * V42[flidx.zy1] * PT[3, 2]
                    )
                )

            ### Yab2 = Vaa Vba2 + Vab2 Vbb + Vac2 Vbc2 + (w -- -w + t)

            X_sum_addY[iuh_local, flidx.xy2, Rij] =
                X_sum_addY[iuh_local, flidx.xy2, Rij] + (
                    (
                        V13[flidx.xx] * V24[flidx.yx2] * P[1, 1] +
                        V13[flidx.xy2] * V24[flidx.yy] * P[2, 2] +
                        V13[flidx.xz2] * V24[flidx.yz2] * P[3, 3]
                    ) + (
                        V31[flidx.xx] * V42[flidx.yx2] * PT[1, 1] +
                        V31[flidx.xy2] * V42[flidx.yy] * PT[2, 2] +
                        V31[flidx.xz2] * V42[flidx.yz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.xz2, Rij] =
                X_sum_addY[iuh_local, flidx.xz2, Rij] + (
                    (
                        V13[flidx.xx] * V24[flidx.zx2] * P[1, 1] +
                        V13[flidx.xz2] * V24[flidx.zz] * P[3, 3] +
                        V13[flidx.xy2] * V24[flidx.zy2] * P[2, 2]
                    ) + (
                        V31[flidx.xx] * V42[flidx.zx2] * PT[1, 1] +
                        V31[flidx.xz2] * V42[flidx.zz] * PT[3, 3] +
                        V31[flidx.xy2] * V42[flidx.zy2] * PT[2, 2]
                    )
                )

            X_sum_addY[iuh_local, flidx.yx2, Rij] =
                X_sum_addY[iuh_local, flidx.yx2, Rij] + (
                    (
                        V13[flidx.yy] * V24[flidx.xy2] * P[2, 2] +
                        V13[flidx.yx2] * V24[flidx.xx] * P[1, 1] +
                        V13[flidx.yz2] * V24[flidx.xz2] * P[3, 3]
                    ) + (
                        V31[flidx.yy] * V42[flidx.xy2] * PT[2, 2] +
                        V31[flidx.yx2] * V42[flidx.xx] * PT[1, 1] +
                        V31[flidx.yz2] * V42[flidx.xz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.yz2, Rij] =
                X_sum_addY[iuh_local, flidx.yz2, Rij] + (
                    (
                        V13[flidx.yy] * V24[flidx.zy2] * P[2, 2] +
                        V13[flidx.yz2] * V24[flidx.zz] * P[3, 3] +
                        V13[flidx.yx2] * V24[flidx.zx2] * P[1, 1]
                    ) + (
                        V31[flidx.yy] * V42[flidx.zy2] * PT[2, 2] +
                        V31[flidx.yz2] * V42[flidx.zz] * PT[3, 3] +
                        V31[flidx.yx2] * V42[flidx.zx2] * PT[1, 1]
                    )
                )

            X_sum_addY[iuh_local, flidx.zx2, Rij] =
                X_sum_addY[iuh_local, flidx.zx2, Rij] + (
                    (
                        V13[flidx.zz] * V24[flidx.xz2] * P[3, 3] +
                        V13[flidx.zx2] * V24[flidx.xx] * P[1, 1] +
                        V13[flidx.zy2] * V24[flidx.xy2] * P[2, 2]
                    ) + (
                        V31[flidx.zz] * V42[flidx.xz2] * PT[3, 3] +
                        V31[flidx.zx2] * V42[flidx.xx] * PT[1, 1] +
                        V31[flidx.zy2] * V42[flidx.xy2] * PT[2, 2]
                    )
                )

            X_sum_addY[iuh_local, flidx.zy2, Rij] =
                X_sum_addY[iuh_local, flidx.zy2, Rij] + (
                    (
                        V13[flidx.zz] * V24[flidx.yz2] * P[3, 3] +
                        V13[flidx.zy2] * V24[flidx.yy] * P[2, 2] +
                        V13[flidx.zx2] * V24[flidx.yx2] * P[1, 1]
                    ) + (
                        V31[flidx.zz] * V42[flidx.yz2] * PT[3, 3] +
                        V31[flidx.zy2] * V42[flidx.yy] * PT[2, 2] +
                        V31[flidx.zx2] * V42[flidx.yx2] * PT[1, 1]
                    )
                )

            ### Yab3 = Vab3 Vba1 + Vab1 Vba3 + (w -- -w + t)

            X_sum_addY[iuh_local, flidx.xy3, Rij] =
                X_sum_addY[iuh_local, flidx.xy3, Rij] + (
                    (
                        V13[flidx.xy3] * V24[flidx.yx1] * P[2, 1] +
                        V13[flidx.xy1] * V24[flidx.yx3] * P[1, 2]
                    ) + (
                        V31[flidx.xy3] * V42[flidx.yx1] * PT[2, 1] +
                        V31[flidx.xy1] * V42[flidx.yx3] * PT[1, 2]
                    )
                )

            X_sum_addY[iuh_local, flidx.xz3, Rij] =
                X_sum_addY[iuh_local, flidx.xz3, Rij] + (
                    (
                        V13[flidx.xz3] * V24[flidx.zx1] * P[3, 1] +
                        V13[flidx.xz1] * V24[flidx.zx3] * P[1, 3]
                    ) + (
                        V31[flidx.xz3] * V42[flidx.zx1] * PT[3, 1] +
                        V31[flidx.xz1] * V42[flidx.zx3] * PT[1, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.yx3, Rij] =
                X_sum_addY[iuh_local, flidx.yx3, Rij] + (
                    (
                        V13[flidx.yx3] * V24[flidx.xy1] * P[1, 2] +
                        V13[flidx.yx1] * V24[flidx.xy3] * P[2, 1]
                    ) + (
                        V31[flidx.yx3] * V42[flidx.xy1] * PT[1, 2] +
                        V31[flidx.yx1] * V42[flidx.xy3] * PT[2, 1]
                    )
                )

            X_sum_addY[iuh_local, flidx.yz3, Rij] =
                X_sum_addY[iuh_local, flidx.yz3, Rij] + (
                    (
                        V13[flidx.yz3] * V24[flidx.zy1] * P[3, 2] +
                        V13[flidx.yz1] * V24[flidx.zy3] * P[2, 3]
                    ) + (
                        V31[flidx.yz3] * V42[flidx.zy1] * PT[3, 2] +
                        V31[flidx.yz1] * V42[flidx.zy3] * PT[2, 3]
                    )
                )

            X_sum_addY[iuh_local, flidx.zx3, Rij] =
                X_sum_addY[iuh_local, flidx.zx3, Rij] + (
                    (
                        V13[flidx.zx3] * V24[flidx.xz1] * P[1, 3] +
                        V13[flidx.zx1] * V24[flidx.xz3] * P[3, 1]
                    ) + (
                        V31[flidx.zx3] * V42[flidx.xz1] * PT[1, 3] +
                        V31[flidx.zx1] * V42[flidx.xz3] * PT[3, 1]
                    )
                )

            X_sum_addY[iuh_local, flidx.zy3, Rij] =
                X_sum_addY[iuh_local, flidx.zy3, Rij] + (
                    (
                        V13[flidx.zy3] * V24[flidx.yz1] * P[2, 3] +
                        V13[flidx.zy1] * V24[flidx.yz3] * P[3, 2]
                    ) + (
                        V31[flidx.zy3] * V42[flidx.yz1] * PT[2, 3] +
                        V31[flidx.zy1] * V42[flidx.yz3] * PT[3, 2]
                    )
                )

        end
    end
end


function Xh_from_X(X::Array{T,5}) where {T}
    Nflavours, NUnique, N, _, _ = size(X)

    zeros(Int(N / 2), Nflavours, NUnique, N, N)
end

function X_from_Xh!(X::Array{T,5}, Xh::Array{T,5}) where {T}
    Nflavours, NUnique, N, _, _ = size(X)
    for is = 1:N, it = 1:N, Rij = 1:NUnique, n = 1:Nflavours, iu = 1:N
        if (is + it + iu) % 2 == 0
            X[n, Rij, is, it, iu] = Xh[(iu-1)÷2+1, n, Rij, it, is]
        else
            X[n, Rij, is, it, iu] = zero(T)
        end
    end
end

"""
NOTE: This functions needs 120k allocations for lattice_size=16, square lattice.
      which is quite weird.
"""
function split_sitesum_v2_noblocking(siteSum, Npairs, Nsum)
    unrolled_kj_m_xk = Vector{@NamedTuple{kj::Int8, m::Int8, xk::Int8}}(undef, sum(Nsum))

    counter = 1
    for Rij = 1:Npairs, ki = 1:Npairs, x in (@view siteSum[:, Rij])
        if x.ki == ki
            unrolled_kj_m_xk[counter] = (kj = x.kj, m = x.m, xk = x.xk)
            counter += 1
        end
    end

    ns = Vector{Int8}(undef, (Npairs * Npairs))

    counter = 1
    for Rij = 1:Npairs, ki = 1:Npairs
        ns[counter] = count(x -> x.ki == ki, (@view siteSum[:, Rij]))
        counter += 1
    end


    return unrolled_kj_m_xk, cumsum(ns)


end


function split_sitesum(siteSum, max_blocksize, Npairs, Nsum)
    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m


    nblocks = Int(ceil(Npairs / max_blocksize))
    blocksize = Npairs / nblocks

    limits(iblock) = (
        round(Int32, blocksize * (iblock - 1)) + 1,
        round(Int32, min(blocksize * iblock, Npairs)),
    )

    get_S_el() = Array{@NamedTuple{kj::Int8, m::Int8, xk::Int8},3}(
        undef,
        (maximum(Nsum), max_blocksize, max_blocksize),
    )

    S = Vector{
        @NamedTuple{
            S::Base.return_types(get_S_el)[1],
            Nsum_split::Matrix{Int16},
            Rijmin::Int8,
            Rijmax::Int8,
            kimin::Int8,
            kimax::Int8,
        }
    }()
    for Rijblock = 1:nblocks
        Rijmin, Rijmax = limits(Rijblock)
        for kiblock = 1:nblocks
            kimin, kimax = limits(kiblock)
            S_el = get_S_el()
            for i in eachindex(S_el)
                S_el[i] = (kj = 0, m = 0, xk = 0)
            end

            Nsum_el = zeros(Int32, (kimax - kimin + 1, Rijmax - Rijmin + 1))

            for Rij = Rijmin:Rijmax
                for ki_target = kimin:kimax
                    counter = 0
                    for k_spl_old = 1:Nsum[Rij]
                        ki_old = S_ki[k_spl_old, Rij]
                        kj_old = S_kj[k_spl_old, Rij]
                        xk_old = S_xk[k_spl_old, Rij]
                        m_old = S_m[k_spl_old, Rij]

                        if ki_target == ki_old
                            counter += 1
                            S_el[counter, ki_target-kimin+1, Rij-Rijmin+1] =
                                (kj = kj_old, m = m_old, xk = xk_old)
                        end
                    end
                    Nsum_el[ki_target-kimin+1, Rij-Rijmin+1] = counter
                end

            end
            push!(
                S,
                (
                    S = S_el[1:maximum(Nsum_el), :, :],
                    Nsum_split = Nsum_el,
                    Rijmin = Rijmin,
                    Rijmax = Rijmax,
                    kimin = kimin,
                    kimax = kimax,
                ),
            )

        end
    end
    S
end



function getXBubble!(
    Workspace::OneLoopWorkspace,
    T::Real;
    ComputeType::Type{<:AbstractFloat} = Float64,
)
    Par = Workspace.Par
    (; N, lenIntw) = Par.NumericalParams
    (; NUnique, siteSum, Npairs, Nsum) = Par.System

    iSigma = Workspace.State.iSigma
    DiSigma = Workspace.Deriv.iSigma

    iG = SVector{3}([
        (x, nw) -> iG_(iSigma_i, x, nw, T) for iSigma_i in (iSigma.x, iSigma.y, iSigma.z)
    ])
    iSKat = SVector{3}([
        (x, nw) -> iSKat_(iSigma_i, DiSigma_i, x, nw, T) for (iSigma_i, DiSigma_i) in
        zip((iSigma.x, iSigma.y, iSigma.z), (DiSigma.x, DiSigma.y, DiSigma.z))
    ])

    # Convert Gamma to ComputeType if needed
    Gamma =
        eltype(Workspace.State.Gamma) == ComputeType ? Workspace.State.Gamma :
        ComputeType.(Workspace.State.Gamma)

    # Determine block size based on precision
    iuh_blocksize = ComputeType == Float64 ? 4 : 8

    ThreadLocalBuffers = get_ThreadLocalBuffers(N, Par.System, iuh_blocksize, ComputeType)

    #sitesum_split = split_sitesum_v2_noblocking(siteSum, Npairs,Nsum)
    sitesum_split = split_sitesum(siteSum, 15, Npairs, Nsum)

    # Calculate number of blocks (ceiling division)
    iuh_max = div(N, 2)
    num_iublocks = cld(iuh_max, iuh_blocksize)


    # FIXME Threads.@threads :static for is_it = 1:N*N
    for is_it = 1:N*N
        @inbounds begin
            is = (is_it - 1) ÷ N + 1
            it = (is_it - 1) % N + 1
            # WARNING:
            # This works only with :static
            Buffs = ThreadLocalBuffers[Threads.threadid()]
            ns = is - 1
            nt = it - 1

            # Precompute Katanin propagators (convert to ComputeType)
            for Rij = 1:NUnique
                for j = 1:3, i = 1:3
                    Buffs.spropX[i, j, Rij] = -iSKat[i](Rij, 0) * iG[j](Rij, 0)  # nw will be updated later
                end
            end

            for nw = -lenIntw:lenIntw-1 # Matsubara sum

                nw_ns = nw + ns
                nw_nt = nw - nt
                # Update Katanin propagators for current nw (convert to ComputeType)
                for Rij = 1:NUnique
                    for j = 1:3, i = 1:3
                        Buffs.spropX[i, j, Rij] = -iSKat[i](Rij, nw) * iG[j](Rij, nw_ns)
                    end
                end

                for Rij1 = 1:NUnique, Rij2 = 1:NUnique
                    for j = 1:3, i = 1:3
                        Buffs.spropY[i, j, Rij1, Rij2] =
                            -iSKat[i](Rij1, nw) * iG[j](Rij2, nw_nt)
                    end
                end

                # Loop over iuh blocks
                for iublock = 1:num_iublocks
                    iuh_start = (iublock - 1) * iuh_blocksize + 1
                    iuh_end = min(iuh_start + iuh_blocksize - 1, iuh_max)
                    block_length = iuh_end - iuh_start + 1

                    addY!(
                        Buffs.X_sum_addY,
                        Gamma,
                        Workspace.Par.System,
                        N,
                        is,
                        it,
                        nw,
                        Buffs.spropY,
                        Buffs,
                        iuh_start,
                        block_length,
                    ) # add to XTilde-type bubble functions

                    ### If no u--t symmetry, then add all the bubbles
                    ### If use u--t symmetry, then only add for nu smaller then nt (all other obtained by symmetry)
                    addX!(
                        Buffs.X_sum_addX,
                        Gamma,
                        Workspace.Par.System,
                        N,
                        is,
                        it,
                        nw,
                        Buffs.spropX,
                        Buffs,
                        sitesum_split,
                        iuh_start,
                        block_length,
                    )


                    # Copy results back to Workspace.X
                    for Rij = 1:Npairs, n = 1:21, iuh_local = 1:block_length
                        iuh_global = iuh_start + iuh_local - 1
                        iu_parity = (is + it) % 2
                        iu = (iuh_global - 1) * 2 + 1 + (1 - iu_parity)
                        if iu <= N
                            Workspace.X[n, Rij, is, it, iu] +=
                                Buffs.X_sum_addX[iuh_local, n, Rij]
                            Workspace.X[21+n, Rij, is, it, iu] +=
                                Buffs.X_sum_addY[iuh_local, n, Rij]
                        end
                    end

                end
            end
        end
    end
end

function symmetrizeBubble!(X::Array{T,5}, Par) where {T}
    N = Par.NumericalParams.N
    (; Npairs, OnsitePairs) = Par.System
    use_symmetry = Par.Options.use_symmetry
    # use the u <--> t symmetry
    if (use_symmetry)
        # for it in 1:N
        #     for iu in it+1:N, is in 1:N, Rij in 1:Npairs
        #         X.a[Rij,is,it,iu] = -X.a[Rij,is,iu,it]
        #         X.b[Rij,is,it,iu] = -X.b[Rij,is,iu,it]
        #         X.c[Rij,is,it,iu] = (
        #         + X.a[Rij,is,it,iu]+
        #         - X.b[Rij,is,it,iu]+
        #         + X.c[Rij,is,iu,it])
        #     end
        # end
    end
    #local definitions of X.Tilde vertices
    for iu = 1:N
        for it = 1:N, is = 1:N, R in OnsitePairs
            X[21+1, R, is, it, iu] = -X[1, R, it, is, iu]  ###
            X[21+2, R, is, it, iu] = -X[2, R, it, is, iu]  ### Yaa = Xaa
            X[21+3, R, is, it, iu] = -X[3, R, it, is, iu]  ###
            for n = 1:6
                X[21+3+n, R, is, it, iu] = -X[9+n, R, it, is, iu]      ### Yab1 = Xab2
                X[21+9+n, R, is, it, iu] = -X[3+n, R, it, is, iu]      ### Yab2 = Xab1
                X[21+15+n, R, is, it, iu] = -X[15+n, R, it, is, iu]    ### Yab3 = Xab3
            end
        end
    end
end

function addToVertexFromBubble!(Gamma::Array{T,5}, X::Array{T,5}) where {T}
    for iu in axes(Gamma, 5)
        for it in axes(Gamma, 4), is in axes(Gamma, 3), Rij in axes(Gamma, 2)
            for n = 1:9 ### Zaa(s,t,u) = -Yaa(s,u,t) ; Zab1(s,t,u) = -Yab1(s,u,t)
                Gamma[n, Rij, is, it, iu] += (
                    X[n, Rij, is, it, iu] + X[21+n, Rij, is, it, iu] -
                    X[21+n, Rij, is, iu, it]
                )
            end
            for n = 1:6 ### Zab2(s,t,u) = -Yab3(s,u,t) ; Zab3(s,t,u) = -Yab2(s,u,t)
                Gamma[9+n, Rij, is, it, iu] += (
                    X[9+n, Rij, is, it, iu] + X[21+9+n, Rij, is, it, iu] -
                    X[21+15+n, Rij, is, iu, it]
                )
                Gamma[15+n, Rij, is, it, iu] += (
                    X[15+n, Rij, is, it, iu] + X[21+15+n, Rij, is, it, iu] -
                    X[21+9+n, Rij, is, iu, it]
                )
            end
        end
    end
    return Gamma
end

function symmetrizeVertex!(Gamma::Array{T,5}, Par) where {T}
    N = Par.NumericalParams.N
    for iu = 1:N
        for it = 1:N, is = 1:N, R in Par.System.OnsitePairs
            for n = 1:6
                Gamma[9+n, R, is, it, iu] = -Gamma[3+n, R, it, is, iu] ### V^ii_ab2 = -V^ii_ab1
            end
        end
    end
end

######################################################################
######### FLOW EQUATIONS ## FLOW EQUATIONS ## FLOW EQUATIONS #########
######################################################################

function getDFint!(Workspace, T::Real)
    (; State, Deriv, Par) = Workspace
    (; lenIntw_acc) = Par.NumericalParams
    NUnique = Par.System.NUnique

    iSigmax(x, nw) = iSigma_(State.iSigma.x, x, nw)
    iSigmay(x, nw) = iSigma_(State.iSigma.y, x, nw)
    iSigmaz(x, nw) = iSigma_(State.iSigma.z, x, nw)

    iGx(x, nw) = iG_(State.iSigma.x, x, nw, T)
    iGy(x, nw) = iG_(State.iSigma.y, x, nw, T)
    iGz(x, nw) = iG_(State.iSigma.z, x, nw, T)

    iSx(x, nw) = iS_(State.iSigma.x, x, nw, T)
    iSy(x, nw) = iS_(State.iSigma.y, x, nw, T)
    iSz(x, nw) = iS_(State.iSigma.z, x, nw, T)

    for x = 1:NUnique
        sumres = 0.0
        for nw = -lenIntw_acc:lenIntw_acc-1
            w = get_w(nw, T)
            sumres += iSx(x, nw) / iGx(x, nw) * iSigmax(x, nw) / w
            sumres += iSy(x, nw) / iGy(x, nw) * iSigmay(x, nw) / w
            sumres += iSz(x, nw) / iGz(x, nw) * iSigmaz(x, nw) / w
        end
        Deriv.f_int[x] = -0.5 * sumres
    end
end

function get_Self_Energy!(Workspace, T::Real)
    Par = Workspace.Par
    @inline iSx(x, nw) = iS_(Workspace.State.iSigma.x, x, nw, T) / 2
    @inline iSy(x, nw) = iS_(Workspace.State.iSigma.y, x, nw, T) / 2
    @inline iSz(x, nw) = iS_(Workspace.State.iSigma.z, x, nw, T) / 2
    compute1PartBubble!(Workspace.Deriv.iSigma, Workspace.State.Gamma, [iSx, iSy, iSz], Par)
end

function compute1PartBubble!(Dgamma::SigmaType, Gamma::Array{T,5}, Props, Par) where {T}
    invpairs = Par.System.invpairs

    setZero!(Dgamma)
    @inline Gamma_(n, Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, n, s, t, u, isFlavorTransform, Rij, invpairs[Rij], Par.NumericalParams.N)
    addTo1PartBubble!(Dgamma, Gamma_, Props, Par)
end

function addTo1PartBubble!(Dgamma::SigmaType, Gamma_::Function, Props, Par)

    (; N, lenIntw_acc) = Par.NumericalParams
    (; siteSum, Nsum, OnsitePairs) = Par.System

    Threads.@threads for iw1 = 1:N
        nw1 = iw1 - 1
        for (x, Rx) in enumerate(OnsitePairs)
            for nw = -lenIntw_acc:lenIntw_acc-1
                jsum = zeros(3)
                wpw1 = nw1 + nw + 1
                wmw1 = nw - nw1
                for k_spl = 1:Nsum[Rx]
                    (; m, ki, xk) = siteSum[k_spl, Rx]
                    flavTransform = (wmw1 * wpw1 < 0, false, false)
                    gam = @SVector [
                        Gamma_(n, ki, 0, -wmw1, -wpw1, flavTransform) for n = 1:21
                    ]
                    jsum[flidx.xx] +=
                        (
                            gam[flidx.xx] * Props[1](xk, nw) +
                            gam[flidx.yx1] * Props[2](xk, nw) +
                            gam[flidx.zx1] * Props[3](xk, nw)
                        ) * m
                    jsum[flidx.yy] +=
                        (
                            gam[flidx.xy1] * Props[1](xk, nw) +
                            gam[flidx.yy] * Props[2](xk, nw) +
                            gam[flidx.zy1] * Props[3](xk, nw)
                        ) * m
                    jsum[flidx.zz] +=
                        (
                            gam[flidx.xz1] * Props[1](xk, nw) +
                            gam[flidx.yz1] * Props[2](xk, nw) +
                            gam[flidx.zz] * Props[3](xk, nw)
                        ) * m
                end
                Dgamma.x[x, iw1] += -jsum[1]
                Dgamma.y[x, iw1] += -jsum[2]
                Dgamma.z[x, iw1] += -jsum[3]
            end
        end
    end
end

using JLD2
function getDeriv!(Deriv, State, setup, Lam; saveArgs = true)

    (; X, Par) = setup # use pre-allocated X and XTilde to reduce garbage collector time
    Workspace = OneLoopWorkspace(State, Deriv, X, Par)

    getDFint!(Workspace, Lam)
    get_Self_Energy!(Workspace, Lam)
    getXBubble!(Workspace, Lam)
    symmetrizeBubble!(Workspace.X, Par)
    addToVertexFromBubble!(Workspace.Deriv.Gamma, Workspace.X)
    symmetrizeVertex!(Workspace.Deriv.Gamma, Par)

    return
end

####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################

t_to_Lam(t) = exp(t)
Lam_to_t(t) = log(t)

function AllocateSetup(Par::OneLoopParams)
    println("Allocate Setup")
    ## Allocate Memory:
    floattype = _getFloatType(Par)
    return (X = zeros(floattype, getBubbleVDims(Par)), Par = Par)
end

function InitializeState(Par, isotropy)

    N = Par.NumericalParams.N
    (; couplings, NUnique) = Par.System

    VDims = getVDims(Par)
    #floattype = _getFloatType(Par)

    State = ArrayPartition(
        zeros(NUnique),          ### f_int
        zeros(NUnique, N),       ### iSigma_x
        zeros(NUnique, N),       ### iSigma_y
        zeros(NUnique, N),       ### iSigma_z
        zeros(VDims),            ### Gamma
    )

    Gamma = State.x[5]

    setToBareVertex!(Gamma, couplings, isotropy)

    return State

end

function gettMesh(T_min, T_max, npoints)
    t_min = get_t_min(T_min)
    t_max = Lam_to_t(T_max)
    return LinRange(t_min, t_max, npoints)
end

function launchPMFRG!(
    State,
    setup,
    Deriv!::Function;
    method = DP5(),
    npoints = 600,
    save_steps = false,
)
    println("Solving FRG")

    Par = setup[end]
    (; temp_max, temp_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(temp_max)
    tend = get_t_min(temp_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    saved_values = SavedValues(eltype(State), Observables{eltype(State)})

    function save_func(State, t, integrator)
        chi_x = getChi_x(State, t_to_Lam(t), Par)
        chi_y = getChi_y(State, t_to_Lam(t), Par)
        chi_z = getChi_z(State, t_to_Lam(t), Par)

        return Observables(copy(chi_x), copy(chi_y), copy(chi_z))
    end

    ObsSaveat = gettMesh(temp_min, temp_max, npoints)
    saveCB = SavingCallback(
        save_func,
        saved_values,
        save_everystep = false,
        saveat = ObsSaveat,
        tdir = -1,
    )

    problem = ODEProblem(Deriv_subst!, State, (t0, tend), setup) # function, initial state, timespan, ??
    sol = solve(
        problem,
        method,
        reltol = accuracy,
        abstol = accuracy,
        save_everystep = save_steps,
        callback = saveCB,
        dt = Lam_to_t(0.2 * temp_max),
    )

    return sol, saved_values
end

function testPMFRG!(State, setup, Deriv!::Function; loadArgs = false)
    Par = setup[end]
    (; temp_max, temp_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(temp_max)
    tend = get_t_min(temp_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    der = copy(State)
    setZero!(der)

    Deriv_subst!(der, State, setup, t0, s = false)
end

SolveFRG(Par, isotropy; kwargs...) =
    launchPMFRG!(InitializeState(Par, isotropy), AllocateSetup(Par), getDeriv!; kwargs...)
TestFRG(Par, isotropy; kwargs...) =
    testPMFRG!(InitializeState(Par, isotropy), AllocateSetup(Par), getDeriv!; kwargs...)

function get_t_min(Lam)
    Lam < exp(-30) && @warn "temp_min too small! Set to exp(-30) instead."
    max(Lam_to_t(Lam), -30.0)
end

function generateSubstituteDeriv(getDeriv!::Function)

    function DerivSubs!(Deriv, State, setup, t; s = true)
        Lam = t_to_Lam(t)
        a = getDeriv!(Deriv, State, setup, Lam, saveArgs = s)
        Deriv .*= Lam
        a
    end

end

function setToBareVertex!(
    Gamma::AbstractArray{T,5},
    couplings::AbstractVector,
    isotropy::Array{T,2},
) where {T}
    for Rj in axes(Gamma, 2)
        Gamma[flidx.yz2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 1]
        Gamma[flidx.zy2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 1]
        Gamma[flidx.zx2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 2]
        Gamma[flidx.xz2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 2]
        Gamma[flidx.xy2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 3]
        Gamma[flidx.yx2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 3]

        Gamma[flidx.yz3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 1]
        Gamma[flidx.zy3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 1]
        Gamma[flidx.zx3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 2]
        Gamma[flidx.xz3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 2]
        Gamma[flidx.xy3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 3]
        Gamma[flidx.yx3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 3]
    end

    return Gamma
end

#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################

getChi_z(State::ArrayPartition, T::Real, Par) =
    getChi_z(State.x[2], State.x[3], State.x[5], T, Par)
getChi_x(State::ArrayPartition, T::Real, Par) =
    getChi_x(State.x[3], State.x[4], State.x[5], T, Par)
getChi_y(State::ArrayPartition, T::Real, Par) =
    getChi_y(State.x[4], State.x[2], State.x[5], T, Par)

function getChi_z(
    iSigmaX::AbstractArray,
    iSigmaY::AbstractArray,
    Gamma::AbstractArray,
    T::Real,
    Par,
)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGx(x, w) = iG_(iSigmaX, x, w, T)
    iGy(x, w) = iG_(iSigmaY, x, w, T)
    Vxy2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, flidx.xy2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += iGx(xi, nK) * iGy(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGx(xi, nK)^2 * iGy(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += GGGG * Vxy2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end

function getChi_x(
    iSigmaY::AbstractArray,
    iSigmaZ::AbstractArray,
    Gamma::AbstractArray,
    T::Real,
    Par,
)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGy(x, w) = iG_(iSigmaY, x, w, T)
    iGz(x, w) = iG_(iSigmaZ, x, w, T)
    Vyz2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, flidx.yz2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += iGy(xi, nK) * iGz(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGy(xi, nK)^2 * iGz(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += GGGG * Vyz2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end

function getChi_y(
    iSigmaZ::AbstractArray,
    iSigmaX::AbstractArray,
    Gamma::AbstractArray,
    T::Real,
    Par,
)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGz(x, w) = iG_(iSigmaZ, x, w, T)
    iGx(x, w) = iG_(iSigmaX, x, w, T)
    Vzx2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, flidx.zx2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += iGz(xi, nK) * iGx(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGz(xi, nK)^2 * iGx(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += GGGG * Vzx2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end

export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
