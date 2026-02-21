module PMFRG_xyz

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

using SpinFRGLattices, OrdinaryDiffEqLowOrderRK, DiffEqCallbacks, StructArrays
using Unroll
using MuladdMacro
using FastBroadcast

# Include common code shared with pmfrg_xyz (Lambda flow)
include("common.jl")

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

# Similar to Gamma I give X an extra dimension as opposed to create
# A BubbleType struct for it
# For a general Vertex there can be 3^4 = 81 flavor combinations
# In the XYZ model the SO(3) symmetry breaks down to a residual Klein-4 Symmetry
# This means that the Vertex function can only depend on two distinct flavors at most
# This gives 21 different Vertex functions.
# In my convention I dont use X and ̃X but just one big array called X.
# If I need to acces the ̃X part (which in my convention I name Y) I just go X[21 + flavor]
_getFloatType(Par) = typeof(Par.NumericalParams.accuracy)

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


#############################################################
######### PROPAGATORS ## PROPAGATORS ## PROPAGATORS #########
#############################################################

### Propagators will depend on an additional flavor
### Instead of modifying the propagators, I will simply use them
### as V_, by doing iG_(iSigma.x, ...)

function get_w(nw)
    return pi * (2 * nw + 1)
end

function iG_(iSigma::AbstractArray, x::Integer, nw::Integer, T::Real)
    w = get_w(nw)
    return 1.0 / (w * sqrt(T) + iSigma_(iSigma, x, nw))
end

### by differentiating the above inverse by T
function iS_(iSigma::AbstractArray, x::Integer, nw::Integer, T::Real)
    w = get_w(nw)
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
    w = get_w(nw)
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
using LinearAlgebra
using SparseArrays


"Optimized, in-place version of V_ to be used in addX! and addY!"
@inline function FillVBuffer!(
    V::AbstractVector,
    GammaFlavorRow::AbstractVector,
    flavTransf::Tuple{Bool,Bool,Bool},
)
    @unroll for i = 1:3
        V[i] = GammaFlavorRow[i]
    end

    if flavTransf[1]
        @unroll for i = 4:6
            V[i] = GammaFlavorRow[i+3]
        end
        @unroll for i = 7:9
            V[i] = GammaFlavorRow[i-3]
        end
    else
        @unroll for i = 4:9
            V[i] = GammaFlavorRow[i]
        end
    end

    if flavTransf[2]
        @unroll for i = 10:12
            V[i] = GammaFlavorRow[i+3]
        end
        @unroll for i = 13:15
            V[i] = GammaFlavorRow[i-3]
        end
    else
        @unroll for i = 10:15
            V[i] = GammaFlavorRow[i]
        end
    end

    if flavTransf[3]
        @unroll for i = 16:18
            V[i] = GammaFlavorRow[i+3]
        end
        @unroll for i = 19:21
            V[i] = GammaFlavorRow[i-3]
        end
    else
        @unroll for i = 16:21
            V[i] = GammaFlavorRow[i]
        end
    end
end



struct ThreadLocalBuffersT{T}
    V12_addX::Array{T,3}
    V34_addX::Array{T,3}
    X_sum_addX::Array{T,3}
    X_sum_addY::Array{T,3}
    spropX::Array{T,3}
    spropY::Array{T,4}
    Ptm::Matrix{T}
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
            zeros(ComputeType, iuh_blocksize, 21, Npairs),# V12_addX::Array{T,3}
            zeros(ComputeType, iuh_blocksize, 21, Npairs),# V34_addX::Array{T,3}
            zeros(ComputeType, iuh_blocksize, 21, Npairs),# X_sum_addX::Array{T,3}
            zeros(ComputeType, iuh_blocksize, 21, Npairs), # X_sum_addY::Array{T,3}
            zeros(ComputeType, 3, 3, NUnique),            # spropX::Array{T,3}
            zeros(ComputeType, 3, 3, NUnique, NUnique),   # spropY::Array{T,4}
            zeros(ComputeType, 3, 3),                     # Ptm::Matrix{T}
            zeros(ComputeType, 21),                       # V13_addY::Vector{T}
            zeros(ComputeType, 21),                       # V24_addY::Vector{T}
            zeros(ComputeType, 21),                       # V31_addY::Vector{T}
            zeros(ComputeType, 21),                       # V42_addY::Vector{T}
        ) for _ = 1:nbuffers
    ]
end

# The main bottleneck seems to me to be located in the creation of large
# arrays of size 42 and 21 and the continued calling fo the V_ function.
function addX!(
    X_sum_addX::Array{T,3},
    Gamma::Array{T,5},
    System::Geometry,
    N::Integer,
    is::Integer,
    it::Integer,
    nwpr::Integer,
    Props::Array{T,3},
    Buffers::ThreadLocalBuffersT,
    iuh_start::Integer,
    iuh_block_size::Integer,
) where {T}

    (; Npairs, Nsum, siteSum, invpairs) = System
    (; V12_addX, V34_addX) = Buffers
    ns = is - 1
    nt = it - 1

    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m


    for iuh_local = 1:iuh_block_size
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
            FillVBuffer!(
                (@view V12_addX[iuh_local, :, ki]),
                (@view Gamma[:, R12, s1+1, t1+1, u1+1]),
                flavTransf12,
            )
            FillVBuffer!(
                (@view V34_addX[iuh_local, :, ki]),
                (@view Gamma[:, R34, s2+1, t2+1, u2+1]),
                flavTransf34,
            )
        end
    end


    @inbounds @muladd begin
        fill!(X_sum_addX, 0.0)
        for Rij = 1:Npairs
            #loop over all left hand side inequivalent pairs Rij
            for k_spl = 1:Nsum[Rij]
                #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
                ki, kj, m, xk =
                    S_ki[k_spl, Rij], S_kj[k_spl, Rij], S_m[k_spl, Rij], S_xk[k_spl, Rij]

                Ptm = @SMatrix [m * Props[i, j, xk] for i = 1:3, j = 1:3]

                @.. @inbounds @fastmath begin
                    X_sum_addX[1:iuh_block_size, fd.yy, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.yy, ki] *
                        V34_addX[1:iuh_block_size, fd.yy, kj] *
                        Ptm[2, 2] -
                        V12_addX[1:iuh_block_size, fd.yz1, ki] *
                        V34_addX[1:iuh_block_size, fd.zy1, kj] *
                        Ptm[3, 3] -
                        V12_addX[1:iuh_block_size, fd.yx1, ki] *
                        V34_addX[1:iuh_block_size, fd.xy1, kj] *
                        Ptm[1, 1]
                    X_sum_addX[1:iuh_block_size, fd.zz, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.zz, ki] *
                        V34_addX[1:iuh_block_size, fd.zz, kj] *
                        Ptm[3, 3] -
                        V12_addX[1:iuh_block_size, fd.zx1, ki] *
                        V34_addX[1:iuh_block_size, fd.xz1, kj] *
                        Ptm[1, 1] -
                        V12_addX[1:iuh_block_size, fd.zy1, ki] *
                        V34_addX[1:iuh_block_size, fd.yz1, kj] *
                        Ptm[2, 2]
                    X_sum_addX[1:iuh_block_size, fd.xx, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.xx, ki] *
                        V34_addX[1:iuh_block_size, fd.xx, kj] *
                        Ptm[1, 1] -
                        V12_addX[1:iuh_block_size, fd.xy1, ki] *
                        V34_addX[1:iuh_block_size, fd.yx1, kj] *
                        Ptm[2, 2] -
                        V12_addX[1:iuh_block_size, fd.xz1, ki] *
                        V34_addX[1:iuh_block_size, fd.zx1, kj] *
                        Ptm[3, 3]

                    ### Xab1 += -Vaa Vab1 - Vab1 Vbb - Vac1 Vcb1
                    X_sum_addX[1:iuh_block_size, fd.xy1, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.xx, ki] *
                        V34_addX[1:iuh_block_size, fd.xy1, kj] *
                        Ptm[1, 1] -
                        V12_addX[1:iuh_block_size, fd.xy1, ki] *
                        V34_addX[1:iuh_block_size, fd.yy, kj] *
                        Ptm[2, 2] -
                        V12_addX[1:iuh_block_size, fd.xz1, ki] *
                        V34_addX[1:iuh_block_size, fd.zy1, kj] *
                        Ptm[3, 3]
                    X_sum_addX[1:iuh_block_size, fd.xz1, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.xx, ki] *
                        V34_addX[1:iuh_block_size, fd.xz1, kj] *
                        Ptm[1, 1] -
                        V12_addX[1:iuh_block_size, fd.xz1, ki] *
                        V34_addX[1:iuh_block_size, fd.zz, kj] *
                        Ptm[3, 3] -
                        V12_addX[1:iuh_block_size, fd.xy1, ki] *
                        V34_addX[1:iuh_block_size, fd.yz1, kj] *
                        Ptm[2, 2]
                    X_sum_addX[1:iuh_block_size, fd.yx1, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.yy, ki] *
                        V34_addX[1:iuh_block_size, fd.yx1, kj] *
                        Ptm[2, 2] -
                        V12_addX[1:iuh_block_size, fd.yx1, ki] *
                        V34_addX[1:iuh_block_size, fd.xx, kj] *
                        Ptm[1, 1] -
                        V12_addX[1:iuh_block_size, fd.yz1, ki] *
                        V34_addX[1:iuh_block_size, fd.zx1, kj] *
                        Ptm[3, 3]
                    X_sum_addX[1:iuh_block_size, fd.yz1, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.yy, ki] *
                        V34_addX[1:iuh_block_size, fd.yz1, kj] *
                        Ptm[2, 2] -
                        V12_addX[1:iuh_block_size, fd.yz1, ki] *
                        V34_addX[1:iuh_block_size, fd.zz, kj] *
                        Ptm[3, 3] -
                        V12_addX[1:iuh_block_size, fd.yx1, ki] *
                        V34_addX[1:iuh_block_size, fd.xz1, kj] *
                        Ptm[1, 1]
                    X_sum_addX[1:iuh_block_size, fd.zx1, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.zz, ki] *
                        V34_addX[1:iuh_block_size, fd.zx1, kj] *
                        Ptm[3, 3] -
                        V12_addX[1:iuh_block_size, fd.zx1, ki] *
                        V34_addX[1:iuh_block_size, fd.xx, kj] *
                        Ptm[1, 1] -
                        V12_addX[1:iuh_block_size, fd.zy1, ki] *
                        V34_addX[1:iuh_block_size, fd.yx1, kj] *
                        Ptm[2, 2]
                    X_sum_addX[1:iuh_block_size, fd.zy1, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.zz, ki] *
                        V34_addX[1:iuh_block_size, fd.zy1, kj] *
                        Ptm[3, 3] -
                        V12_addX[1:iuh_block_size, fd.zy1, ki] *
                        V34_addX[1:iuh_block_size, fd.yy, kj] *
                        Ptm[2, 2] -
                        V12_addX[1:iuh_block_size, fd.zx1, ki] *
                        V34_addX[1:iuh_block_size, fd.xy1, kj] *
                        Ptm[1, 1]

                    ### Xab2 += -Vab2 Vab2 - Vab3 Vba3
                    X_sum_addX[1:iuh_block_size, fd.xy2, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.xy2, ki] *
                        V34_addX[1:iuh_block_size, fd.xy2, kj] *
                        Ptm[1, 2] -
                        V12_addX[1:iuh_block_size, fd.xy3, ki] *
                        V34_addX[1:iuh_block_size, fd.yx3, kj] *
                        Ptm[2, 1]
                    X_sum_addX[1:iuh_block_size, fd.xz2, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.xz2, ki] *
                        V34_addX[1:iuh_block_size, fd.xz2, kj] *
                        Ptm[1, 3] -
                        V12_addX[1:iuh_block_size, fd.xz3, ki] *
                        V34_addX[1:iuh_block_size, fd.zx3, kj] *
                        Ptm[3, 1]
                    X_sum_addX[1:iuh_block_size, fd.yx2, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.yx2, ki] *
                        V34_addX[1:iuh_block_size, fd.yx2, kj] *
                        Ptm[2, 1] -
                        V12_addX[1:iuh_block_size, fd.yx3, ki] *
                        V34_addX[1:iuh_block_size, fd.xy3, kj] *
                        Ptm[1, 2]
                    X_sum_addX[1:iuh_block_size, fd.yz2, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.yz2, ki] *
                        V34_addX[1:iuh_block_size, fd.yz2, kj] *
                        Ptm[2, 3] -
                        V12_addX[1:iuh_block_size, fd.yz3, ki] *
                        V34_addX[1:iuh_block_size, fd.zy3, kj] *
                        Ptm[3, 2]
                    X_sum_addX[1:iuh_block_size, fd.zx2, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.zx2, ki] *
                        V34_addX[1:iuh_block_size, fd.zx2, kj] *
                        Ptm[3, 1] -
                        V12_addX[1:iuh_block_size, fd.zx3, ki] *
                        V34_addX[1:iuh_block_size, fd.xz3, kj] *
                        Ptm[1, 3]
                    X_sum_addX[1:iuh_block_size, fd.zy2, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.zy2, ki] *
                        V34_addX[1:iuh_block_size, fd.zy2, kj] *
                        Ptm[3, 2] -
                        V12_addX[1:iuh_block_size, fd.zy3, ki] *
                        V34_addX[1:iuh_block_size, fd.yz3, kj] *
                        Ptm[2, 3]

                    ### Xab3 += -Vab2 Vab3 - Vab3 Vba2
                    X_sum_addX[1:iuh_block_size, fd.xy3, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.xy2, ki] *
                        V34_addX[1:iuh_block_size, fd.xy3, kj] *
                        Ptm[1, 2] -
                        V12_addX[1:iuh_block_size, fd.xy3, ki] *
                        V34_addX[1:iuh_block_size, fd.yx2, kj] *
                        Ptm[2, 1]
                    X_sum_addX[1:iuh_block_size, fd.xz3, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.xz2, ki] *
                        V34_addX[1:iuh_block_size, fd.xz3, kj] *
                        Ptm[1, 3] -
                        V12_addX[1:iuh_block_size, fd.xz3, ki] *
                        V34_addX[1:iuh_block_size, fd.zx2, kj] *
                        Ptm[3, 1]
                    X_sum_addX[1:iuh_block_size, fd.yx3, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.yx2, ki] *
                        V34_addX[1:iuh_block_size, fd.yx3, kj] *
                        Ptm[2, 1] -
                        V12_addX[1:iuh_block_size, fd.yx3, ki] *
                        V34_addX[1:iuh_block_size, fd.xy2, kj] *
                        Ptm[1, 2]
                    X_sum_addX[1:iuh_block_size, fd.yz3, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.yz2, ki] *
                        V34_addX[1:iuh_block_size, fd.yz3, kj] *
                        Ptm[2, 3] -
                        V12_addX[1:iuh_block_size, fd.yz3, ki] *
                        V34_addX[1:iuh_block_size, fd.zy2, kj] *
                        Ptm[3, 2]
                    X_sum_addX[1:iuh_block_size, fd.zx3, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.zx2, ki] *
                        V34_addX[1:iuh_block_size, fd.zx3, kj] *
                        Ptm[3, 1] -
                        V12_addX[1:iuh_block_size, fd.zx3, ki] *
                        V34_addX[1:iuh_block_size, fd.xz2, kj] *
                        Ptm[1, 3]
                    X_sum_addX[1:iuh_block_size, fd.zy3, Rij] +=
                        -V12_addX[1:iuh_block_size, fd.zy2, ki] *
                        V34_addX[1:iuh_block_size, fd.zy3, kj] *
                        Ptm[3, 2] -
                        V12_addX[1:iuh_block_size, fd.zy3, ki] *
                        V34_addX[1:iuh_block_size, fd.yz2, kj] *
                        Ptm[2, 3]
                end
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

            FillVBuffer!(V13, (@view Gamma[:, R13, s13+1, t13+1, u13+1]), flavTransf13)
            FillVBuffer!(V24, (@view Gamma[:, R24, s24+1, t24+1, u24+1]), flavTransf24)
            FillVBuffer!(V31, (@view Gamma[:, R31, s31+1, t31+1, u31+1]), flavTransf31)
            FillVBuffer!(V42, (@view Gamma[:, R42, s42+1, t42+1, u42+1]), flavTransf42)


            P = @SMatrix [Props[i, j, xi, xj] for i = 1:3, j = 1:3]
            PT = @SMatrix [Props[j, i, xj, xi] for i = 1:3, j = 1:3]

            ### Yaa = Vaa Vaa + Vab2 Vab2 + Vac2 Vac2 + (w -- -w + t)

            X_sum_addY[iuh_local, fd.xx, Rij] =
                X_sum_addY[iuh_local, fd.xx, Rij] + (
                    (
                        V13[fd.xx] * V24[fd.xx] * P[1, 1] +
                        V13[fd.xy2] * V24[fd.xy2] * P[2, 2] +
                        V13[fd.xz2] * V24[fd.xz2] * P[3, 3]
                    ) + (
                        V31[fd.xx] * V42[fd.xx] * PT[1, 1] +
                        V31[fd.xy2] * V42[fd.xy2] * PT[2, 2] +
                        V31[fd.xz2] * V42[fd.xz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.yy, Rij] =
                X_sum_addY[iuh_local, fd.yy, Rij] + (
                    (
                        V13[fd.yy] * V24[fd.yy] * P[2, 2] +
                        V13[fd.yx2] * V24[fd.yx2] * P[1, 1] +
                        V13[fd.yz2] * V24[fd.yz2] * P[3, 3]
                    ) + (
                        V31[fd.yy] * V42[fd.yy] * PT[2, 2] +
                        V31[fd.yx2] * V42[fd.yx2] * PT[1, 1] +
                        V31[fd.yz2] * V42[fd.yz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.zz, Rij] =
                X_sum_addY[iuh_local, fd.zz, Rij] + (
                    (
                        V13[fd.zz] * V24[fd.zz] * P[3, 3] +
                        V13[fd.zx2] * V24[fd.zx2] * P[1, 1] +
                        V13[fd.zy2] * V24[fd.zy2] * P[2, 2]
                    ) + (
                        V31[fd.zz] * V42[fd.zz] * PT[3, 3] +
                        V31[fd.zx2] * V42[fd.zx2] * PT[1, 1] +
                        V31[fd.zy2] * V42[fd.zy2] * PT[2, 2]
                    )
                )

            ### Yab1 = Vab3 Vab3 + Vab1 Vab1 + (w -- -w + t)

            X_sum_addY[iuh_local, fd.xy1, Rij] =
                X_sum_addY[iuh_local, fd.xy1, Rij] + (
                    (
                        V13[fd.xy3] * V24[fd.xy3] * P[2, 1] +
                        V13[fd.xy1] * V24[fd.xy1] * P[1, 2]
                    ) + (
                        V31[fd.xy3] * V42[fd.xy3] * PT[2, 1] +
                        V31[fd.xy1] * V42[fd.xy1] * PT[1, 2]
                    )
                )

            X_sum_addY[iuh_local, fd.xz1, Rij] =
                X_sum_addY[iuh_local, fd.xz1, Rij] + (
                    (
                        V13[fd.xz3] * V24[fd.xz3] * P[3, 1] +
                        V13[fd.xz1] * V24[fd.xz1] * P[1, 3]
                    ) + (
                        V31[fd.xz3] * V42[fd.xz3] * PT[3, 1] +
                        V31[fd.xz1] * V42[fd.xz1] * PT[1, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.yx1, Rij] =
                X_sum_addY[iuh_local, fd.yx1, Rij] + (
                    (
                        V13[fd.yx3] * V24[fd.yx3] * P[1, 2] +
                        V13[fd.yx1] * V24[fd.yx1] * P[2, 1]
                    ) + (
                        V31[fd.yx3] * V42[fd.yx3] * PT[1, 2] +
                        V31[fd.yx1] * V42[fd.yx1] * PT[2, 1]
                    )
                )

            X_sum_addY[iuh_local, fd.yz1, Rij] =
                X_sum_addY[iuh_local, fd.yz1, Rij] + (
                    (
                        V13[fd.yz3] * V24[fd.yz3] * P[3, 2] +
                        V13[fd.yz1] * V24[fd.yz1] * P[2, 3]
                    ) + (
                        V31[fd.yz3] * V42[fd.yz3] * PT[3, 2] +
                        V31[fd.yz1] * V42[fd.yz1] * PT[2, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.zx1, Rij] =
                X_sum_addY[iuh_local, fd.zx1, Rij] + (
                    (
                        V13[fd.zx3] * V24[fd.zx3] * P[1, 3] +
                        V13[fd.zx1] * V24[fd.zx1] * P[3, 1]
                    ) + (
                        V31[fd.zx3] * V42[fd.zx3] * PT[1, 3] +
                        V31[fd.zx1] * V42[fd.zx1] * PT[3, 1]
                    )
                )

            X_sum_addY[iuh_local, fd.zy1, Rij] =
                X_sum_addY[iuh_local, fd.zy1, Rij] + (
                    (
                        V13[fd.zy3] * V24[fd.zy3] * P[2, 3] +
                        V13[fd.zy1] * V24[fd.zy1] * P[3, 2]
                    ) + (
                        V31[fd.zy3] * V42[fd.zy3] * PT[2, 3] +
                        V31[fd.zy1] * V42[fd.zy1] * PT[3, 2]
                    )
                )

            ### Yab2 = Vaa Vba2 + Vab2 Vbb + Vac2 Vbc2 + (w -- -w + t)

            X_sum_addY[iuh_local, fd.xy2, Rij] =
                X_sum_addY[iuh_local, fd.xy2, Rij] + (
                    (
                        V13[fd.xx] * V24[fd.yx2] * P[1, 1] +
                        V13[fd.xy2] * V24[fd.yy] * P[2, 2] +
                        V13[fd.xz2] * V24[fd.yz2] * P[3, 3]
                    ) + (
                        V31[fd.xx] * V42[fd.yx2] * PT[1, 1] +
                        V31[fd.xy2] * V42[fd.yy] * PT[2, 2] +
                        V31[fd.xz2] * V42[fd.yz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.xz2, Rij] =
                X_sum_addY[iuh_local, fd.xz2, Rij] + (
                    (
                        V13[fd.xx] * V24[fd.zx2] * P[1, 1] +
                        V13[fd.xz2] * V24[fd.zz] * P[3, 3] +
                        V13[fd.xy2] * V24[fd.zy2] * P[2, 2]
                    ) + (
                        V31[fd.xx] * V42[fd.zx2] * PT[1, 1] +
                        V31[fd.xz2] * V42[fd.zz] * PT[3, 3] +
                        V31[fd.xy2] * V42[fd.zy2] * PT[2, 2]
                    )
                )

            X_sum_addY[iuh_local, fd.yx2, Rij] =
                X_sum_addY[iuh_local, fd.yx2, Rij] + (
                    (
                        V13[fd.yy] * V24[fd.xy2] * P[2, 2] +
                        V13[fd.yx2] * V24[fd.xx] * P[1, 1] +
                        V13[fd.yz2] * V24[fd.xz2] * P[3, 3]
                    ) + (
                        V31[fd.yy] * V42[fd.xy2] * PT[2, 2] +
                        V31[fd.yx2] * V42[fd.xx] * PT[1, 1] +
                        V31[fd.yz2] * V42[fd.xz2] * PT[3, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.yz2, Rij] =
                X_sum_addY[iuh_local, fd.yz2, Rij] + (
                    (
                        V13[fd.yy] * V24[fd.zy2] * P[2, 2] +
                        V13[fd.yz2] * V24[fd.zz] * P[3, 3] +
                        V13[fd.yx2] * V24[fd.zx2] * P[1, 1]
                    ) + (
                        V31[fd.yy] * V42[fd.zy2] * PT[2, 2] +
                        V31[fd.yz2] * V42[fd.zz] * PT[3, 3] +
                        V31[fd.yx2] * V42[fd.zx2] * PT[1, 1]
                    )
                )

            X_sum_addY[iuh_local, fd.zx2, Rij] =
                X_sum_addY[iuh_local, fd.zx2, Rij] + (
                    (
                        V13[fd.zz] * V24[fd.xz2] * P[3, 3] +
                        V13[fd.zx2] * V24[fd.xx] * P[1, 1] +
                        V13[fd.zy2] * V24[fd.xy2] * P[2, 2]
                    ) + (
                        V31[fd.zz] * V42[fd.xz2] * PT[3, 3] +
                        V31[fd.zx2] * V42[fd.xx] * PT[1, 1] +
                        V31[fd.zy2] * V42[fd.xy2] * PT[2, 2]
                    )
                )

            X_sum_addY[iuh_local, fd.zy2, Rij] =
                X_sum_addY[iuh_local, fd.zy2, Rij] + (
                    (
                        V13[fd.zz] * V24[fd.yz2] * P[3, 3] +
                        V13[fd.zy2] * V24[fd.yy] * P[2, 2] +
                        V13[fd.zx2] * V24[fd.yx2] * P[1, 1]
                    ) + (
                        V31[fd.zz] * V42[fd.yz2] * PT[3, 3] +
                        V31[fd.zy2] * V42[fd.yy] * PT[2, 2] +
                        V31[fd.zx2] * V42[fd.yx2] * PT[1, 1]
                    )
                )

            ### Yab3 = Vab3 Vba1 + Vab1 Vba3 + (w -- -w + t)

            X_sum_addY[iuh_local, fd.xy3, Rij] =
                X_sum_addY[iuh_local, fd.xy3, Rij] + (
                    (
                        V13[fd.xy3] * V24[fd.yx1] * P[2, 1] +
                        V13[fd.xy1] * V24[fd.yx3] * P[1, 2]
                    ) + (
                        V31[fd.xy3] * V42[fd.yx1] * PT[2, 1] +
                        V31[fd.xy1] * V42[fd.yx3] * PT[1, 2]
                    )
                )

            X_sum_addY[iuh_local, fd.xz3, Rij] =
                X_sum_addY[iuh_local, fd.xz3, Rij] + (
                    (
                        V13[fd.xz3] * V24[fd.zx1] * P[3, 1] +
                        V13[fd.xz1] * V24[fd.zx3] * P[1, 3]
                    ) + (
                        V31[fd.xz3] * V42[fd.zx1] * PT[3, 1] +
                        V31[fd.xz1] * V42[fd.zx3] * PT[1, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.yx3, Rij] =
                X_sum_addY[iuh_local, fd.yx3, Rij] + (
                    (
                        V13[fd.yx3] * V24[fd.xy1] * P[1, 2] +
                        V13[fd.yx1] * V24[fd.xy3] * P[2, 1]
                    ) + (
                        V31[fd.yx3] * V42[fd.xy1] * PT[1, 2] +
                        V31[fd.yx1] * V42[fd.xy3] * PT[2, 1]
                    )
                )

            X_sum_addY[iuh_local, fd.yz3, Rij] =
                X_sum_addY[iuh_local, fd.yz3, Rij] + (
                    (
                        V13[fd.yz3] * V24[fd.zy1] * P[3, 2] +
                        V13[fd.yz1] * V24[fd.zy3] * P[2, 3]
                    ) + (
                        V31[fd.yz3] * V42[fd.zy1] * PT[3, 2] +
                        V31[fd.yz1] * V42[fd.zy3] * PT[2, 3]
                    )
                )

            X_sum_addY[iuh_local, fd.zx3, Rij] =
                X_sum_addY[iuh_local, fd.zx3, Rij] + (
                    (
                        V13[fd.zx3] * V24[fd.xz1] * P[1, 3] +
                        V13[fd.zx1] * V24[fd.xz3] * P[3, 1]
                    ) + (
                        V31[fd.zx3] * V42[fd.xz1] * PT[1, 3] +
                        V31[fd.zx1] * V42[fd.xz3] * PT[3, 1]
                    )
                )

            X_sum_addY[iuh_local, fd.zy3, Rij] =
                X_sum_addY[iuh_local, fd.zy3, Rij] + (
                    (
                        V13[fd.zy3] * V24[fd.yz1] * P[2, 3] +
                        V13[fd.zy1] * V24[fd.yz3] * P[3, 2]
                    ) + (
                        V31[fd.zy3] * V42[fd.yz1] * PT[2, 3] +
                        V31[fd.zy1] * V42[fd.yz3] * PT[3, 2]
                    )
                )

        end
    end
end


function set_spropX!(spropX, NUnique, iSigma, DiSigma, T, nw1, nw2, ComputeType)

    iG = SVector{3}([
        (x, nw) -> iG_(iSigma_i, x, nw, T) for iSigma_i in (iSigma.x, iSigma.y, iSigma.z)
    ])
    iSKat = SVector{3}([
        (x, nw) -> iSKat_(iSigma_i, DiSigma_i, x, nw, T) for (iSigma_i, DiSigma_i) in
        zip((iSigma.x, iSigma.y, iSigma.z), (DiSigma.x, DiSigma.y, DiSigma.z))
    ])

    for Rij = 1:NUnique
        for j = 1:3, i = 1:3
            spropX[i, j, Rij] = ComputeType(-iSKat[i](Rij, nw1) * iG[j](Rij, nw2))
        end
    end

end

function set_spropY!(spropY, NUnique, iSigma, DiSigma, T, nw1, nw2, ComputeType)

    iG = SVector{3}([
        (x, nw) -> iG_(iSigma_i, x, nw, T) for iSigma_i in (iSigma.x, iSigma.y, iSigma.z)
    ])
    iSKat = SVector{3}([
        (x, nw) -> iSKat_(iSigma_i, DiSigma_i, x, nw, T) for (iSigma_i, DiSigma_i) in
        zip((iSigma.x, iSigma.y, iSigma.z), (DiSigma.x, DiSigma.y, DiSigma.z))
    ])


    for Rij1 = 1:NUnique, Rij2 = 1:NUnique
        for j = 1:3, i = 1:3
            spropY[i, j, Rij1, Rij2] = ComputeType(-iSKat[i](Rij1, nw1) * iG[j](Rij2, nw2))
        end
    end


end



function getXBubble!(
    Workspace::OneLoopWorkspace,
    FlowParameter::Real;
    ComputeType::Type = Float64,
)
    T = FlowParameter
    Par = Workspace.Par
    (; N, lenIntw) = Par.NumericalParams
    (; NUnique, Npairs) = Par.System

    iSigma = Workspace.State.iSigma
    DiSigma = Workspace.Deriv.iSigma

    # Convert Gamma to ComputeType if needed
    Gamma =
        eltype(Workspace.State.Gamma) == ComputeType ? Workspace.State.Gamma :
        ComputeType.(Workspace.State.Gamma)

    # Determine block size based on precision
    # iuh_blocksize = ComputeType == Float64 ? 4 : 8 # 256b - based, avx2
    optimal_iuh_blocksize = ComputeType == Float64 ? 8 : 16 # 64B - based, cache line  

    ThreadLocalBuffers =
        get_ThreadLocalBuffers(N, Par.System, optimal_iuh_blocksize, ComputeType)

    Threads.@threads :static for is_it = 1:N*N
        @inbounds begin
            is = (is_it - 1) ÷ N + 1
            it = (is_it - 1) % N + 1
            # WARNING:
            # This works only with :static
            Buffs = ThreadLocalBuffers[Threads.threadid()]
            ns = is - 1
            nt = it - 1

            for nw = -lenIntw:lenIntw-1 # Matsubara sum
                # Update Katanin propagators for current nw (convert to ComputeType)
                set_spropX!(
                    Buffs.spropX,
                    NUnique,
                    iSigma,
                    DiSigma,
                    T,
                    nw,
                    nw + ns,
                    ComputeType,
                )
                set_spropY!(
                    Buffs.spropY,
                    NUnique,
                    iSigma,
                    DiSigma,
                    T,
                    nw,
                    nw - nt,
                    ComputeType,
                )
                #for Rij1 = 1:NUnique, Rij2 = 1:NUnique
                #    for j = 1:3, i = 1:3
                #        Buffs.spropY[i, j, Rij1, Rij2] =
                #            ComputeType(-iSKat[i](Rij1, nw) * iG[j](Rij2, nw_nt))
                #    end
                #end

                # Calculate number of blocks (ceiling division)
                iuhalf_max = div(N, 2)
                # Add some clarification...
                num_iuh_blocks = cld(iuhalf_max, optimal_iuh_blocksize)

                # Loop over iuh blocks
                for iblock = 1:num_iuh_blocks
                    iuh_start = (iblock - 1) * optimal_iuh_blocksize + 1
                    iuh_end = min(iuh_start + optimal_iuh_blocksize - 1, iuhalf_max)
                    iuh_block_size = iuh_end - iuh_start + 1

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
                        iuh_block_size,
                    )

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
                        iuh_start,
                        iuh_block_size,
                    )


                    # Copy results back to Workspace.X
                    for Rij = 1:Npairs, iuh_local = 1:iuh_block_size
                        iuh_global = iuh_start + iuh_local - 1
                        iu_parity = (is + it) % 2
                        iu = (iuh_global - 1) * 2 + 1 + (1 - iu_parity)
                        if iu <= N
                            (@view Workspace.X[1:21, Rij, is, it, iu]) .+=
                                (@view Buffs.X_sum_addX[iuh_local, :, Rij])
                        end
                    end
                    for Rij = 1:Npairs, iuh_local = 1:iuh_block_size
                        iuh_global = iuh_start + iuh_local - 1
                        iu_parity = (is + it) % 2
                        iu = (iuh_global - 1) * 2 + 1 + (1 - iu_parity)
                        if iu <= N
                            (@view Workspace.X[22:42, Rij, is, it, iu]) .+=
                                (@view Buffs.X_sum_addY[iuh_local, :, Rij])
                        end
                    end

                end
            end
        end
    end
end


######################################################################
######### FLOW EQUATIONS ## FLOW EQUATIONS ## FLOW EQUATIONS #########
######################################################################

function get_iS(FlowParam::Real, iSigma::SigmaType, _::NumericalParams)
    T = FlowParam

    @inline iSx(x, nw) = iS_(iSigma.x, x, nw, T) / 2
    @inline iSy(x, nw) = iS_(iSigma.y, x, nw, T) / 2
    @inline iSz(x, nw) = iS_(iSigma.z, x, nw, T) / 2

    return iSx, iSy, iSz
end


function get_iG(FlowParam::Real, iSigma::SigmaType, _::NumericalParams)
    T = FlowParam

    @inline iGx(x, nw) = iG_(iSigma.x, x, nw, T)
    @inline iGy(x, nw) = iG_(iSigma.y, x, nw, T)
    @inline iGz(x, nw) = iG_(iSigma.z, x, nw, T)

    return iGx, iGy, iGz
end

function get_Theta(_::Real, _::NumericalParams)

    return _ -> 1
end

function get_f_int_factor(_::NumericalParams)
    return 1
end

function get_get_w(_::NumericalParams)
    return nw -> get_w(nw)
end

function get_Dgamma_factor(_::NumericalParams)
    return 1
end


####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################


function flow_parameter_max_min(NumParams::NumericalParams)

    (; temp_max, temp_min) = NumParams
    return temp_max, temp_min
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
        V_(Gamma, fd.xy2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

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
        V_(Gamma, fd.yz2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

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
        V_(Gamma, fd.zx2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

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


## 
getChi_z_refactored(State::ArrayPartition, T::Real, Par) =
    getChi_3(State.x[2], State.x[3], State.x[5], T, fd.xy2, Par)
getChi_x_refactored(State::ArrayPartition, T::Real, Par) =
    getChi_3(State.x[3], State.x[4], State.x[5], T, fd.yz2, Par)
getChi_y_refactored(State::ArrayPartition, T::Real, Par) =
    getChi_3(State.x[4], State.x[2], State.x[5], T, fd.zx2, Par)

function getChi_3(
    iSigma1::AbstractArray,
    iSigma2::AbstractArray,
    Gamma::AbstractArray,
    T::Real,
    fd_idx,
    Par,
)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iG1(x, w) = iG_(iSigma1, x, w, T)
    iG2(x, w) = iG_(iSigma2, x, w, T)
    V12_2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, fd_idx, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += iG1(xi, nK) * iG2(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iG1(xi, nK)^2 * iG2(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += GGGG * V12_2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end


export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
