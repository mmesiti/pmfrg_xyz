module PMFRG_xyz

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

using RecursiveArrayTools
using SpinFRGLattices, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools, StructArrays
using SpinFRGLattices.StaticArrays
using Unroll
# using LoopVectorization

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
            # The layout of the fd-module is hence *very* important

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
module fd
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
    V12_addX::Matrix{T}
    V34_addX::Matrix{T}
    X_sum::Vector{T}
    spropX::Array{T,3}
    spropY::Array{T,4}
    Ptm::Matrix{T}
    V13_addY::Vector{T}
    V24_addY::Vector{T}
    V31_addY::Vector{T}
    V42_addY::Vector{T}

end

function get_ThreadLocalBuffers(System)::Vector{ThreadLocalBuffersT{Float64}}
    (; Npairs, NUnique) = System

    [
        ThreadLocalBuffersT(
            zeros(21, Npairs),
            zeros(21, Npairs),
            zeros(21),
            zeros(3, 3, NUnique),
            zeros(3, 3, NUnique, NUnique),
            zeros(3, 3),
            zeros(21),
            zeros(21),
            zeros(21),
            zeros(21),
        ) for _ = 1:Threads.nthreads()
    ]
end





# The main bottleneck seems to me to be located in the creation of large
# arrays of size 42 and 21 and the continued calling fo the V_ function.
function addX!(
    X::Array{T,5},
    Gamma::Array{T,5},
    System::Geometry,
    N::Int64,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Props::Array{T,3},
    Buffers::ThreadLocalBuffersT,
) where {T}
    (; Npairs, Nsum, siteSum, invpairs) = System

    Vert(n, Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, n, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    ns = is - 1
    nt = it - 1
    nu = iu - 1

    wpw1, wpw2, _, _, _, _, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    flavTransf12 = (wpw1 * wpw2 > 0, ns * wpw2 > 0, ns * wpw1 < 0)
    flavTransf34 = (wmw3 * wmw4 < 0, ns * wmw4 > 0, ns * wmw3 > 0)

    # get fields of siteSum struct as Matrices for better use of LoopVectorization
    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    (; V12_addX, V34_addX, X_sum, Ptm) = Buffers
    V12 = V12_addX
    V34 = V34_addX
    @inbounds for ki = 1:Npairs
        for n = 1:21
            V12[n, ki] = Vert(n, ki, ns, wpw1, -wpw2, flavTransf12)
            V34[n, ki] = Vert(n, ki, ns, -wmw3, -wmw4, flavTransf34)
        end
    end



    @inbounds for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        fill!(X_sum, 0.0)
        for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki, kj, m, xk =
                S_ki[k_spl, Rij], S_kj[k_spl, Rij], S_m[k_spl, Rij], S_xk[k_spl, Rij]

            Ptm = @SMatrix [m * Props[i, j, xk] for i = 1:3, j = 1:3]

            X_sum[fd.yy] +=
                -V12[fd.yy, ki] * V34[fd.yy, kj] * Ptm[2, 2] -
                V12[fd.yz1, ki] * V34[fd.zy1, kj] * Ptm[3, 3] -
                V12[fd.yx1, ki] * V34[fd.xy1, kj] * Ptm[1, 1]
            X_sum[fd.zz] +=
                -V12[fd.zz, ki] * V34[fd.zz, kj] * Ptm[3, 3] -
                V12[fd.zx1, ki] * V34[fd.xz1, kj] * Ptm[1, 1] -
                V12[fd.zy1, ki] * V34[fd.yz1, kj] * Ptm[2, 2]
            X_sum[fd.xx] +=
                -V12[fd.xx, ki] * V34[fd.xx, kj] * Ptm[1, 1] -
                V12[fd.xy1, ki] * V34[fd.yx1, kj] * Ptm[2, 2] -
                V12[fd.xz1, ki] * V34[fd.zx1, kj] * Ptm[3, 3]

            ### Xab1 = -Vaa Vab1 - Vab1 Vbb - Vac1 Vcb1
            X_sum[fd.xy1] +=
                -V12[fd.xx, ki] * V34[fd.xy1, kj] * Ptm[1, 1] -
                V12[fd.xy1, ki] * V34[fd.yy, kj] * Ptm[2, 2] -
                V12[fd.xz1, ki] * V34[fd.zy1, kj] * Ptm[3, 3]
            X_sum[fd.xz1] +=
                -V12[fd.xx, ki] * V34[fd.xz1, kj] * Ptm[1, 1] -
                V12[fd.xz1, ki] * V34[fd.zz, kj] * Ptm[3, 3] -
                V12[fd.xy1, ki] * V34[fd.yz1, kj] * Ptm[2, 2]
            X_sum[fd.yx1] +=
                -V12[fd.yy, ki] * V34[fd.yx1, kj] * Ptm[2, 2] -
                V12[fd.yx1, ki] * V34[fd.xx, kj] * Ptm[1, 1] -
                V12[fd.yz1, ki] * V34[fd.zx1, kj] * Ptm[3, 3]
            X_sum[fd.yz1] +=
                -V12[fd.yy, ki] * V34[fd.yz1, kj] * Ptm[2, 2] -
                V12[fd.yz1, ki] * V34[fd.zz, kj] * Ptm[3, 3] -
                V12[fd.yx1, ki] * V34[fd.xz1, kj] * Ptm[1, 1]
            X_sum[fd.zx1] +=
                -V12[fd.zz, ki] * V34[fd.zx1, kj] * Ptm[3, 3] -
                V12[fd.zx1, ki] * V34[fd.xx, kj] * Ptm[1, 1] -
                V12[fd.zy1, ki] * V34[fd.yx1, kj] * Ptm[2, 2]
            X_sum[fd.zy1] +=
                -V12[fd.zz, ki] * V34[fd.zy1, kj] * Ptm[3, 3] -
                V12[fd.zy1, ki] * V34[fd.yy, kj] * Ptm[2, 2] -
                V12[fd.zx1, ki] * V34[fd.xy1, kj] * Ptm[1, 1]

            ### Xab2 = -Vab2 Vab2 - Vab3 Vba3
            X_sum[fd.xy2] +=
                -V12[fd.xy2, ki] * V34[fd.xy2, kj] * Ptm[1, 2] -
                V12[fd.xy3, ki] * V34[fd.yx3, kj] * Ptm[2, 1]
            X_sum[fd.xz2] +=
                -V12[fd.xz2, ki] * V34[fd.xz2, kj] * Ptm[1, 3] -
                V12[fd.xz3, ki] * V34[fd.zx3, kj] * Ptm[3, 1]
            X_sum[fd.yx2] +=
                -V12[fd.yx2, ki] * V34[fd.yx2, kj] * Ptm[2, 1] -
                V12[fd.yx3, ki] * V34[fd.xy3, kj] * Ptm[1, 2]
            X_sum[fd.yz2] +=
                -V12[fd.yz2, ki] * V34[fd.yz2, kj] * Ptm[2, 3] -
                V12[fd.yz3, ki] * V34[fd.zy3, kj] * Ptm[3, 2]
            X_sum[fd.zx2] +=
                -V12[fd.zx2, ki] * V34[fd.zx2, kj] * Ptm[3, 1] -
                V12[fd.zx3, ki] * V34[fd.xz3, kj] * Ptm[1, 3]
            X_sum[fd.zy2] +=
                -V12[fd.zy2, ki] * V34[fd.zy2, kj] * Ptm[3, 2] -
                V12[fd.zy3, ki] * V34[fd.yz3, kj] * Ptm[2, 3]

            ### Xab3 = -Vab2 Vab3 - Vab3 Vba2
            X_sum[fd.xy3] +=
                -V12[fd.xy2, ki] * V34[fd.xy3, kj] * Ptm[1, 2] -
                V12[fd.xy3, ki] * V34[fd.yx2, kj] * Ptm[2, 1]
            X_sum[fd.xz3] +=
                -V12[fd.xz2, ki] * V34[fd.xz3, kj] * Ptm[1, 3] -
                V12[fd.xz3, ki] * V34[fd.zx2, kj] * Ptm[3, 1]
            X_sum[fd.yx3] +=
                -V12[fd.yx2, ki] * V34[fd.yx3, kj] * Ptm[2, 1] -
                V12[fd.yx3, ki] * V34[fd.xy2, kj] * Ptm[1, 2]
            X_sum[fd.yz3] +=
                -V12[fd.yz2, ki] * V34[fd.yz3, kj] * Ptm[2, 3] -
                V12[fd.yz3, ki] * V34[fd.zy2, kj] * Ptm[3, 2]
            X_sum[fd.zx3] +=
                -V12[fd.zx2, ki] * V34[fd.zx3, kj] * Ptm[3, 1] -
                V12[fd.zx3, ki] * V34[fd.xz2, kj] * Ptm[1, 3]
            X_sum[fd.zy3] +=
                -V12[fd.zy2, ki] * V34[fd.zy3, kj] * Ptm[3, 2] -
                V12[fd.zy3, ki] * V34[fd.yz2, kj] * Ptm[2, 3]
        end

        (@view X[1:21, Rij, is, it, iu]) .+= X_sum
    end
    return
end

function addY!(
    X::Array{T,5},
    Gamma::Array{T,5},
    System::Geometry,
    N::Int64,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Props::Array{T,4},
    Buffers::ThreadLocalBuffersT,
) where {T}
    (; Npairs, invpairs, PairTypes, OnsitePairs) = System

    Vert(n, Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, n, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    _, wpw2, _, wpw4, wmw1, _, wmw3, _ = mixedFrequencies(ns, nt, nu, nwpr)
    flavTransf13 = (nt * wmw3 < 0, wmw1 * wmw3 > 0, wmw1 * nt > 0)
    flavTransf24 = (nt * wpw4 < 0, wpw2 * wpw4 > 0, wpw2 * nt > 0)
    flavTransf31 = (nt * wmw1 > 0, wmw3 * wmw1 > 0, wmw3 * nt < 0)
    flavTransf42 = (nt * wpw2 > 0, wpw4 * wpw2 > 0, wpw4 * nt < 0)

    (; X_sum, V13_addY, V24_addY, V31_addY, V42_addY) = Buffers

    V13 = V13_addY
    V24 = V24_addY
    V31 = V31_addY
    V42 = V42_addY


    # Xtilde only defined for nonlocal pairs Rij != Rii
    @inbounds for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        # loop over all left hand side inequivalent pairs Rij
        #Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij) 
        (; xi, xj) = PairTypes[Rij]


        for n = 1:21
            V13[n] = Vert(n, Rij, -wmw1, nt, wmw3, flavTransf13)
            V24[n] = Vert(n, Rij, wpw2, -nt, -wpw4, flavTransf24)
            V31[n] = Vert(n, Rij, wmw3, nt, -wmw1, flavTransf31)
            V42[n] = Vert(n, Rij, -wpw4, -nt, wpw2, flavTransf42)
        end

        fill!(X_sum, 0.0)

        P = @SMatrix [Props[i, j, xi, xj] for i = 1:3, j = 1:3]
        PT = @SMatrix [Props[j, i, xj, xi] for i = 1:3, j = 1:3]

        ### Yaa = Vaa Vaa + Vab2 Vab2 + Vac2 Vac2 + (w -- -w + t)

        X_sum[fd.xx] += (
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

        X_sum[fd.yy] += (
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

        X_sum[fd.zz] += (
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

        X_sum[fd.xy1] += (
            (V13[fd.xy3] * V24[fd.xy3] * P[2, 1] + V13[fd.xy1] * V24[fd.xy1] * P[1, 2]) + (
                V31[fd.xy3] * V42[fd.xy3] * PT[2, 1] + V31[fd.xy1] * V42[fd.xy1] * PT[1, 2]
            )
        )

        X_sum[fd.xz1] += (
            (V13[fd.xz3] * V24[fd.xz3] * P[3, 1] + V13[fd.xz1] * V24[fd.xz1] * P[1, 3]) + (
                V31[fd.xz3] * V42[fd.xz3] * PT[3, 1] + V31[fd.xz1] * V42[fd.xz1] * PT[1, 3]
            )
        )

        X_sum[fd.yx1] += (
            (V13[fd.yx3] * V24[fd.yx3] * P[1, 2] + V13[fd.yx1] * V24[fd.yx1] * P[2, 1]) + (
                V31[fd.yx3] * V42[fd.yx3] * PT[1, 2] + V31[fd.yx1] * V42[fd.yx1] * PT[2, 1]
            )
        )

        X_sum[fd.yz1] += (
            (V13[fd.yz3] * V24[fd.yz3] * P[3, 2] + V13[fd.yz1] * V24[fd.yz1] * P[2, 3]) + (
                V31[fd.yz3] * V42[fd.yz3] * PT[3, 2] + V31[fd.yz1] * V42[fd.yz1] * PT[2, 3]
            )
        )

        X_sum[fd.zx1] += (
            (V13[fd.zx3] * V24[fd.zx3] * P[1, 3] + V13[fd.zx1] * V24[fd.zx1] * P[3, 1]) + (
                V31[fd.zx3] * V42[fd.zx3] * PT[1, 3] + V31[fd.zx1] * V42[fd.zx1] * PT[3, 1]
            )
        )

        X_sum[fd.zy1] += (
            (V13[fd.zy3] * V24[fd.zy3] * P[2, 3] + V13[fd.zy1] * V24[fd.zy1] * P[3, 2]) + (
                V31[fd.zy3] * V42[fd.zy3] * PT[2, 3] + V31[fd.zy1] * V42[fd.zy1] * PT[3, 2]
            )
        )

        ### Yab2 = Vaa Vba2 + Vab2 Vbb + Vac2 Vbc2 + (w -- -w + t)

        X_sum[fd.xy2] += (
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

        X_sum[fd.xz2] += (
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

        X_sum[fd.yx2] += (
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

        X_sum[fd.yz2] += (
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

        X_sum[fd.zx2] += (
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

        X_sum[fd.zy2] += (
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

        X_sum[fd.xy3] += (
            (V13[fd.xy3] * V24[fd.yx1] * P[2, 1] + V13[fd.xy1] * V24[fd.yx3] * P[1, 2]) + (
                V31[fd.xy3] * V42[fd.yx1] * PT[2, 1] + V31[fd.xy1] * V42[fd.yx3] * PT[1, 2]
            )
        )

        X_sum[fd.xz3] += (
            (V13[fd.xz3] * V24[fd.zx1] * P[3, 1] + V13[fd.xz1] * V24[fd.zx3] * P[1, 3]) + (
                V31[fd.xz3] * V42[fd.zx1] * PT[3, 1] + V31[fd.xz1] * V42[fd.zx3] * PT[1, 3]
            )
        )

        X_sum[fd.yx3] += (
            (V13[fd.yx3] * V24[fd.xy1] * P[1, 2] + V13[fd.yx1] * V24[fd.xy3] * P[2, 1]) + (
                V31[fd.yx3] * V42[fd.xy1] * PT[1, 2] + V31[fd.yx1] * V42[fd.xy3] * PT[2, 1]
            )
        )

        X_sum[fd.yz3] += (
            (V13[fd.yz3] * V24[fd.zy1] * P[3, 2] + V13[fd.yz1] * V24[fd.zy3] * P[2, 3]) + (
                V31[fd.yz3] * V42[fd.zy1] * PT[3, 2] + V31[fd.yz1] * V42[fd.zy3] * PT[2, 3]
            )
        )

        X_sum[fd.zx3] += (
            (V13[fd.zx3] * V24[fd.xz1] * P[1, 3] + V13[fd.zx1] * V24[fd.xz3] * P[3, 1]) + (
                V31[fd.zx3] * V42[fd.xz1] * PT[1, 3] + V31[fd.zx1] * V42[fd.xz3] * PT[3, 1]
            )
        )

        X_sum[fd.zy3] += (
            (V13[fd.zy3] * V24[fd.yz1] * P[2, 3] + V13[fd.zy1] * V24[fd.yz3] * P[3, 2]) + (
                V31[fd.zy3] * V42[fd.yz1] * PT[2, 3] + V31[fd.zy1] * V42[fd.yz3] * PT[3, 2]
            )
        )

        (@view X[22:end, Rij, is, it, iu]) .+= X_sum
    end
end


function getXBubble!(Workspace::OneLoopWorkspace, T::Real)
    Par = Workspace.Par
    (; N, lenIntw) = Par.NumericalParams
    (; NUnique) = Par.System

    iSigma = Workspace.State.iSigma
    DiSigma = Workspace.Deriv.iSigma

    iG = SVector{3}([
        (x, nw) -> iG_(iSigma_i, x, nw, T) for iSigma_i in (iSigma.x, iSigma.y, iSigma.z)
    ])
    iSKat = SVector{3}([
        (x, nw) -> iSKat_(iSigma_i, DiSigma_i, x, nw, T) for (iSigma_i, DiSigma_i) in
        zip((iSigma.x, iSigma.y, iSigma.z), (DiSigma.x, DiSigma.y, DiSigma.z))
    ])

    ThreadLocalBuffers = get_ThreadLocalBuffers(Par.System)

    Threads.@threads :static for is_it = 1:N*N
        @inbounds begin
            is = div(is_it - 1, N) + 1
            it = (is_it - 1) % N + 1
            # WARNING: 
            # This works only with :static
            Buffs = ThreadLocalBuffers[Threads.threadid()]
            ns = is - 1
            nt = it - 1

            # Precompute Katanin propagators
            for Rij = 1:NUnique
                for j = 1:3, i = 1:3
                    Buffs.spropX[i, j, Rij] = -iSKat[i](Rij, 0) * iG[j](Rij, 0)  # nw will be updated later
                end
            end

            for nw = -lenIntw:lenIntw-1 # Matsubara sum
                nw_ns = nw + ns
                nw_nt = nw - nt
                # Update Katanin propagators for current nw
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

                for iu = 1:N
                    nu = iu - 1
                    if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                        continue
                    end
                    addY!(
                        Workspace.X,
                        Workspace.State.Gamma,
                        Workspace.Par.System,
                        N,
                        is,
                        it,
                        iu,
                        nw,
                        Buffs.spropY,
                        Buffs,
                    ) # add to XTilde-type bubble functions

                    ### If no u--t symmetry, then add all the bubbles
                    ### If use u--t symmetry, then only add for nu smaller then nt (all other obtained by symmetry)
                    # if(!Par.Options.use_symmetry || nu<=nt)
                    addX!(
                        Workspace.X,
                        Workspace.State.Gamma,
                        Workspace.Par.System,
                        N,
                        is,
                        it,
                        iu,
                        nw,
                        Buffs.spropX,
                        Buffs,
                    )
                    # end
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
                    jsum[fd.xx] +=
                        (
                            gam[fd.xx] * Props[1](xk, nw) +
                            gam[fd.yx1] * Props[2](xk, nw) +
                            gam[fd.zx1] * Props[3](xk, nw)
                        ) * m
                    jsum[fd.yy] +=
                        (
                            gam[fd.xy1] * Props[1](xk, nw) +
                            gam[fd.yy] * Props[2](xk, nw) +
                            gam[fd.zy1] * Props[3](xk, nw)
                        ) * m
                    jsum[fd.zz] +=
                        (
                            gam[fd.xz1] * Props[1](xk, nw) +
                            gam[fd.yz1] * Props[2](xk, nw) +
                            gam[fd.zz] * Props[3](xk, nw)
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
        Gamma[fd.yz2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 1]
        Gamma[fd.zy2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 1]
        Gamma[fd.zx2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 2]
        Gamma[fd.xz2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 2]
        Gamma[fd.xy2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 3]
        Gamma[fd.yx2, Rj, :, :, :] .= -couplings[Rj] * isotropy[Rj, 3]

        Gamma[fd.yz3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 1]
        Gamma[fd.zy3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 1]
        Gamma[fd.zx3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 2]
        Gamma[fd.xz3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 2]
        Gamma[fd.xy3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 3]
        Gamma[fd.yx3, Rj, :, :, :] .= couplings[Rj] * isotropy[Rj, 3]
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

export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
