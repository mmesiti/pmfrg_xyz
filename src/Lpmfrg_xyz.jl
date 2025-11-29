module pmfrg_xyz

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

using RecursiveArrayTools
using SpinFRGLattices, OrdinaryDiffEq, DiffEqCallbacks, RecursiveArrayTools, StructArrays
using SpinFRGLattices.StaticArrays

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

struct SigmaType{T}
    x::Array{T,2}
    y::Array{T,2}
    z::Array{T,2}
end

struct StateType{T}
    f_int::Vector{T}
    iSigma::SigmaType{T}
    Gamma::Array{T,5}
end

struct Observables{T}
    Chi_x::Vector{T}
    Chi_y::Vector{T}
    Chi_z::Vector{T}
end

struct NumericalParams{T<:Real}
    T::T
    N::Int

    accuracy::T
    lambda_min::T
    lambda_max::T

    lenIntw::Int
    lenIntw_acc::Int
end

struct OptionParams
    use_symmetry::Bool
    minimal_output::Bool
end

struct OneLoopParams_1{T,SType}
    System::SType
    NumericalParams::NumericalParams{T}
    Options::OptionParams
end

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
_getFloatType(Par) = typeof(Par.NumericalParams.T)

function SigmaType(NUnique::Int, N::Int, type = Float64)
    return SigmaType(
        zeros(type, NUnique, N),
        zeros(type, NUnique, N),
        zeros(type, NUnique, N),
    )
end
SigmaType(Par) = SigmaType(Par.System.NPairs, Par.NumericalParams.N)

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

function NumericalParams(;
    T::Real = 0.5,
    N::Integer = 24,
    accuracy = 1e-6,
    lambda_min = exp(-10.0),
    lambda_max = exp(10.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = 2 * maximum((N, lenIntw)),
)

    return NumericalParams(T, N, accuracy, lambda_min, lambda_max, lenIntw, lenIntw_acc)
end

OptionParams(; use_symmetry::Bool = true, MinimalOutput::Bool = false, kwargs...) =
    OptionParams(use_symmetry, MinimalOutput)
Params(System; kwargs...) =
    OneLoopParams_1(System, NumericalParams(; kwargs...), OptionParams(; kwargs...))

#############################################################
######### PROPAGATORS ## PROPAGATORS ## PROPAGATORS #########
#############################################################

### Propagators will depend on an additional flavor
### Instead of modifying the propagators, I will simply use them
### as V_, by doing iG_(iSigma.x, ...)

function get_w(nw, T)
    return pi * T * (2 * nw + 1)
end

function get_sign_iw(nw::Integer, N::Integer)
    # s = sign(nw)
    nw_bounds = min(nw, N - 1)  ### used to be min(abs(nw),...), but nw is set positive in gamma
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

function iG_(iSigma::AbstractArray, x::Integer, Lam::Real, nw::Integer, T::Real)
    w = get_w(nw, T)
    return w / (w^2 + w * iSigma_(iSigma, x, nw) + Lam^2)
end

function iS_(iSigma::AbstractArray, x::Integer, Lam::Real, nw::Integer, T::Real)
    w = get_w(nw, T)
    return -iG_(iSigma, x, Lam, nw, T)^2 * 2 * Lam / w
end

function iSKat_(
    iSigma::AbstractArray,
    DSigma::AbstractArray,
    x::Integer,
    Lam::Real,
    nw::Integer,
    T::Real,
)
    w = get_w(nw, T)
    return -iG_(iSigma, x, Lam, nw, T)^2 * (2 * Lam / w + iSigma_(DSigma, x, nw))
end

####################################################
######### VERTICES ## VERTICES ## VERTICES #########
####################################################

### Symmetries: 
###     s <--> -s
###     t <--> -t, i <--> j
###     u <--> -u, i <--> j
function ConvertFreqArgs(ns, nt, nu, Nw)
    swapsites = nt * nu < 0
    ns, nt, nu = abs.((ns, nt, nu))

    ns = min(ns, Nw - 1 - (ns + Nw - 1) % 2) ### weird cutoff, idk why
    nt = min(nt, Nw - 1 - (nt + Nw - 1) % 2)
    nu = min(nu, Nw - 1 - (nu + Nw - 1) % 2)

    return ns, nt, nu, swapsites
end

using LinearAlgebra
using SparseArrays

function V_(
    Vertex::AbstractArray,
    n::Int,
    ns::Int,
    nt::Int,
    nu::Int,
    Rij::Integer,
    Rji::Integer,
    N::Integer;
    isX = false,
)
    isFlavorTransform = (nt * nu < 0, ns * nu < 0, ns * nt < 0)

    block = div(n + 2, 6)

    n_transf = n
    if (block != 0)
        if (isFlavorTransform[block])
            n_transf = ((n - 3) - (block - 1) * 6 + 2) % 6 + 1 + 3 + (block - 1) * 6
        end
    end

    ns, nt, nu, swapsites = ConvertFreqArgs(ns, nt, nu, N)
    Rij = ifelse(swapsites, Rji, Rij)
    return Vertex[n_transf, Rij, ns+1, nt+1, nu+1]
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

function addX!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
    (; State, X, Par) = Workspace
    N = Par.NumericalParams.N
    (; Npairs, Nsum, siteSum, invpairs) = Par.System

    Vert(n, Rij, s, t, u) = V_(State.Gamma, n, s, t, u, Rij, invpairs[Rij], N; isX = true)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    # get fields of siteSum struct as Matrices for better use of LoopVectorization
    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    X_sum = zeros(42)
    for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        fill!(X_sum, 0.0)
        sumsum = 0
        for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki, kj, m, xk =
                S_ki[k_spl, Rij], S_kj[k_spl, Rij], S_m[k_spl, Rij], S_xk[k_spl, Rij]
            Ptm = @SMatrix [
                Props[xk, xk, 1, 1] Props[xk, xk, 1, 2] Props[xk, xk, 1, 3]
                Props[xk, xk, 2, 1] Props[xk, xk, 2, 2] Props[xk, xk, 2, 3]
                Props[xk, xk, 3, 1] Props[xk, xk, 3, 2] Props[xk, xk, 3, 3]
            ]
            Ptm = Ptm * m ### Props now contains two flavor indices

            V12(n) = Vert(n, ki, ns, wpw1, -wpw2)
            V34(n) = Vert(n, kj, ns, -wmw3, -wmw4)

            X_sum[fd.yy] +=
                -V12(fd.yy) * V34(fd.yy) * Ptm[2, 2] -
                V12(fd.yz1) * V34(fd.zy1) * Ptm[3, 3] -
                V12(fd.yx1) * V34(fd.xy1) * Ptm[1, 1]
            X_sum[fd.zz] +=
                -V12(fd.zz) * V34(fd.zz) * Ptm[3, 3] -
                V12(fd.zx1) * V34(fd.xz1) * Ptm[1, 1] -
                V12(fd.zy1) * V34(fd.yz1) * Ptm[2, 2]
            X_sum[fd.xx] +=
                -V12(fd.xx) * V34(fd.xx) * Ptm[1, 1] -
                V12(fd.xy1) * V34(fd.yx1) * Ptm[2, 2] -
                V12(fd.xz1) * V34(fd.zx1) * Ptm[3, 3]

            ### Xab1 = -Vaa Vab1 - Vab1 Vbb - Vac1 Vcb1
            X_sum[fd.xy1] +=
                -V12(fd.xx) * V34(fd.xy1) * Ptm[1, 1] -
                V12(fd.xy1) * V34(fd.yy) * Ptm[2, 2] - V12(fd.xz1) * V34(fd.zy1) * Ptm[3, 3]
            X_sum[fd.xz1] +=
                -V12(fd.xx) * V34(fd.xz1) * Ptm[1, 1] -
                V12(fd.xz1) * V34(fd.zz) * Ptm[3, 3] - V12(fd.xy1) * V34(fd.yz1) * Ptm[2, 2]
            X_sum[fd.yx1] +=
                -V12(fd.yy) * V34(fd.yx1) * Ptm[2, 2] -
                V12(fd.yx1) * V34(fd.xx) * Ptm[1, 1] - V12(fd.yz1) * V34(fd.zx1) * Ptm[3, 3]
            X_sum[fd.yz1] +=
                -V12(fd.yy) * V34(fd.yz1) * Ptm[2, 2] -
                V12(fd.yz1) * V34(fd.zz) * Ptm[3, 3] - V12(fd.yx1) * V34(fd.xz1) * Ptm[1, 1]
            X_sum[fd.zx1] +=
                -V12(fd.zz) * V34(fd.zx1) * Ptm[3, 3] -
                V12(fd.zx1) * V34(fd.xx) * Ptm[1, 1] - V12(fd.zy1) * V34(fd.yx1) * Ptm[2, 2]
            X_sum[fd.zy1] +=
                -V12(fd.zz) * V34(fd.zy1) * Ptm[3, 3] -
                V12(fd.zy1) * V34(fd.yy) * Ptm[2, 2] - V12(fd.zx1) * V34(fd.xy1) * Ptm[1, 1]

            ### Xab2 = -Vab2 Vab2 - Vab3 Vba3
            X_sum[fd.xy2] +=
                -V12(fd.xy2) * V34(fd.xy2) * Ptm[1, 2] -
                V12(fd.xy3) * V34(fd.yx3) * Ptm[2, 1]
            X_sum[fd.xz2] +=
                -V12(fd.xz2) * V34(fd.xz2) * Ptm[1, 3] -
                V12(fd.xz3) * V34(fd.zx3) * Ptm[3, 1]
            X_sum[fd.yx2] +=
                -V12(fd.yx2) * V34(fd.yx2) * Ptm[2, 1] -
                V12(fd.yx3) * V34(fd.xy3) * Ptm[1, 2]
            X_sum[fd.yz2] +=
                -V12(fd.yz2) * V34(fd.yz2) * Ptm[2, 3] -
                V12(fd.yz3) * V34(fd.zy3) * Ptm[3, 2]
            X_sum[fd.zx2] +=
                -V12(fd.zx2) * V34(fd.zx2) * Ptm[3, 1] -
                V12(fd.zx3) * V34(fd.xz3) * Ptm[1, 3]
            X_sum[fd.zy2] +=
                -V12(fd.zy2) * V34(fd.zy2) * Ptm[3, 2] -
                V12(fd.zy3) * V34(fd.yz3) * Ptm[2, 3]

            ### Xab3 = -Vab2 Vab3 - Vab3 Vba2
            X_sum[fd.xy3] +=
                -V12(fd.xy2) * V34(fd.xy3) * Ptm[1, 2] -
                V12(fd.xy3) * V34(fd.yx2) * Ptm[2, 1]
            X_sum[fd.xz3] +=
                -V12(fd.xz2) * V34(fd.xz3) * Ptm[1, 3] -
                V12(fd.xz3) * V34(fd.zx2) * Ptm[3, 1]
            X_sum[fd.yx3] +=
                -V12(fd.yx2) * V34(fd.yx3) * Ptm[2, 1] -
                V12(fd.yx3) * V34(fd.xy2) * Ptm[1, 2]
            X_sum[fd.yz3] +=
                -V12(fd.yz2) * V34(fd.yz3) * Ptm[2, 3] -
                V12(fd.yz3) * V34(fd.zy2) * Ptm[3, 2]
            X_sum[fd.zx3] +=
                -V12(fd.zx2) * V34(fd.zx3) * Ptm[3, 1] -
                V12(fd.zx3) * V34(fd.xz2) * Ptm[1, 3]
            X_sum[fd.zy3] +=
                -V12(fd.zy2) * V34(fd.zy3) * Ptm[3, 2] -
                V12(fd.zy3) * V34(fd.yz2) * Ptm[2, 3]
        end

        X[:, Rij, is, it, iu] .+= X_sum
    end
    return
end

function addY!(
    Workspace,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Props;
    _l = 1.0,
)
    (; State, X, Par) = Workspace
    N = Par.NumericalParams.N
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    Vert(n, Rij, s, t, u) = V_(State.Gamma, n, s, t, u, Rij, invpairs[Rij], N)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)

    X_sum = zeros(42)

    # Xtilde only defined for nonlocal pairs Rij != Rii
    for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        # loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij) 
        (; xi, xj) = PairTypes[Rij]

        function P_(n::Int, m::Int)
            return Props[xi, xj, n, m]
        end

        function PT_(n::Int, m::Int)
            return Props[xj, xi, m, n]
        end

        V13(n) = Vert(n, Rij, -wmw1, nt, wmw3)
        V24(n) = Vert(n, Rij, wpw2, -nt, -wpw4)

        V31(n) = Vert(n, Rij, wmw3, nt, -wmw1)
        V42(n) = Vert(n, Rij, -wpw4, -nt, wpw2)

        fill!(X_sum, 0.0)

        ### Yaa = Vaa Vaa + Vab2 Vab2 + Vac2 Vac2 + (w -- -w + t)

        X_sum[21+fd.xx] += (
            (
                V13(fd.xx) * V24(fd.xx) * P_(1, 1) +
                V13(fd.xy2) * V24(fd.xy2) * P_(2, 2) +
                V13(fd.xz2) * V24(fd.xz2) * P_(3, 3)
            ) + (
                V31(fd.xx) * V42(fd.xx) * PT_(1, 1) +
                V31(fd.xy2) * V42(fd.xy2) * PT_(2, 2) +
                V31(fd.xz2) * V42(fd.xz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.yy] += (
            (
                V13(fd.yy) * V24(fd.yy) * P_(2, 2) +
                V13(fd.yx2) * V24(fd.yx2) * P_(1, 1) +
                V13(fd.yz2) * V24(fd.yz2) * P_(3, 3)
            ) + (
                V31(fd.yy) * V42(fd.yy) * PT_(2, 2) +
                V31(fd.yx2) * V42(fd.yx2) * PT_(1, 1) +
                V31(fd.yz2) * V42(fd.yz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.zz] += (
            (
                V13(fd.zz) * V24(fd.zz) * P_(3, 3) +
                V13(fd.zx2) * V24(fd.zx2) * P_(1, 1) +
                V13(fd.zy2) * V24(fd.zy2) * P_(2, 2)
            ) + (
                V31(fd.zz) * V42(fd.zz) * PT_(3, 3) +
                V31(fd.zx2) * V42(fd.zx2) * PT_(1, 1) +
                V31(fd.zy2) * V42(fd.zy2) * PT_(2, 2)
            )
        )

        ### Yab1 = Vab3 Vab3 + Vab1 Vab1 + (w -- -w + t)

        X_sum[21+fd.xy1] += (
            (V13(fd.xy3) * V24(fd.xy3) * P_(2, 1) + V13(fd.xy1) * V24(fd.xy1) * P_(1, 2)) + (
                V31(fd.xy3) * V42(fd.xy3) * PT_(2, 1) +
                V31(fd.xy1) * V42(fd.xy1) * PT_(1, 2)
            )
        )

        X_sum[21+fd.xz1] += (
            (V13(fd.xz3) * V24(fd.xz3) * P_(3, 1) + V13(fd.xz1) * V24(fd.xz1) * P_(1, 3)) + (
                V31(fd.xz3) * V42(fd.xz3) * PT_(3, 1) +
                V31(fd.xz1) * V42(fd.xz1) * PT_(1, 3)
            )
        )

        X_sum[21+fd.yx1] += (
            (V13(fd.yx3) * V24(fd.yx3) * P_(1, 2) + V13(fd.yx1) * V24(fd.yx1) * P_(2, 1)) + (
                V31(fd.yx3) * V42(fd.yx3) * PT_(1, 2) +
                V31(fd.yx1) * V42(fd.yx1) * PT_(2, 1)
            )
        )

        X_sum[21+fd.yz1] += (
            (V13(fd.yz3) * V24(fd.yz3) * P_(3, 2) + V13(fd.yz1) * V24(fd.yz1) * P_(2, 3)) + (
                V31(fd.yz3) * V42(fd.yz3) * PT_(3, 2) +
                V31(fd.yz1) * V42(fd.yz1) * PT_(2, 3)
            )
        )

        X_sum[21+fd.zx1] += (
            (V13(fd.zx3) * V24(fd.zx3) * P_(1, 3) + V13(fd.zx1) * V24(fd.zx1) * P_(3, 1)) + (
                V31(fd.zx3) * V42(fd.zx3) * PT_(1, 3) +
                V31(fd.zx1) * V42(fd.zx1) * PT_(3, 1)
            )
        )

        X_sum[21+fd.zy1] += (
            (V13(fd.zy3) * V24(fd.zy3) * P_(2, 3) + V13(fd.zy1) * V24(fd.zy1) * P_(3, 2)) + (
                V31(fd.zy3) * V42(fd.zy3) * PT_(2, 3) +
                V31(fd.zy1) * V42(fd.zy1) * PT_(3, 2)
            )
        )

        ### Yab2 = Vaa Vba2 + Vab2 Vbb + Vac2 Vbc2 + (w -- -w + t)

        X_sum[21+fd.xy2] += (
            (
                V13(fd.xx) * V24(fd.yx2) * P_(1, 1) +
                V13(fd.xy2) * V24(fd.yy) * P_(2, 2) +
                V13(fd.xz2) * V24(fd.yz2) * P_(3, 3)
            ) + (
                V31(fd.xx) * V42(fd.yx2) * PT_(1, 1) +
                V31(fd.xy2) * V42(fd.yy) * PT_(2, 2) +
                V31(fd.xz2) * V42(fd.yz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.xz2] += (
            (
                V13(fd.xx) * V24(fd.zx2) * P_(1, 1) +
                V13(fd.xz2) * V24(fd.zz) * P_(3, 3) +
                V13(fd.xy2) * V24(fd.zy2) * P_(2, 2)
            ) + (
                V31(fd.xx) * V42(fd.zx2) * PT_(1, 1) +
                V31(fd.xz2) * V42(fd.zz) * PT_(3, 3) +
                V31(fd.xy2) * V42(fd.zy2) * PT_(2, 2)
            )
        )

        X_sum[21+fd.yx2] += (
            (
                V13(fd.yy) * V24(fd.xy2) * P_(2, 2) +
                V13(fd.yx2) * V24(fd.xx) * P_(1, 1) +
                V13(fd.yz2) * V24(fd.xz2) * P_(3, 3)
            ) + (
                V31(fd.yy) * V42(fd.xy2) * PT_(2, 2) +
                V31(fd.yx2) * V42(fd.xx) * PT_(1, 1) +
                V31(fd.yz2) * V42(fd.xz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.yz2] += (
            (
                V13(fd.yy) * V24(fd.zy2) * P_(2, 2) +
                V13(fd.yz2) * V24(fd.zz) * P_(3, 3) +
                V13(fd.yx2) * V24(fd.zx2) * P_(1, 1)
            ) + (
                V31(fd.yy) * V42(fd.zy2) * PT_(2, 2) +
                V31(fd.yz2) * V42(fd.zz) * PT_(3, 3) +
                V31(fd.yx2) * V42(fd.zx2) * PT_(1, 1)
            )
        )

        X_sum[21+fd.zx2] += (
            (
                V13(fd.zz) * V24(fd.xz2) * P_(3, 3) +
                V13(fd.zx2) * V24(fd.xx) * P_(1, 1) +
                V13(fd.zy2) * V24(fd.xy2) * P_(2, 2)
            ) + (
                V31(fd.zz) * V42(fd.xz2) * PT_(3, 3) +
                V31(fd.zx2) * V42(fd.xx) * PT_(1, 1) +
                V31(fd.zy2) * V42(fd.xy2) * PT_(2, 2)
            )
        )

        X_sum[21+fd.zy2] += (
            (
                V13(fd.zz) * V24(fd.yz2) * P_(3, 3) +
                V13(fd.zy2) * V24(fd.yy) * P_(2, 2) +
                V13(fd.zx2) * V24(fd.yx2) * P_(1, 1)
            ) + (
                V31(fd.zz) * V42(fd.yz2) * PT_(3, 3) +
                V31(fd.zy2) * V42(fd.yy) * PT_(2, 2) +
                V31(fd.zx2) * V42(fd.yx2) * PT_(1, 1)
            )
        )

        ### Yab3 = Vab3 Vba1 + Vab1 Vba3 + (w -- -w + t)

        X_sum[21+fd.xy3] += (
            (V13(fd.xy3) * V24(fd.yx1) * P_(2, 1) + V13(fd.xy1) * V24(fd.yx3) * P_(1, 2)) + (
                V31(fd.xy3) * V42(fd.yx1) * PT_(2, 1) +
                V31(fd.xy1) * V42(fd.yx3) * PT_(1, 2)
            )
        )

        X_sum[21+fd.xz3] += (
            (V13(fd.xz3) * V24(fd.zx1) * P_(3, 1) + V13(fd.xz1) * V24(fd.zx3) * P_(1, 3)) + (
                V31(fd.xz3) * V42(fd.zx1) * PT_(3, 1) +
                V31(fd.xz1) * V42(fd.zx3) * PT_(1, 3)
            )
        )

        X_sum[21+fd.yx3] += (
            (V13(fd.yx3) * V24(fd.xy1) * P_(1, 2) + V13(fd.yx1) * V24(fd.xy3) * P_(2, 1)) + (
                V31(fd.yx3) * V42(fd.xy1) * PT_(1, 2) +
                V31(fd.yx1) * V42(fd.xy3) * PT_(2, 1)
            )
        )

        X_sum[21+fd.yz3] += (
            (V13(fd.yz3) * V24(fd.zy1) * P_(3, 2) + V13(fd.yz1) * V24(fd.zy3) * P_(2, 3)) + (
                V31(fd.yz3) * V42(fd.zy1) * PT_(3, 2) +
                V31(fd.yz1) * V42(fd.zy3) * PT_(2, 3)
            )
        )

        X_sum[21+fd.zx3] += (
            (V13(fd.zx3) * V24(fd.xz1) * P_(1, 3) + V13(fd.zx1) * V24(fd.xz3) * P_(3, 1)) + (
                V31(fd.zx3) * V42(fd.xz1) * PT_(1, 3) +
                V31(fd.zx1) * V42(fd.xz3) * PT_(3, 1)
            )
        )

        X_sum[21+fd.zy3] += (
            (V13(fd.zy3) * V24(fd.yz1) * P_(2, 3) + V13(fd.zy1) * V24(fd.yz3) * P_(3, 2)) + (
                V31(fd.zy3) * V42(fd.yz1) * PT_(2, 3) +
                V31(fd.zy1) * V42(fd.yz3) * PT_(3, 2)
            )
        )

        X[:, Rij, is, it, iu] .+= X_sum
    end
end

function getXBubble!(Workspace, Lam)
    Par = Workspace.Par
    (; T, N, lenIntw) = Par.NumericalParams
    (; NUnique) = Par.System

    iGx(x, nw) = iG_(Workspace.State.iSigma.x, x, Lam, nw, T)
    iGy(x, nw) = iG_(Workspace.State.iSigma.y, x, Lam, nw, T)
    iGz(x, nw) = iG_(Workspace.State.iSigma.z, x, Lam, nw, T)

    iSKatx(x, nw) =
        iSKat_(Workspace.State.iSigma.x, Workspace.Deriv.iSigma.x, x, Lam, nw, T)
    iSKaty(x, nw) =
        iSKat_(Workspace.State.iSigma.y, Workspace.Deriv.iSigma.y, x, Lam, nw, T)
    iSKatz(x, nw) =
        iSKat_(Workspace.State.iSigma.z, Workspace.Deriv.iSigma.z, x, Lam, nw, T)

    function getKataninProp!(BubbleProp, nw1, nw2)
        for i = 1:Par.System.NUnique, j = 1:Par.System.NUnique
            BubbleProp[i, j, 1, 1] = iSKatx(i, nw1) * iGx(j, nw2) * T
            BubbleProp[i, j, 1, 2] = iSKatx(i, nw1) * iGy(j, nw2) * T
            BubbleProp[i, j, 1, 3] = iSKatx(i, nw1) * iGz(j, nw2) * T
            BubbleProp[i, j, 2, 1] = iSKaty(i, nw1) * iGx(j, nw2) * T
            BubbleProp[i, j, 2, 2] = iSKaty(i, nw1) * iGy(j, nw2) * T
            BubbleProp[i, j, 2, 3] = iSKaty(i, nw1) * iGz(j, nw2) * T
            BubbleProp[i, j, 3, 1] = iSKatz(i, nw1) * iGx(j, nw2) * T
            BubbleProp[i, j, 3, 2] = iSKatz(i, nw1) * iGy(j, nw2) * T
            BubbleProp[i, j, 3, 3] = iSKatz(i, nw1) * iGz(j, nw2) * T
        end

        ### Relative minus sign between paper & Nils' thesis
        return -BubbleProp
        # return SMatrix{NUnique, NUnique, 3, 3}(BubbleProp)
        ### SMatrix can only create 2d array (according to ChatGPT). Use SArray instead
    end

    for is = 1:N, it = 1:N
        BubbleProp = zeros(NUnique, NUnique, 3, 3)
        ns = is - 1
        nt = it - 1
        for nw = -lenIntw:lenIntw-1 # Matsubara sum
            spropX = getKataninProp!(BubbleProp, nw, nw + ns)
            spropY = getKataninProp!(BubbleProp, nw, nw - nt)
            for iu = 1:N
                nu = iu - 1
                if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                    continue
                end
                addY!(Workspace, is, it, iu, nw, spropY, _l = Lam) # add to XTilde-type bubble functions

                ### If no u--t symmetry, then add all the bubbles
                ### If use u--t symmetry, then only add for nu smaller then nt (all other obtained by symmetry)
                # if(!Par.Options.use_symmetry || nu<=nt)

                addX!(Workspace, is, it, iu, nw, spropX)
                # end
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

function getDFint!(Workspace, Lam::Real)
    (; State, Deriv, Par) = Workspace
    (; T, lenIntw_acc) = Par.NumericalParams
    NUnique = Par.System.NUnique

    iSigmax(x, nw) = iSigma_(State.iSigma.x, x, nw)
    iSigmay(x, nw) = iSigma_(State.iSigma.y, x, nw)
    iSigmaz(x, nw) = iSigma_(State.iSigma.z, x, nw)

    iGx(x, nw) = iG_(State.iSigma.x, x, Lam, nw, T)
    iGy(x, nw) = iG_(State.iSigma.y, x, Lam, nw, T)
    iGz(x, nw) = iG_(State.iSigma.z, x, Lam, nw, T)

    iSx(x, nw) = iS_(State.iSigma.x, x, Lam, nw, T)
    iSy(x, nw) = iS_(State.iSigma.y, x, Lam, nw, T)
    iSz(x, nw) = iS_(State.iSigma.z, x, Lam, nw, T)

    Theta(Lam, w) = w^2 / (w^2 + Lam^2)

    for x = 1:NUnique
        sumres = 0.0
        for nw = -lenIntw_acc:lenIntw_acc-1
            w = get_w(nw, T)
            sumres += iSx(x, nw) / iGy(x, nw) * Theta(Lam, w) * iSigmax(x, nw) / w
            sumres += iSy(x, nw) / iGy(x, nw) * Theta(Lam, w) * iSigmay(x, nw) / w
            sumres += iSz(x, nw) / iGz(x, nw) * Theta(Lam, w) * iSigmaz(x, nw) / w
        end
        Deriv.f_int[x] = -0.5 * T * sumres
    end
end

function get_Self_Energy!(Workspace, Lam)
    Par = Workspace.Par
    @inline iSx(x, nw) =
        iS_(Workspace.State.iSigma.x, x, Lam, nw, Par.NumericalParams.T) / 2
    @inline iSy(x, nw) =
        iS_(Workspace.State.iSigma.y, x, Lam, nw, Par.NumericalParams.T) / 2
    @inline iSz(x, nw) =
        iS_(Workspace.State.iSigma.z, x, Lam, nw, Par.NumericalParams.T) / 2
    compute1PartBubble!(Workspace.Deriv.iSigma, Workspace.State.Gamma, [iSx, iSy, iSz], Par)
end

function compute1PartBubble!(Dgamma::SigmaType, Gamma::Array{T,5}, Props, Par) where {T}
    invpairs = Par.System.invpairs

    setZero!(Dgamma)
    @inline Gamma_(n, Rij, s, t, u) =
        V_(Gamma, n, s, t, u, Rij, invpairs[Rij], Par.NumericalParams.N)
    addTo1PartBubble!(Dgamma, Gamma_, Props, Par)
end

function addTo1PartBubble!(Dgamma::SigmaType, Gamma_::Function, Props, Par)

    (; T, N, lenIntw_acc) = Par.NumericalParams
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
                    gam(n) = Gamma_(n, ki, 0, -wmw1, -wpw1)
                    jsum[fd.xx] +=
                        (
                            gam(fd.xx) * Props[1](xk, nw) +
                            gam(fd.yx1) * Props[2](xk, nw) +
                            gam(fd.zx1) * Props[3](xk, nw)
                        ) * m
                    jsum[fd.yy] +=
                        (
                            gam(fd.xy1) * Props[1](xk, nw) +
                            gam(fd.yy) * Props[2](xk, nw) +
                            gam(fd.zy1) * Props[3](xk, nw)
                        ) * m
                    jsum[fd.zz] +=
                        (
                            gam(fd.xz1) * Props[1](xk, nw) +
                            gam(fd.yz1) * Props[2](xk, nw) +
                            gam(fd.zz) * Props[3](xk, nw)
                        ) * m
                end
                Dgamma.x[x, iw1] += -T * jsum[1]
                Dgamma.y[x, iw1] += -T * jsum[2]
                Dgamma.z[x, iw1] += -T * jsum[3]
            end
        end
    end
end

struct OneLoopWorkspace{T,ParType}
    State::StateType{T}
    Deriv::StateType{T}
    X::Array{T,5}
    Par::ParType
end

function OneLoopWorkspace(State, Deriv, X, Par)
    setZero!(Deriv)
    setZero!(X)

    return OneLoopWorkspace(StateType(State.x...), StateType(Deriv.x...), X, Par)
end

using JLD2
function getDeriv!(Deriv, State, setup, Lam; saveArgs = true)

    (X, Par) = setup # use pre-allocated X and XTilde to reduce garbage collector time

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

function AllocateSetup(Par::OneLoopParams_1)
    println("One Loop: T= ", Par.NumericalParams.T)
    ## Allocate Memory:
    floattype = _getFloatType(Par)
    X = zeros(floattype, getBubbleVDims(Par))
    return (X, Par)
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

function launchPMFRG!(State, setup, Deriv!::Function; method = DP5(), npoints = 600)

    Par = setup[end]
    (; lambda_max, lambda_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(lambda_max)
    tend = get_t_min(lambda_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    saved_values = SavedValues(eltype(State), Observables{eltype(State)})

    function save_func(State, t, integrator)
        chi_x = getChi_x(State, t_to_Lam(t), Par)
        chi_y = getChi_y(State, t_to_Lam(t), Par)
        chi_z = getChi_z(State, t_to_Lam(t), Par)

        return Observables(copy(chi_x), copy(chi_y), copy(chi_z))
    end

    ObsSaveat = gettMesh(lambda_min, lambda_max, npoints)
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
        save_everystep = false,
        callback = saveCB,
        dt = Lam_to_t(0.2 * lambda_max),
    )

    return sol, saved_values
end

function testPMFRG!(State, setup, Deriv!::Function; loadArgs = false)
    Par = setup[end]
    (; lambda_max, lambda_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(lambda_max)
    tend = get_t_min(lambda_min)
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
    Lam < exp(-30) && @warn "lambda_min too small! Set to exp(-30) instead."
    max(Lam_to_t(Lam), -30.0)
end

function generateSubstituteDeriv(getDeriv!::Function)

    function DerivSubs!(Deriv, State, par, t; s = true)
        Lam = t_to_Lam(t)
        a = getDeriv!(Deriv, State, par, Lam, saveArgs = s)
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

getChi_z(State::ArrayPartition, Lam::Real, Par) =
    getChi_z(State.x[2], State.x[3], State.x[5], Lam, Par)
getChi_x(State::ArrayPartition, Lam::Real, Par) =
    getChi_x(State.x[3], State.x[4], State.x[5], Lam, Par)
getChi_y(State::ArrayPartition, Lam::Real, Par) =
    getChi_y(State.x[4], State.x[2], State.x[5], Lam, Par)

function getChi_z(
    iSigmaX::AbstractArray,
    iSigmaY::AbstractArray,
    Gamma::AbstractArray,
    Lam::Real,
    Par,
)
    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGx(x, w) = iG_(iSigmaX, x, Lam, w, T)
    iGy(x, w) = iG_(iSigmaY, x, Lam, w, T)
    Vxy2(Rij, s, t, u) = V_(Gamma, fd.xy2, s, t, u, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += T * iGx(xi, nK) * iGy(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGx(xi, nK)^2 * iGy(xj, nK2)^2
                Chi[Rij] += T^2 * GGGG * Vxy2(Rij, 0, npwpw2, -w2mw)
            end
        end
    end
    return (Chi)
end

function getChi_x(
    iSigmaY::AbstractArray,
    iSigmaZ::AbstractArray,
    Gamma::AbstractArray,
    Lam::Real,
    Par,
)
    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGy(x, w) = iG_(iSigmaY, x, Lam, w, T)
    iGz(x, w) = iG_(iSigmaZ, x, Lam, w, T)
    Vyz2(Rij, s, t, u) = V_(Gamma, fd.yz2, s, t, u, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += T * iGy(xi, nK) * iGz(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGy(xi, nK)^2 * iGz(xj, nK2)^2
                Chi[Rij] += T^2 * GGGG * Vyz2(Rij, 0, npwpw2, -w2mw)
            end
        end
    end
    return (Chi)
end

function getChi_y(
    iSigmaZ::AbstractArray,
    iSigmaX::AbstractArray,
    Gamma::AbstractArray,
    Lam::Real,
    Par,
)
    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGz(x, w) = iG_(iSigmaZ, x, Lam, w, T)
    iGx(x, w) = iG_(iSigmaX, x, Lam, w, T)
    Vzx2(Rij, s, t, u) = V_(Gamma, fd.zx2, s, t, u, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += T * iGz(xi, nK) * iGx(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGz(xi, nK)^2 * iGx(xj, nK2)^2
                Chi[Rij] += T^2 * GGGG * Vzx2(Rij, 0, npwpw2, -w2mw)
            end
        end
    end
    return (Chi)
end

export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
