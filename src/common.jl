# Common code between Lambda Flow (Lpmfrg_xyz) and T Flow (PMFRG_xyz)
# This file is included by both modules
# Note: Each module handles its own using statements

#################################################
######### UTILITIES #############################
#################################################

using SpinFRGLattices, OrdinaryDiffEqLowOrderRK, DiffEqCallbacks, StructArrays
using RecursiveArrayTools
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

#################################################
######### STRUCTS ###############################
#################################################

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

struct OptionParams
    use_symmetry::Bool
    minimal_output::Bool
end

#################################################
######### DIMENSIONS AND CONSTRUCTORS ###########
#################################################

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

function SigmaType(NUnique::Int, N::Int, type = Float64)
    return SigmaType(
        zeros(type, NUnique, N),
        zeros(type, NUnique, N),
        zeros(type, NUnique, N),
    )
end
SigmaType(Par) = SigmaType(Par.System.Npairs, Par.NumericalParams.N)

function StateType(NUnique::Int, N::Int, VDims::Tuple, type = Float64)
    return StateType(zeros(type, NUnique), SigmaType(type, NUnique, N), zeros(type, VDims))
end
StateType(Par) =
    StateType(Par.System.NUnique, Par.NumericalParams.N, getVDims(Par), _getFloatType(Par))
StateType(f_int, iSigma_x, iSigma_y, iSigma_z, Gamma) =
    StateType(f_int, SigmaType(iSigma_x, iSigma_y, iSigma_z), Gamma)
RecursiveArrayTools.ArrayPartition(x) =
    ArrayPartition(x.f_int, x.iSigma.x, x.iSigma.y, x.iSigma.z, x.Gamma)
StateType(Arr::ArrayPartition) = StateType(Arr.x...)

OptionParams(; use_symmetry::Bool = true, MinimalOutput::Bool = false, kwargs...) =
    OptionParams(use_symmetry, MinimalOutput)
Params(System; kwargs...) =
    OneLoopParams(System, NumericalParams(; kwargs...), OptionParams(; kwargs...))

#################################################
######### PROPAGATOR HELPERS ####################
#################################################

function get_sign_iw(nw::Integer, N::Integer)
    # s = sign(nw)
    nw_bounds = min(nw, N - 1)
    return nw_bounds + 1
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

function ConvertFreqArgs(ns, nt, nu, Nw)
    ns, nt, nu = abs.((ns, nt, nu))

    ns = min(ns, Nw - 1 - (ns + Nw - 1) % 2) ### weird cutoff, idk why
    nt = min(nt, Nw - 1 - (nt + Nw - 1) % 2)
    nu = min(nu, Nw - 1 - (nu + Nw - 1) % 2)

    return ns, nt, nu
end

#################################################
######### VERTEX HELPERS ########################
#################################################

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

# Flavor dimension constants
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

#################################################
######### VERTEX ACCESS #########################
#################################################

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
        block = div(n + 2, 6)

        if (block != 0 && isFlavorTransform[block])
            block_start = 4 + (block - 1) * 6
            offset = n - block_start
            new_offset = (offset + 3) % 6

            n_transf = block_start + new_offset
        else
            n_transf = n
        end

        ns, nt, nu = ConvertFreqArgs(ns, nt, nu, N)
        Rij = ifelse(isFlavorTransform[1], Rji, Rij)
        return Vertex[n_transf, Rij, ns+1, nt+1, nu+1]
    end
end

#################################################
######### SYMMETRY FUNCTIONS ####################
#################################################

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

#################################################
######### 1-PARTICLE BUBBLE #####################
#################################################

function compute1PartBubble!(Dgamma::SigmaType, Gamma::Array{T,5}, Props, Par) where {T}
    invpairs = Par.System.invpairs

    setZero!(Dgamma)
    @inline Gamma_(n, Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, n, s, t, u, isFlavorTransform, Rij, invpairs[Rij], Par.NumericalParams.N)
    addTo1PartBubble!(Dgamma, Gamma_, Props, Par)
end

#################################################
######### FLOW PARAMETER FUNCTIONS ##############
#################################################

t_to_Lam(t) = exp(t)
Lam_to_t(t) = log(t)

function gettMesh(T_min, T_max, npoints)
    t_min = get_t_min(T_min)
    t_max = Lam_to_t(T_max)
    return LinRange(t_min, t_max, npoints)
end

function get_t_min(Lam)
    Lam < exp(-30) && @warn "Flow parameter minimum too small! Set to exp(-30) instead."
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

#################################################
######### INITIALIZATION FUNCTIONS ##############
#################################################

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

struct OneLoopParams{SType,NumericalParams}
    System::SType
    NumericalParams::NumericalParams
    Options::OptionParams
end



function AllocateSetup(Par::OneLoopParams)
    println("Allocate Setup")
    ## Allocate Memory:
    floattype = _getFloatType(Par)
    return (X = zeros(floattype, getBubbleVDims(Par)), Par = Par)
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
