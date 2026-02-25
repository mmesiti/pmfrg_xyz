module PMFRG_xyz



using SpinFRGLattices, OrdinaryDiffEqLowOrderRK, DiffEqCallbacks, StructArrays
using RecursiveArrayTools
using SpinFRGLattices.StaticArrays
using Unroll
using MuladdMacro
using FastBroadcast

#################################################
######### UTILITIES #############################
#################################################


_getFloatType(Par) = typeof(Par.NumericalParams.accuracy)

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

abstract type AbstractNumericalParams end
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

Params(System, NumParams::AbstractNumericalParams; kwargs...) =
    OneLoopParams(System, NumParams, OptionParams(; kwargs...))


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

function get_Self_Energy!(Workspace, FlowParam::Real)
    Par = Workspace.Par
    (; iSigma, Gamma) = Workspace.State
    DiSigma = Workspace.Deriv.iSigma

    props = get_iS(FlowParam, iSigma, Par.NumericalParams)
    compute1PartBubble!(DiSigma, Gamma, props, Par)

end

function getDFint!(Workspace, FlowParam::Real)
    (; State, Deriv, Par) = Workspace
    (; lenIntw_acc) = Par.NumericalParams
    NUnique = Par.System.NUnique

    iSigmax(x, nw) = iSigma_(State.iSigma.x, x, nw)
    iSigmay(x, nw) = iSigma_(State.iSigma.y, x, nw)
    iSigmaz(x, nw) = iSigma_(State.iSigma.z, x, nw)

    iGx, iGy, iGz = get_iGs(FlowParam, State.iSigma, Par.NumericalParams)

    iSx, iSy, iSz = get_iS(FlowParam, State.iSigma, Par.NumericalParams)

    Theta = get_Theta(FlowParam, Par.NumericalParams)

    f = T_Dimension(Par.NumericalParams)

    _get_w = get_get_w(Par.NumericalParams)

    for x = 1:NUnique
        sumres = 0.0
        for nw = (-lenIntw_acc):(lenIntw_acc-1)
            w = _get_w(nw)

            sumres += iSx(x, nw) / iGx(x, nw) * Theta(w) * iSigmax(x, nw) / w
            sumres += iSy(x, nw) / iGy(x, nw) * Theta(w) * iSigmay(x, nw) / w
            sumres += iSz(x, nw) / iGz(x, nw) * Theta(w) * iSigmaz(x, nw) / w
        end
        Deriv.f_int[x] = -f * sumres
    end
end
function addTo1PartBubble!(Dgamma::SigmaType, Gamma_::Function, Props, Par)

    (; N, lenIntw_acc) = Par.NumericalParams
    (; siteSum, Nsum, OnsitePairs) = Par.System

    f = T_Dimension(Par.NumericalParams)

    Threads.@threads for iw1 = 1:N
        nw1 = iw1 - 1
        for (x, Rx) in enumerate(OnsitePairs)
            for nw = (-lenIntw_acc):(lenIntw_acc-1)
                jsum = zeros(3)
                wpw1 = nw1 + nw + 1
                wmw1 = nw - nw1
                flavTransform = (wmw1 * wpw1 < 0, false, false)
                for k_spl = 1:Nsum[Rx]
                    (; m, ki, xk) = siteSum[k_spl, Rx]
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
                Dgamma.x[x, iw1] += -f * jsum[1]
                Dgamma.y[x, iw1] += -f * jsum[2]
                Dgamma.z[x, iw1] += -f * jsum[3]
            end
        end
    end
end
using JLD2
function getDeriv!(Deriv, State, setup, FlowParameter; saveArgs = true)

    (; X, Par) = setup # use pre-allocated X and XTilde to reduce garbage collector time

    Workspace = OneLoopWorkspace(State, Deriv, X, Par)

    getDFint!(Workspace, FlowParameter)
    get_Self_Energy!(Workspace, FlowParameter)
    getXBubble!(Workspace, FlowParameter)

    symmetrizeBubble!(Workspace.X, Par)
    addToVertexFromBubble!(Workspace.Deriv.Gamma, Workspace.X)
    symmetrizeVertex!(Workspace.Deriv.Gamma, Par)

    return
end

####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################



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
    (; accuracy) = Par.NumericalParams
    flow_parameter_max, flow_parameter_min = flow_parameter_max_min(Par.NumericalParams)

    t0 = Lam_to_t(flow_parameter_max)
    tend = get_t_min(flow_parameter_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    saved_values = SavedValues(eltype(State), Observables{eltype(State)})

    function save_func(State, t, integrator)
        chi_x = getChi_x(State, t_to_Lam(t), Par)
        chi_y = getChi_y(State, t_to_Lam(t), Par)
        chi_z = getChi_z(State, t_to_Lam(t), Par)

        return Observables(copy(chi_x), copy(chi_y), copy(chi_z))
    end

    ObsSaveat = gettMesh(flow_parameter_min, flow_parameter_max, npoints)
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
        dt = Lam_to_t(0.2 * flow_parameter_max),
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

#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################


getChi_z(State::ArrayPartition, T::Real, Par) = getChi_3(
    State.x[2], # Sigma x
    State.x[3], # Sitma y
    State.x[5],
    T,
    fd.xy2,
    Par,
)
getChi_x(State::ArrayPartition, T::Real, Par) = getChi_3(
    State.x[3], # Sigma y
    State.x[4], # Sigma z
    State.x[5],
    T,
    fd.yz2,
    Par,
)
getChi_y(State::ArrayPartition, T::Real, Par) = getChi_3(
    State.x[4], # Sigma z
    State.x[2], # Sigma x
    State.x[5],
    T,
    fd.zx2,
    Par,
)

function getChi_3(
    iSigma1::AbstractArray,
    iSigma2::AbstractArray,
    Gamma::AbstractArray,
    FlowParam::Real,
    fd_idx,
    Par,
)
    (; N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iG1 = get_iG_i(FlowParam, iSigma1, Par.NumericalParams)
    iG2 = get_iG_i(FlowParam, iSigma2, Par.NumericalParams)
    V12_2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, fd_idx, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)
    f = T_Dimension(Par.NumericalParams)
    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = (-lenIntw_acc):(lenIntw_acc-1)
            if Rij in OnsitePairs
                Chi[Rij, 1] += f * iG1(xi, nK) * iG2(xi, nK)
            end
            for nK2 = (-lenIntw_acc):(lenIntw_acc-1)
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iG1(xi, nK)^2 * iG2(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += f^2 * GGGG * V12_2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end

function get_iGs(FlowParam::Real, iSigma::SigmaType, NumPar::AbstractNumericalParams)

    iGx = get_iG_i(FlowParam, iSigma.x, NumPar)
    iGy = get_iG_i(FlowParam, iSigma.y, NumPar)
    iGz = get_iG_i(FlowParam, iSigma.z, NumPar)

    return iGx, iGy, iGz

end
using LinearAlgebra
using SparseArrays
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


function set_spropX!(
    spropX,
    NUnique,
    iSigma,
    DiSigma,
    FlowParam,
    nw1,
    nw2,
    ComputeType,
    NumPar::AbstractNumericalParams,
)

    iGs = get_iGs(FlowParam, iSigma, NumPar)
    iSKat = get_iSKat(iSigma, DiSigma, FlowParam, NumPar)

    f = T_Dimension(NumPar)
    for Rij = 1:NUnique
        for j = 1:3, i = 1:3
            spropX[i, j, Rij] = ComputeType(-iSKat[i](Rij, nw1) * iGs[j](Rij, nw2) * f)
        end
    end

end


function set_spropY!(
    spropY,
    NUnique,
    iSigma,
    DiSigma,
    FlowParam,
    nw1,
    nw2,
    ComputeType,
    NumPar::AbstractNumericalParams,
)

    iGs = get_iGs(FlowParam, iSigma, NumPar)
    iSKat = get_iSKat(iSigma, DiSigma, FlowParam, NumPar)

    f = T_Dimension(NumPar)
    for Rij1 = 1:NUnique, Rij2 = 1:NUnique
        for j = 1:3, i = 1:3
            spropY[i, j, Rij1, Rij2] =
                ComputeType(-iSKat[i](Rij1, nw1) * iGs[j](Rij2, nw2) * f)
        end
    end


end

function getXBubble!(
    Workspace::OneLoopWorkspace,
    FlowParameter::Real;
    ComputeType::Type = Float64,
)
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

    Threads.@threads :static for is_it = 1:(N*N)
        @inbounds begin
            is = (is_it - 1) ÷ N + 1
            it = (is_it - 1) % N + 1
            # WARNING:
            # This works only with :static
            Buffs = ThreadLocalBuffers[Threads.threadid()]
            ns = is - 1
            nt = it - 1

            for nw = (-lenIntw):(lenIntw-1) # Matsubara sum
                # Update Katanin propagators for current nw (convert to ComputeType)
                set_spropX!(
                    Buffs.spropX,
                    NUnique,
                    iSigma,
                    DiSigma,
                    FlowParameter,
                    nw,
                    nw + ns,
                    ComputeType,
                    Par.NumericalParams,
                )
                set_spropY!(
                    Buffs.spropY,
                    NUnique,
                    iSigma,
                    DiSigma,
                    FlowParameter,
                    nw,
                    nw - nt,
                    ComputeType,
                    Par.NumericalParams,
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


##############################
########### T-FLOW ###########
##############################
struct TFlowNumericalParams{T<:Real} <: AbstractNumericalParams
    N::Int

    accuracy::T
    temp_min::T
    temp_max::T

    lenIntw::Int
    lenIntw_acc::Int
end

function TFlowNumericalParams(;
    N::Integer = 24,
    accuracy = 1e-6,
    temp_min = exp(-10.0),
    temp_max = exp(10.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = 2 * maximum((N, lenIntw)),
)

    return TFlowNumericalParams(N, accuracy, temp_min, temp_max, lenIntw, lenIntw_acc)
end

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


function get_iS(FlowParam::Real, iSigma::SigmaType, _::TFlowNumericalParams)
    T = FlowParam

    @inline iSx(x, nw) = iS_(iSigma.x, x, nw, T) / 2
    @inline iSy(x, nw) = iS_(iSigma.y, x, nw, T) / 2
    @inline iSz(x, nw) = iS_(iSigma.z, x, nw, T) / 2

    return iSx, iSy, iSz
end


function get_iG_i(FlowParam::Real, iSigma_i::AbstractArray, ::TFlowNumericalParams)
    T = FlowParam
    @inline iG_i(x, nw) = iG_(iSigma_i, x, nw, T)
    return iG_i
end

function get_iSKat(iSigma, DiSigma, FlowParam::Real, _::TFlowNumericalParams)

    iSKat_x(x, nw) = iSKat_(iSigma.x, DiSigma.x, x, nw, FlowParam)
    iSKat_y(x, nw) = iSKat_(iSigma.y, DiSigma.y, x, nw, FlowParam)
    iSKat_z(x, nw) = iSKat_(iSigma.z, DiSigma.z, x, nw, FlowParam)

    return iSKat_x, iSKat_y, iSKat_z

end



function get_Theta(_::Real, ::TFlowNumericalParams)
    return _ -> 1
end
"""
The factor of T that is there in the Lambda-flow version
disappears in the T-flow version on dimensional grounds, See:

https://journals.aps.org/prb/pdf/10.1103/PhysRevB.109.195109

Sec II, after Eq. 3.
"""
function T_Dimension(::TFlowNumericalParams)
    return 1
end

function get_get_w(::TFlowNumericalParams)
    return nw -> get_w(nw)
end

function flow_parameter_max_min(NumParams::TFlowNumericalParams)

    (; temp_max, temp_min) = NumParams
    return temp_max, temp_min
end



##############################
########### L-FLOW ###########
##############################


struct LFlowNumericalParams{T<:Real} <: AbstractNumericalParams
    T::T
    N::Int

    accuracy::T
    lambda_min::T
    lambda_max::T

    lenIntw::Int
    lenIntw_acc::Int
end


function LFlowNumericalParams(;
    T::Real = 0.5,
    N::Integer = 24,
    accuracy = 1e-6,
    lambda_min = exp(-10.0),
    lambda_max = exp(10.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = 2 * maximum((N, lenIntw)),
)

    return LFlowNumericalParams(
        T,
        N,
        accuracy,
        lambda_min,
        lambda_max,
        lenIntw,
        lenIntw_acc,
    )
end


#############################################################
######### PROPAGATORS ## PROPAGATORS ## PROPAGATORS #########
#############################################################

### Propagators will depend on an additional flavor
### Instead of modifying the propagators, I will simply use them
### as V_, by doing iG_(iSigma.x, ...)

function get_w(nw, T)
    return pi * T * (2 * nw + 1)
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

######################################################################
######### FLOW EQUATIONS ## FLOW EQUATIONS ## FLOW EQUATIONS #########
######################################################################

function get_iS(FlowParam::Real, iSigma::SigmaType, NumParams::LFlowNumericalParams)
    Lam = FlowParam
    T = NumParams.T

    @inline iSx(x, nw) = iS_(iSigma.x, x, Lam, nw, T) / 2
    @inline iSy(x, nw) = iS_(iSigma.y, x, Lam, nw, T) / 2
    @inline iSz(x, nw) = iS_(iSigma.z, x, Lam, nw, T) / 2

    return iSx, iSy, iSz

end


function get_iG_i(FlowParam::Real, iSigma_i::AbstractArray, NumParams::LFlowNumericalParams)
    Lam = FlowParam
    T = NumParams.T
    @inline iG_i(x, nw) = iG_(iSigma_i, x, Lam, nw, T)
    return iG_i
end

function get_iSKat(iSigma, DiSigma, FlowParam::Real, NumParams::LFlowNumericalParams)

    T = NumParams.T
    iSKat_x(x, nw) = iSKat_(iSigma.x, DiSigma.x, x, FlowParam, nw, T)
    iSKat_y(x, nw) = iSKat_(iSigma.y, DiSigma.y, x, FlowParam, nw, T)
    iSKat_z(x, nw) = iSKat_(iSigma.z, DiSigma.z, x, FlowParam, nw, T)

    return iSKat_x, iSKat_y, iSKat_z

end

function get_Theta(FlowParam::Real, _::LFlowNumericalParams)
    Lam = FlowParam
    Theta(w) = w^2 / (w^2 + Lam^2)
    return Theta
end

function T_Dimension(NumParams::LFlowNumericalParams)
    return NumParams.T
end

function get_get_w(NumParams::LFlowNumericalParams)
    return nw -> get_w(nw, NumParams.T)
end



function flow_parameter_max_min(NumParams::LFlowNumericalParams)

    (; lambda_max, lambda_min) = NumParams
    return lambda_max, lambda_min
end





export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
