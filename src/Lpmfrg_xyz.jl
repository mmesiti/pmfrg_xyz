module pmfrg_xyz

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

# Include common code shared with PMFRG_xyz
include("common.jl")

struct NumericalParams{T<:Real}
    T::T
    N::Int

    accuracy::T
    lambda_min::T
    lambda_max::T

    lenIntw::Int
    lenIntw_acc::Int
end

_getFloatType(Par) = typeof(Par.NumericalParams.T)

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

function get_iS(FlowParam::Real, iSigma::SigmaType, NumParams::NumericalParams)
    Lam = FlowParam
    T = NumParams.T

    @inline iSx(x, nw) = iS_(iSigma.x, x, Lam, nw, T) / 2
    @inline iSy(x, nw) = iS_(iSigma.y, x, Lam, nw, T) / 2
    @inline iSz(x, nw) = iS_(iSigma.z, x, Lam, nw, T) / 2

    return iSx, iSy, iSz

end


function get_iG_i(FlowParam::Real, iSigma_i::AbstractArray, NumParams::NumericalParams)
    Lam = FlowParam
    T = NumParams.T
    @inline iG_i(x, nw) = iG_(iSigma_i, x, Lam, nw, T)
    return iG_i
end

function get_iSKat(iSigma, DiSigma, FlowParam::Real, NumParams::NumericalParams)

    T = NumParams.T
    iSKat_x(x, nw) = iSKat_(iSigma.x, DiSigma.x, x, FlowParam, nw, T)
    iSKat_y(x, nw) = iSKat_(iSigma.y, DiSigma.y, x, FlowParam, nw, T)
    iSKat_z(x, nw) = iSKat_(iSigma.z, DiSigma.z, x, FlowParam, nw, T)

    return iSKat_x, iSKat_y, iSKat_z

end

function get_Theta(FlowParam::Real, _::NumericalParams)
    Lam = FlowParam
    Theta(w) = w^2 / (w^2 + Lam^2)
    return Theta
end

function get_f_int_factor(NumParams::NumericalParams)
    return NumParams.T
end

function get_get_w(NumParams::NumericalParams)
    return nw -> get_w(nw, NumParams.T)
end

function get_Dgamma_factor(NumParams::NumericalParams)
    return NumParams.T
end
function get_chi_factor(NumParams::NumericalParams)
    return NumParams.T
end

function get_props_factor(NumParams::NumericalParams)
    return NumParams.T
end


####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################


function flow_parameter_max_min(NumParams::NumericalParams)

    (; lambda_max, lambda_min) = NumParams
    return lambda_max, lambda_min
end



#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################


export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
