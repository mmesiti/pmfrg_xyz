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


function get_iG_i(FlowParam::Real, iSigma_i::AbstractArray, ::NumericalParams)
    T = FlowParam
    @inline iG_i(x, nw) = iG_(iSigma_i, x, nw, T)
    return iG_i
end

function get_iSKat(iSigma, DiSigma, FlowParam::Real, _::NumericalParams)

    iSKat_x(x, nw) = iSKat_(iSigma.x, DiSigma.x, x, nw, FlowParam)
    iSKat_y(x, nw) = iSKat_(iSigma.y, DiSigma.y, x, nw, FlowParam)
    iSKat_z(x, nw) = iSKat_(iSigma.z, DiSigma.z, x, nw, FlowParam)

    return iSKat_x, iSKat_y, iSKat_z

end



function get_Theta(_::Real, ::NumericalParams)
    return _ -> 1
end

function get_f_int_factor(::NumericalParams)
    return 1
end

function get_get_w(::NumericalParams)
    return nw -> get_w(nw)
end

function get_Dgamma_factor(::NumericalParams)
    return 1
end

function get_chi_factor(::NumericalParams)
    return 1
end

function get_props_factor(::NumericalParams)
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


export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
