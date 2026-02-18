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

struct ObservablesDynamic{T}
    Chi_x::Array{T,2}
    Chi_y::Array{T,2}
    Chi_z::Array{T,2}
end

struct NumericalParams{T<:Real}
    N::Int

    accuracy::T
    temp_min::T
    temp_max::T

    lenIntw::Int
    lenIntw_acc::Int
end

struct OptionParams
    use_symmetry::Bool
    minimal_output::Bool
end

struct OneLoopParams{T,SType}
    System::SType
    NumericalParams::NumericalParams{T}
    Options::OptionParams
end

struct OneLoopWorkspace{T,ParType}
    State::StateType{T}
    Deriv::StateType{T}
    X::Array{T,5}
    Par::ParType
end

# A general vertex can have 3^4 = 81 flavor combinations
# The XYZ model possesses Klein-4 Symmetry reducing the amount to 21.
# The two distinct bubbles X and Y are stored in one 42-dimensional array.
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
