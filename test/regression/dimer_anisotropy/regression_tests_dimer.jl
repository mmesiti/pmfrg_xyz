using PMFRG_xyz
using Test
using JLD2
function recursive_value_equality(strA, strB)
    strA == strB ||
        typeof(strA) == typeof(strB) && if fieldnames(typeof(strA)) == ()
            strA == strB
        else
            all((
                getfield(strA, field) == getfield(strB, field) ||
                recursive_value_equality(getfield(strA, field), getfield(strB, field)) for
                field in fieldnames(typeof(strA))
            ))
        end
end
function recursive_value_equality(strA::Vector, strB::Vector)
    strA == strB ||
        typeof(strA) == typeof(strB) && (
            length(strA) == length(strB) && all((
                strA[i] == strB[i] || recursive_value_equality(strA[i], strB[i]) for
                i in eachindex(strA)
            ))
        )
end
function (recursive_value_equality(
    strA::Vector{T},
    strB::Vector{U},
) where {T<:Number}) where {U<:Number}
    length(strA) == length(strB) && (
        strA ≈ strB || all((
            strA[i] ≈ strB[i] || recursive_value_equality(strA[i], strB[i]) for
            i in eachindex(strA)
        ))
    )
end

function run_getXbubble_regression_tests()
        @testset verbose = true "Tests for PMFRG_xyz.getXBubble!" begin
            Core.@doc "You might need to modify this function!" function compare_return_values(
                rvexp,
                rv,
            )
                recursive_value_equality(rvexp, rv)
            end
            Core.@doc "You might need to modify this function!" function compare_arguments_post(
                args_post_exp,
                arg_post,
            )
                recursive_value_equality(args_post_exp, arg_post)
            end
            data = load_object(
                joinpath(@__DIR__(), "regression_tests_dimer-PMFRG_xyz.getXBubble!.data"),
            )
            @testset for i = 1:length(data["return_value"])
                return_value = (data["return_value"])[i]
                arguments = (data["arguments"])[i]
                arguments_post = (data["arguments_post"])[i]
                Workspace = arguments[1]
                ThreadLocalBuffers = PMFRG_xyz.get_ThreadLocalBuffers(Workspace.Par.System)
                @test compare_return_values(return_value, (PMFRG_xyz.getXBubble!)(arguments...)) &&
                      compare_arguments_post(arguments, arguments_post)
            end
        end
        end
 

function run_SolveFRG_regression_tests()
       @testset verbose = true "Tests for PMFRG_xyz.SolveFRG" begin
            Core.@doc "You might need to modify this function!" function compare_return_values(
                rvexp,
                rv,
            )
                recursive_value_equality(rvexp, rv)
            end
            Core.@doc "You might need to modify this function!" function compare_arguments_post(
                args_post_exp,
                arg_post,
            )
                recursive_value_equality(args_post_exp, arg_post)
            end
            data =
                load_object(joinpath(@__DIR__(), "regression_tests_dimer-PMFRG_xyz.SolveFRG.data"))
            @testset for i = 1:length(data["return_value"])
                return_value = (data["return_value"])[i]
                arguments = (data["arguments"])[i]
                arguments_post = (data["arguments_post"])[i]
                @test compare_return_values(return_value, (PMFRG_xyz.SolveFRG)(arguments...)) &&
                      compare_arguments_post(arguments, arguments_post)
            end
        end
    end

