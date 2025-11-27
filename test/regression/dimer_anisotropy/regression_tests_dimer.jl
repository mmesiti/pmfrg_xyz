using PMFRG_xyz
using Test
using JLD2
function recursive_value_test(strA, strB, name)::Nothing
    if strA == strB
        @test strA == strB # true
    else
        if fieldnames(typeof(strA)) == ()
            @test strA == strB # ralse
        else
            @testset "$name" for field in fieldnames(typeof(strA))

                recursive_value_test(getfield(strA, field),
                    getfield(strB, field),
                    "$field")
            end
        end
    end
    nothing
end

recursive_value_test(::Nothing, ::Nothing, name)::Nothing = nothing

function recursive_value_test(strA::Vector, strB::Vector, name)::Nothing
    if strA == strB
        @test strA == strB
    else
        @test length(strA) == length(strB)
        if length(strA) == length(strB)
            @testset "$name" for i in eachindex(strA)
                recursive_value_test(strA[i], strB[i], "$i")
            end
        end
    end
    nothing
end

function recursive_value_equality(
    strA::Vector{T},
    strB::Vector{U},
    _ = nothing)::Nothing where {T<:Number} where {U<:Number}
    @test length(strA) == length(strB)
    if length(strA) == length(strB)
        @test strA ≈ strB
    end
    nothing
end

function run_getXbubble_regression_tests()
    @testset verbose = true "Tests for PMFRG_xyz.getXBubble!" begin
        data = load_object(
            joinpath(@__DIR__(), "regression_tests_dimer-PMFRG_xyz.getXBubble!.data"),
        )
        @testset for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            # Workspace = arguments[1]
            # ThreadLocalBuffers = PMFRG_xyz.get_ThreadLocalBuffers(Workspace.Par.System)
            recursive_value_test(return_value,
                (PMFRG_xyz.getXBubble!)(arguments...),
                "return values - case $i")
            recursive_value_test(arguments,
                arguments_post,
                "arguments - case $i")
        end
    end
    nothing
end


function run_SolveFRG_regression_tests()
    @testset verbose = true "Tests for PMFRG_xyz.SolveFRG" begin
        data =
            load_object(joinpath(@__DIR__(), "regression_tests_dimer-PMFRG_xyz.SolveFRG.data"))
        @testset for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            recursive_value_test(return_value,
                (PMFRG_xyz.SolveFRG)(arguments...),
                "return values - case $i")
            recursive_value_test(arguments,
                arguments_post,
                "arguments - case $i")
        end
    end
end

