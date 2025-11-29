using PMFRG_xyz
using Test
using JLD2
function recursive_value_test(strA, strB, name, verbose)::Nothing
    @testset verbose = verbose "$name" begin
        if strA == strB
            @test strA == strB # true
        else
            if fieldnames(typeof(strA)) == ()
                @test strA == strB # false
            else
                for field in union(fieldnames(typeof(strA)), fieldnames(typeof(strB)))

                    if field in fieldnames(typeof(strA)) &&
                       field in fieldnames(typeof(strB))
                        recursive_value_test(
                            getfield(strA, field),
                            getfield(strB, field),
                            "$field",
                            false,
                        )
                    else
                        @test field in fieldnames(typeof(strA))
                        @test field in fieldnames(typeof(strB))
                    end
                end
            end
        end
    end
    nothing
end

recursive_value_test(::Nothing, ::Nothing, _, _)::Nothing = nothing

function recursive_value_test(strA::Vector, strB::Vector, name, verbose)::Nothing
    @testset verbose = verbose "$name" begin
        if strA == strB
            @test strA == strB
        else
            @test length(strA) == length(strB)
            if length(strA) == length(strB)
                for i in eachindex(strA)
                    recursive_value_test(strA[i], strB[i], "$i", false)
                end
            end
        end
    end
    nothing
end

function recursive_value_equality(
    strA::Vector{T},
    strB::Vector{U},
    _,
    _,
)::Nothing where {T<:Number} where {U<:Number}
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
        @testset verbose = false for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            # Workspace = arguments[1]
            # ThreadLocalBuffers = PMFRG_xyz.get_ThreadLocalBuffers(Workspace.Par.System)
            recursive_value_test(
                return_value,
                (PMFRG_xyz.getXBubble!)(arguments...),
                "return values - case $i",
                true,
            )
            @testset "arguments" begin
                for i in eachindex(arguments)
                    recursive_value_test(arguments[i], arguments_post[i], "idx = $i", false)
                end
            end

        end
    end
    nothing
end


function run_SolveFRG_regression_tests()
    @testset verbose = true "Tests for PMFRG_xyz.SolveFRG" begin
        data = load_object(
            joinpath(@__DIR__(), "regression_tests_dimer-PMFRG_xyz.SolveFRG.data"),
        )
        @testset verbose = true for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            recursive_value_test(
                return_value,
                (PMFRG_xyz.SolveFRG)(arguments...),
                "return values - case $i",
                true,
            )
            @testset verbose = true for i in eachindex(arguments)
                recursive_value_test(
                    arguments[i],
                    arguments_post[i],
                    "arguments - case $i",
                    false,
                )
            end
        end
    end
end
