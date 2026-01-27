# Generic regression test utilities
# These functions do not depend on any specific package (PMFRG_xyz or pmfrg_xyz)
# and can be used for both T Flow and Lambda Flow regression tests

using Test
using JLD2
import SciMLBase

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

# Skip ODEFunction comparisons (not relevant for regression tests)
recursive_value_test(_, ::SciMLBase.ODEFunction, _, _)::Nothing = nothing

function recursive_value_test(strA::Array, strB::Array, name, verbose)::Nothing
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

function recursive_value_test(
    strA::Array{T},
    strB::Array{U},
    _,
    _,
)::Nothing where {T<:Number} where {U<:Number}
    @testset verbose = true "lengths" begin
        @test length(strA) == length(strB)
    end
    if length(strA) == length(strB)
        @test strA ≈ strB
    end
    nothing
end
