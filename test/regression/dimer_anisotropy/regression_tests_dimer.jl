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


function convert_X_Xnew(X)
    nf, npairs, N, _, _ = size(X)
    Xnew = zeros(nf, npairs, N, N, N)
    for is = 1:N, it = 1:N, iu = 1:N
        Xnew[:, :, iu, it, is] = X[:, :, is, it, iu]
    end
    Xnew
end

function run_getXbubble_regression_tests()
    @testset verbose = true "Tests for PMFRG_xyz.getXBubble!" begin
        data = load_object(
            joinpath(@__DIR__(), "regression_tests_dimer-PMFRG_xyz.getXBubble!.data"),
        )
        @testset verbose = false for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            Workspace, T = arguments = (data["arguments"])[i]

            Workspace_post, T_post = (data["arguments_post"])[i]

            Workspace.X[:, :, :, :, :] = convert_X_Xnew(Workspace.X)
            Workspace.State.Gamma[:, :, :, :, :] = convert_X_Xnew(Workspace.State.Gamma)


            recursive_value_test(
                return_value,
                (PMFRG_xyz.getXBubble!)(Workspace, T),
                "return values - case $i",
                true,
            )

            Workspace_post.X[:, :, :, :, :] = convert_X_Xnew(Workspace_post.X)
            Workspace_post.State.Gamma[:, :, :, :, :] =
                convert_X_Xnew(Workspace_post.State.Gamma)


            @testset "arguments" begin
                for i in eachindex(arguments)
                    recursive_value_test(Workspace, Workspace_post, "idx = ws", false)
                    recursive_value_test(T, T_post, "idx = T", false)
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
            exp_return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            return_value = (PMFRG_xyz.SolveFRG)(arguments...)

            begin
                let expsol = exp_return_value[1]
                    for thing in (expsol.prob.p.X, expsol.u[2].x[5])
                        thing[:, :, :, :, :] = convert_X_Xnew(thing)
                    end
                end
            end

            recursive_value_test(
                exp_return_value,
                return_value,
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
