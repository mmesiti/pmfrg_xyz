using JLD2

function get_data()
    (; State, FlowParam, Par) =
        load_object(joinpath(@__DIR__, "get_Chi_dimer_random_input.jld2"))
    (; chi_x, chi_y, chi_z) = load_object(joinpath(@__DIR__, "get_Chi_dimer_output.jld2"))
    return (; State, FlowParam, Par, chi_x, chi_y, chi_z)
end


function test_chi_Tflow()
    (; State, FlowParam, Par, chi_x, chi_y, chi_z) = get_data()

    @testset "T-Flow" begin
        @test getChi_z(State, FlowParam, Par) ≈ chi_z
        @test getChi_x(State, FlowParam, Par) ≈ chi_x
        @test getChi_y(State, FlowParam, Par) ≈ chi_y
    end
end
