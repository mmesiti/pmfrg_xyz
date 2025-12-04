using Test

include("../../src/Xtype.jl")

"""
    test_mapping_bijectivity(mapping::AbstractXIndexMapping)

Test that the mapping between linear and multi-indices is bijective.
This is the core requirement for any valid index mapping.
"""
function test_mapping_bijectivity(mapping::AbstractXIndexMapping)
    n_flavors = mapping.n_flavors
    n_pairs = mapping.n_pairs
    N = mapping.N

    @testset "Bijectivity tests" begin
        # Test 1: All linear indices map to valid multi-indices
        total_len = total_length(mapping)
        @testset "All linear indices valid" begin
            for i = 1:total_len
                n, Rij, is, it, iu = multi_index(mapping, i)

                # Check bounds
                @test 1 <= n <= n_flavors
                @test 1 <= Rij <= n_pairs
                @test 1 <= is <= N
                @test 1 <= it <= N
                @test 1 <= iu <= N

                # Check validity condition
                @test is_valid_multi_index(is, it, iu)
            end
        end

        # Test 2: Round-trip linear -> multi -> linear
        @testset "Round-trip: linear -> multi -> linear" begin
            for i = 1:total_len
                n, Rij, is, it, iu = multi_index(mapping, i)
                i_roundtrip = linear_index(mapping, n, Rij, is, it, iu)
                @test i_roundtrip == i
            end
        end

        # Test 3: Round-trip multi -> linear -> multi for all valid combinations
        @testset "Round-trip: multi -> linear -> multi" begin
            for iu = 1:N, it = 1:N, is = 1:N
                if !is_valid_multi_index(is, it, iu)
                    continue
                end

                for Rij = 1:n_pairs, n = 1:n_flavors
                    i = linear_index(mapping, n, Rij, is, it, iu)
                    n_rt, Rij_rt, is_rt, it_rt, iu_rt = multi_index(mapping, i)

                    @test n_rt == n
                    @test Rij_rt == Rij
                    @test is_rt == is
                    @test it_rt == it
                    @test iu_rt == iu
                end
            end
        end

        # Test 4: All linear indices are unique (injectivity)
        @testset "All linear indices unique" begin
            seen_indices = Set{Int}()

            for iu = 1:N, it = 1:N, is = 1:N
                if !is_valid_multi_index(is, it, iu)
                    continue
                end

                for Rij = 1:n_pairs, n = 1:n_flavors
                    i = linear_index(mapping, n, Rij, is, it, iu)

                    # Check uniqueness
                    @test !(i in seen_indices)
                    push!(seen_indices, i)

                    # Check bounds
                    @test 1 <= i <= total_len
                end
            end

            # Check we found all indices
            @test length(seen_indices) == total_len
        end

        # Test 5: Verify total_length is correct
        @testset "Total length correct" begin
            # Count valid (is, it, iu) combinations
            num_valid_freq = 0
            for iu = 1:N, it = 1:N, is = 1:N
                if is_valid_multi_index(is, it, iu)
                    num_valid_freq += 1
                end
            end

            expected_length = num_valid_freq * n_flavors * n_pairs
            @test total_length(mapping) == expected_length
        end
    end
end

"""
    test_xvector_operations(mapping::AbstractXIndexMapping)

Test XVector operations including access, modification, and array interface.
"""
function test_xvector_operations(mapping::AbstractXIndexMapping)
    n_flavors = mapping.n_flavors
    n_pairs = mapping.n_pairs
    N = mapping.N

    @testset "XVector operations" begin
        # Create XVector
        xvec = XVector{Float64}(mapping)

        @testset "Basic properties" begin
            @test length(xvec) == total_length(mapping)
            @test eltype(xvec) == Float64
            @test all(xvec.data .== 0.0)
        end

        @testset "Multi-index access" begin
            # Set some values
            for iu = 1:N, it = 1:N, is = 1:N
                if !is_valid_multi_index(is, it, iu)
                    continue
                end

                for Rij = 1:min(2, n_pairs), n = 1:min(3, n_flavors)
                    val = Float64(n + 10 * Rij + 100 * is + 1000 * it + 10000 * iu)
                    xvec[n, Rij, is, it, iu] = val
                end
            end

            # Check values
            for iu = 1:N, it = 1:N, is = 1:N
                if !is_valid_multi_index(is, it, iu)
                    continue
                end

                for Rij = 1:min(2, n_pairs), n = 1:min(3, n_flavors)
                    expected = Float64(n + 10 * Rij + 100 * is + 1000 * it + 10000 * iu)
                    @test xvec[n, Rij, is, it, iu] == expected
                end
            end
        end

        @testset "Fill operations" begin
            fill!(xvec, 42.0)
            @test all(xvec.data .== 42.0)

            for iu = 1:N, it = 1:N, is = 1:N
                if !is_valid_multi_index(is, it, iu)
                    continue
                end
                @test xvec[1, 1, is, it, iu] == 42.0
            end
        end

        @testset "Similar and zero" begin
            xvec2 = similar(xvec)
            @test typeof(xvec2) == typeof(xvec)
            @test length(xvec2) == length(xvec)

            xvec3 = zero(xvec)
            @test all(xvec3.data .== 0.0)
        end
    end
end

"""
    test_validity_condition()

Test the is_valid_multi_index function.
"""
function test_validity_condition()
    @testset "Validity condition" begin
        # Valid cases: (is-1) + (it-1) + (iu-1) is odd
        @test is_valid_multi_index(1, 1, 2)  # 0+0+1 = 1 (odd)
        @test is_valid_multi_index(2, 2, 2)  # 1+1+1 = 3 (odd)
        @test is_valid_multi_index(1, 2, 1)  # 0+1+0 = 1 (odd)

        # Invalid cases: (is-1) + (it-1) + (iu-1) is even
        @test !is_valid_multi_index(1, 1, 1)  # 0+0+0 = 0 (even)
        @test !is_valid_multi_index(2, 3, 4)  # 1+2+3 = 6 (even)
        @test !is_valid_multi_index(1, 2, 2)  # 0+1+1 = 2 (even)
    end
end

"""
    run_all_tests(mapping::AbstractXIndexMapping)

Run all tests for a given mapping implementation.
"""
function run_all_tests(mapping::AbstractXIndexMapping)
    @testset "XType tests for $(typeof(mapping))" begin
        test_validity_condition()
        test_mapping_bijectivity(mapping)
        test_xvector_operations(mapping)
    end
end

# Export test functions
export test_mapping_bijectivity,
    test_xvector_operations, test_validity_condition, run_all_tests
