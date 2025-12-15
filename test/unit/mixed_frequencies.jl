using Test

using .PMFRG_xyz:
    mixedFrequenciesConverted,
    mixedFrequenciesConvertedInverse,
    mixedFrequenciesX,
    ConvertFreqArgs

"""
Property test for mixedFrequenciesConverted.
Tests that the function correctly combines mixedFrequenciesX, ConvertFreqArgs, and flavor transformations.
"""
function test_mixedFrequenciesConverted()
    @testset "mixedFrequenciesConverted property tests" begin
        # Test with various N values
        for N in [8, 10, 16, 24]
            @testset "N=$N" begin
                # Test with various frequency combinations
                for ns = 0:N-1, nt = 0:N-1, nu = 0:N-1
                    # Only test valid frequency combinations (is+it+iu must be odd)
                    if (ns + nt + nu) % 2 == 0
                        continue
                    end

                    for nwpr = -N:N-1
                        # Call the new combined function
                        t1, u1, t2, u2, flavTransf12, flavTransf34 =
                            mixedFrequenciesConverted(ns, nt, nu, nwpr, N)

                        # Verify against existing code
                        wpw1, wpw2, wmw3, wmw4 = mixedFrequenciesX(ns, nt, nu, nwpr)

                        # Check flavor transformations
                        expected_flavTransf12 =
                            (wpw1 * wpw2 > 0, ns * wpw2 > 0, ns * wpw1 < 0)
                        expected_flavTransf34 =
                            (wmw3 * wmw4 < 0, ns * wmw4 > 0, ns * wmw3 > 0)
                        @test flavTransf12 == expected_flavTransf12
                        @test flavTransf34 == expected_flavTransf34

                        # Check converted frequency arguments
                        expected_s1, expected_t1, expected_u1 =
                            ConvertFreqArgs(ns, wpw1, -wpw2, N)
                        expected_s2, expected_t2, expected_u2 =
                            ConvertFreqArgs(ns, -wmw3, -wmw4, N)

                        @test t1 == expected_t1
                        @test u1 == expected_u1
                        @test t2 == expected_t2
                        @test u2 == expected_u2
                    end
                end
            end
        end
    end
end

"""
Property test for mixedFrequenciesConvertedInverse.
Tests round-trip: (ns,nt,nu,nwpr,N) -> mixedFrequenciesConverted -> mixedFrequenciesConvertedInverse
The inverse returns all possible inputs that map to the same output.
"""
function test_mixedFrequenciesConvertedInverse()
    @testset "mixedFrequenciesConvertedInverse round-trip tests" begin
        @testset for N in (8, 16, 24, 32)
            # Test with full domain cases (including ns = 0)
            # Note: ns+nt+nu must be odd
            test_cases = [
                (ns, nt, nu, nwpr) for ns = 0:N-1 for nt = 0:N-1 for nu = 0:N-1 for
                nwpr = -N:N-1 if (ns + nt + nu) % 2 == 1
            ]

            @testset "N=$N" begin
                for (ns, nt, nu, nwpr) in test_cases
                    # Forward direction
                    t1, u1, t2, u2, flavTransf12, flavTransf34 =
                        mixedFrequenciesConverted(ns, nt, nu, nwpr, N)

                    # Inverse direction - returns all possible inputs
                    recovered_inputs = mixedFrequenciesConvertedInverse(
                        ns,
                        t1,
                        u1,
                        t2,
                        u2,
                        flavTransf12,
                        flavTransf34,
                        N,
                    )

                    # Test that the original input is in the list of recovered inputs
                    @test (nt, nu, nwpr) in recovered_inputs

                    # Verify all returned inputs actually produce the same output
                    for (nt_check, nu_check, nwpr_check) in recovered_inputs
                        t1_check, u1_check, t2_check, u2_check, flav12_check, flav34_check =
                            mixedFrequenciesConverted(ns, nt_check, nu_check, nwpr_check, N)

                        @test t1_check == t1
                        @test u1_check == u1
                        @test t2_check == t2
                        @test u2_check == u2
                        @test flav12_check == flavTransf12
                        @test flav34_check == flavTransf34
                    end
                end
            end
        end
    end
end
