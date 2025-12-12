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
                        s, t1, u1, t2, u2, flavTransf12, flavTransf34 =
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

                        @test s == expected_s1
                        @test s == expected_s2  # s1 == s2 always
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
Tests round-trip: (ns,nt,nu,nwpr,N) -> mixedFrequenciesConverted -> mixedFrequenciesConvertedInverse -> (ns,nt,nu,nwpr,N)
"""
function test_mixedFrequenciesConvertedInverse()
    @testset "mixedFrequenciesConvertedInverse round-trip tests" begin
        N = 16
        # Test with a few specific cases with ns strictly positive
        # Note: ns+nt+nu must be odd
        test_cases = [
            (3, 5, 7, 2),    # 3+5+7 = 15 (odd)
            (1, 2, 4, -5),   # 1+2+4 = 7 (odd)
            (7, 1, 3, 3),    # 7+1+3 = 11 (odd)
            (5, 4, 2, -2),   # 5+4+2 = 11 (odd)
            (2, 1, 4, 1),    # 2+1+4 = 7 (odd)
        ]

        @testset "N=$N" begin
            for (ns, nt, nu, nwpr) in test_cases
                # Forward direction
                s, t1, u1, t2, u2, flavTransf12, flavTransf34 =
                    mixedFrequenciesConverted(ns, nt, nu, nwpr, N)

                # Inverse direction
                ns_recovered, nt_recovered, nu_recovered, nwpr_recovered =
                    mixedFrequenciesConvertedInverse(
                        s,
                        t1,
                        u1,
                        t2,
                        u2,
                        flavTransf12,
                        flavTransf34,
                        N,
                    )

                # Test round-trip
                @test ns_recovered == ns
                @test nt_recovered == nt
                @test nu_recovered == nu
                @test nwpr_recovered == nwpr
            end
        end
    end
end
