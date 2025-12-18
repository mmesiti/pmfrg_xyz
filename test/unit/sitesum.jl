using Test
using PMFRG_xyz

function test_site_sum_split()

    (; siteSum, Npairs, Nsum) = getSquareLattice(16)

    sitesum_split = PMFRG_xyz.split_sitesum(siteSum, 16, Npairs, Nsum)

    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m



    @testset "sitesum split blocking" begin
        for (; S, Nsum_split, Rijmin, Rijmax, kimin, kimax) in sitesum_split,
            Rij = Rijmin:Rijmax

            for ki = kimin:kimax

                Rij_iblock = Rij - Rijmin + 1
                ki_iblock = ki - kimin + 1

                for k_spl = 1:Nsum_split[ki_iblock, Rij_iblock]
                    (; kj, m, xk) = S[k_spl, ki_iblock, Rij_iblock]

                    for k_spl_old = 1:Nsum[Rij]
                        ki_old = S_ki[k_spl_old, Rij]
                        kj_old = S_kj[k_spl_old, Rij]
                        xk_old = S_xk[k_spl_old, Rij]
                        m_old = S_m[k_spl_old, Rij]
                        if ki == ki_old
                            if kj == kj_old
                                @test m_old == m
                                @test xk_old == xk
                            end
                        end
                    end

                end
            end
        end

        Nsum_check = zeros(Npairs)
        for (; Rijmin, Rijmax, Nsum_split) in sitesum_split
            Rijextent = Rijmax - Rijmin + 1
            Nsum_check[Rijmin:Rijmax] .+=
                [sum(Nsum_split[:, Rij_iblock]) for Rij_iblock = 1:Rijextent]
        end
        @test Nsum_check == Nsum
    end

end

function test_site_sum_split_v2_noblocking()
    (; siteSum, Npairs, Nsum) = getSquareLattice(16)

    unrolled_kj_m_xk, cumsum_ns = PMFRG_xyz.split_sitesum_v2_noblocking(siteSum, Npairs)

    start(Rij, ki) =
        let i = ki + (Rij - 1) * Npairs
            if i > 1
                cumsum_ns[i-1] + 1
            else
                1
            end
        end
    stop(Rij, ki) = cumsum_ns[ki+(Rij-1)*Npairs]

    @testset "sitesum split v2 no blocking " begin
        for Rij = 1:Npairs
            for ki = 1:Npairs
                for i = start(Rij, ki):stop(Rij, ki)
                    (kj, m, xk) = unrolled_kj_m_xk[i]
                    for k_spl = 1:Nsum[Rij]
                        if (siteSum[k_spl, Rij].ki == ki && siteSum[k_spl, Rij].kj == kj)
                            @test m == siteSum[k_spl, Rij].m
                            @test xk == siteSum[k_spl, Rij].xk
                        end
                    end
                end
            end
        end
        for Rij = 1:Npairs
            @test Nsum[Rij] == stop(Rij, Npairs) - start(Rij, 1) + 1
        end
    end
end
