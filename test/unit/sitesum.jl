
using PMFRG_xyz

function test_site_sum_split()
    (; siteSum, Npairs, Nsum) = getSquareLattice(16)

    sitesum_split = PMFRG_xyz.split_sitesum(siteSum, 16, Npairs, Nsum)

    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m



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
