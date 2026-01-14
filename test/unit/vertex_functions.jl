using Test
using SpinFRGLattices

import PMFRG_xyz: V_, FillVBuffer!, ConvertFreqArgs, fd

function test_fillvbuffer_v_equivalence()
    N = 8
    Npairs = 6
    Gamma = randn(21, Npairs, N, N, N)
    invpairs = [mod1(Npairs - i + 1, Npairs) for i = 1:Npairs]

    @testset for R = 1:Npairs,
        raw_s in rand(-N:N, 10),
        raw_t in rand(-N:N, 10),
        raw_u in rand(-N:N, 10),
        ft1 in [false, true],
        ft2 in [false, true],
        ft3 in [false, true]

        flavTransf = (ft1, ft2, ft3)
        s, t, u = ConvertFreqArgs(raw_s, raw_t, raw_u, N)
        R_swapped = flavTransf[1] ? invpairs[R] : R

        V_from_Vert = zeros(21)
        FillVBuffer!(V_from_Vert, (@view Gamma[:, R_swapped, s+1, t+1, u+1]), flavTransf)

        for n = 1:21
            V_from_V_ = V_(Gamma, n, raw_s, raw_t, raw_u, flavTransf, R, invpairs[R], N)
            @test V_from_Vert[n] == V_from_V_
        end
    end
end
