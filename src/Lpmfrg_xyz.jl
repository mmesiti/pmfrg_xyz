module pmfrg_xyz

#################################################
######### STRUCTS ## STRUCTS ## STRUCTS #########
#################################################

# Include common code shared with PMFRG_xyz
include("common.jl")

struct NumericalParams{T<:Real}
    T::T
    N::Int

    accuracy::T
    lambda_min::T
    lambda_max::T

    lenIntw::Int
    lenIntw_acc::Int
end

_getFloatType(Par) = typeof(Par.NumericalParams.T)

function NumericalParams(;
    T::Real = 0.5,
    N::Integer = 24,
    accuracy = 1e-6,
    lambda_min = exp(-10.0),
    lambda_max = exp(10.0),
    lenIntw::Int = N,
    lenIntw_acc::Int = 2 * maximum((N, lenIntw)),
)

    return NumericalParams(T, N, accuracy, lambda_min, lambda_max, lenIntw, lenIntw_acc)
end


#############################################################
######### PROPAGATORS ## PROPAGATORS ## PROPAGATORS #########
#############################################################

### Propagators will depend on an additional flavor
### Instead of modifying the propagators, I will simply use them
### as V_, by doing iG_(iSigma.x, ...)

function get_w(nw, T)
    return pi * T * (2 * nw + 1)
end

function iG_(iSigma::AbstractArray, x::Integer, Lam::Real, nw::Integer, T::Real)
    w = get_w(nw, T)
    return w / (w^2 + w * iSigma_(iSigma, x, nw) + Lam^2)
end

function iS_(iSigma::AbstractArray, x::Integer, Lam::Real, nw::Integer, T::Real)
    w = get_w(nw, T)
    return -iG_(iSigma, x, Lam, nw, T)^2 * 2 * Lam / w
end

function iSKat_(
    iSigma::AbstractArray,
    DSigma::AbstractArray,
    x::Integer,
    Lam::Real,
    nw::Integer,
    T::Real,
)
    w = get_w(nw, T)
    return -iG_(iSigma, x, Lam, nw, T)^2 * (2 * Lam / w + iSigma_(DSigma, x, nw))
end

####################################################
######### VERTICES ## VERTICES ## VERTICES #########
####################################################

### Symmetries: 
###     s <--> -s
###     t <--> -t, i <--> j
###     u <--> -u, i <--> j

function addX!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, spropX)
    (; State, X, Par) = Workspace
    N = Par.NumericalParams.N
    (; Npairs, Nsum, siteSum, invpairs) = Par.System

    Vert(n, Rij, s, t, u, flavTransf) =
        V_(State.Gamma, n, s, t, u, flavTransf, Rij, invpairs[Rij], N)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)
    flavTransf12 = (wpw1 * wpw2 > 0, ns * wpw2 > 0, ns * wpw1 < 0)
    flavTransf34 = (wmw3 * wmw4 < 0, ns * wmw4 > 0, ns * wmw3 > 0)

    # get fields of siteSum struct as Matrices for better use of LoopVectorization
    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    X_sum = zeros(42)
    for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij
        fill!(X_sum, 0.0)
        sumsum = 0
        for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort! 
            ki, kj, m, xk =
                S_ki[k_spl, Rij], S_kj[k_spl, Rij], S_m[k_spl, Rij], S_xk[k_spl, Rij]
            Ptm = @SMatrix [
                spropX[1, 1, xk] spropX[1, 2, xk] spropX[1, 3, xk]
                spropX[2, 1, xk] spropX[2, 2, xk] spropX[2, 3, xk]
                spropX[3, 1, xk] spropX[3, 2, xk] spropX[3, 3, xk]
            ]
            Ptm = Ptm * m ### Props now contains two flavor indices

            V12(n) = Vert(n, ki, ns, wpw1, -wpw2, flavTransf12)
            V34(n) = Vert(n, kj, ns, -wmw3, -wmw4, flavTransf34)

            X_sum[fd.yy] +=
                -V12(fd.yy) * V34(fd.yy) * Ptm[2, 2] -
                V12(fd.yz1) * V34(fd.zy1) * Ptm[3, 3] -
                V12(fd.yx1) * V34(fd.xy1) * Ptm[1, 1]
            X_sum[fd.zz] +=
                -V12(fd.zz) * V34(fd.zz) * Ptm[3, 3] -
                V12(fd.zx1) * V34(fd.xz1) * Ptm[1, 1] -
                V12(fd.zy1) * V34(fd.yz1) * Ptm[2, 2]
            X_sum[fd.xx] +=
                -V12(fd.xx) * V34(fd.xx) * Ptm[1, 1] -
                V12(fd.xy1) * V34(fd.yx1) * Ptm[2, 2] -
                V12(fd.xz1) * V34(fd.zx1) * Ptm[3, 3]

            ### Xab1 = -Vaa Vab1 - Vab1 Vbb - Vac1 Vcb1
            X_sum[fd.xy1] +=
                -V12(fd.xx) * V34(fd.xy1) * Ptm[1, 1] -
                V12(fd.xy1) * V34(fd.yy) * Ptm[2, 2] - V12(fd.xz1) * V34(fd.zy1) * Ptm[3, 3]
            X_sum[fd.xz1] +=
                -V12(fd.xx) * V34(fd.xz1) * Ptm[1, 1] -
                V12(fd.xz1) * V34(fd.zz) * Ptm[3, 3] - V12(fd.xy1) * V34(fd.yz1) * Ptm[2, 2]
            X_sum[fd.yx1] +=
                -V12(fd.yy) * V34(fd.yx1) * Ptm[2, 2] -
                V12(fd.yx1) * V34(fd.xx) * Ptm[1, 1] - V12(fd.yz1) * V34(fd.zx1) * Ptm[3, 3]
            X_sum[fd.yz1] +=
                -V12(fd.yy) * V34(fd.yz1) * Ptm[2, 2] -
                V12(fd.yz1) * V34(fd.zz) * Ptm[3, 3] - V12(fd.yx1) * V34(fd.xz1) * Ptm[1, 1]
            X_sum[fd.zx1] +=
                -V12(fd.zz) * V34(fd.zx1) * Ptm[3, 3] -
                V12(fd.zx1) * V34(fd.xx) * Ptm[1, 1] - V12(fd.zy1) * V34(fd.yx1) * Ptm[2, 2]
            X_sum[fd.zy1] +=
                -V12(fd.zz) * V34(fd.zy1) * Ptm[3, 3] -
                V12(fd.zy1) * V34(fd.yy) * Ptm[2, 2] - V12(fd.zx1) * V34(fd.xy1) * Ptm[1, 1]

            ### Xab2 = -Vab2 Vab2 - Vab3 Vba3
            X_sum[fd.xy2] +=
                -V12(fd.xy2) * V34(fd.xy2) * Ptm[1, 2] -
                V12(fd.xy3) * V34(fd.yx3) * Ptm[2, 1]
            X_sum[fd.xz2] +=
                -V12(fd.xz2) * V34(fd.xz2) * Ptm[1, 3] -
                V12(fd.xz3) * V34(fd.zx3) * Ptm[3, 1]
            X_sum[fd.yx2] +=
                -V12(fd.yx2) * V34(fd.yx2) * Ptm[2, 1] -
                V12(fd.yx3) * V34(fd.xy3) * Ptm[1, 2]
            X_sum[fd.yz2] +=
                -V12(fd.yz2) * V34(fd.yz2) * Ptm[2, 3] -
                V12(fd.yz3) * V34(fd.zy3) * Ptm[3, 2]
            X_sum[fd.zx2] +=
                -V12(fd.zx2) * V34(fd.zx2) * Ptm[3, 1] -
                V12(fd.zx3) * V34(fd.xz3) * Ptm[1, 3]
            X_sum[fd.zy2] +=
                -V12(fd.zy2) * V34(fd.zy2) * Ptm[3, 2] -
                V12(fd.zy3) * V34(fd.yz3) * Ptm[2, 3]

            ### Xab3 = -Vab2 Vab3 - Vab3 Vba2
            X_sum[fd.xy3] +=
                -V12(fd.xy2) * V34(fd.xy3) * Ptm[1, 2] -
                V12(fd.xy3) * V34(fd.yx2) * Ptm[2, 1]
            X_sum[fd.xz3] +=
                -V12(fd.xz2) * V34(fd.xz3) * Ptm[1, 3] -
                V12(fd.xz3) * V34(fd.zx2) * Ptm[3, 1]
            X_sum[fd.yx3] +=
                -V12(fd.yx2) * V34(fd.yx3) * Ptm[2, 1] -
                V12(fd.yx3) * V34(fd.xy2) * Ptm[1, 2]
            X_sum[fd.yz3] +=
                -V12(fd.yz2) * V34(fd.yz3) * Ptm[2, 3] -
                V12(fd.yz3) * V34(fd.zy2) * Ptm[3, 2]
            X_sum[fd.zx3] +=
                -V12(fd.zx2) * V34(fd.zx3) * Ptm[3, 1] -
                V12(fd.zx3) * V34(fd.xz2) * Ptm[1, 3]
            X_sum[fd.zy3] +=
                -V12(fd.zy2) * V34(fd.zy3) * Ptm[3, 2] -
                V12(fd.zy3) * V34(fd.yz2) * Ptm[2, 3]
        end

        X[:, Rij, is, it, iu] .+= X_sum
    end
    return
end

function addY!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, spropY)
    (; State, X, Par) = Workspace
    N = Par.NumericalParams.N
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    Vert(n, Rij, s, t, u, flavTransf) =
        V_(State.Gamma, n, s, t, u, flavTransf, Rij, invpairs[Rij], N)
    ns = is - 1
    nt = it - 1
    nu = iu - 1
    wpw1, wpw2, wpw3, wpw4, wmw1, wmw2, wmw3, wmw4 = mixedFrequencies(ns, nt, nu, nwpr)
    flavTransf13 = (nt * wmw3 < 0, wmw1 * wmw3 > 0, wmw1 * nt > 0)
    flavTransf24 = (nt * wpw4 < 0, wpw2 * wpw4 > 0, wpw2 * nt > 0)
    flavTransf31 = (nt * wmw1 > 0, wmw3 * wmw1 > 0, wmw3 * nt < 0)
    flavTransf42 = (nt * wpw2 > 0, wpw4 * wpw2 > 0, wpw4 * nt < 0)

    X_sum = zeros(42)

    # Xtilde only defined for nonlocal pairs Rij != Rii
    for Rij = 1:Npairs
        Rij in OnsitePairs && continue
        # loop over all left hand side inequivalent pairs Rij
        Rji = invpairs[Rij] # store pair corresponding to Rji (easiest case: Rji = Rij) 
        (; xi, xj) = PairTypes[Rij]

        function P_(n::Int, m::Int)
            return spropY[n, m, xi, xj]
        end

        function PT_(n::Int, m::Int)
            return spropY[m, n, xj, xi]
        end

        V13(n) = Vert(n, Rij, -wmw1, nt, wmw3, flavTransf13)
        V24(n) = Vert(n, Rij, wpw2, -nt, -wpw4, flavTransf24)

        V31(n) = Vert(n, Rij, wmw3, nt, -wmw1, flavTransf31)
        V42(n) = Vert(n, Rij, -wpw4, -nt, wpw2, flavTransf42)

        fill!(X_sum, 0.0)

        ### Yaa = Vaa Vaa + Vab2 Vab2 + Vac2 Vac2 + (w -- -w + t)

        X_sum[21+fd.xx] += (
            (
                V13(fd.xx) * V24(fd.xx) * P_(1, 1) +
                V13(fd.xy2) * V24(fd.xy2) * P_(2, 2) +
                V13(fd.xz2) * V24(fd.xz2) * P_(3, 3)
            ) + (
                V31(fd.xx) * V42(fd.xx) * PT_(1, 1) +
                V31(fd.xy2) * V42(fd.xy2) * PT_(2, 2) +
                V31(fd.xz2) * V42(fd.xz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.yy] += (
            (
                V13(fd.yy) * V24(fd.yy) * P_(2, 2) +
                V13(fd.yx2) * V24(fd.yx2) * P_(1, 1) +
                V13(fd.yz2) * V24(fd.yz2) * P_(3, 3)
            ) + (
                V31(fd.yy) * V42(fd.yy) * PT_(2, 2) +
                V31(fd.yx2) * V42(fd.yx2) * PT_(1, 1) +
                V31(fd.yz2) * V42(fd.yz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.zz] += (
            (
                V13(fd.zz) * V24(fd.zz) * P_(3, 3) +
                V13(fd.zx2) * V24(fd.zx2) * P_(1, 1) +
                V13(fd.zy2) * V24(fd.zy2) * P_(2, 2)
            ) + (
                V31(fd.zz) * V42(fd.zz) * PT_(3, 3) +
                V31(fd.zx2) * V42(fd.zx2) * PT_(1, 1) +
                V31(fd.zy2) * V42(fd.zy2) * PT_(2, 2)
            )
        )

        ### Yab1 = Vab3 Vab3 + Vab1 Vab1 + (w -- -w + t)

        X_sum[21+fd.xy1] += (
            (V13(fd.xy3) * V24(fd.xy3) * P_(2, 1) + V13(fd.xy1) * V24(fd.xy1) * P_(1, 2)) + (
                V31(fd.xy3) * V42(fd.xy3) * PT_(2, 1) +
                V31(fd.xy1) * V42(fd.xy1) * PT_(1, 2)
            )
        )

        X_sum[21+fd.xz1] += (
            (V13(fd.xz3) * V24(fd.xz3) * P_(3, 1) + V13(fd.xz1) * V24(fd.xz1) * P_(1, 3)) + (
                V31(fd.xz3) * V42(fd.xz3) * PT_(3, 1) +
                V31(fd.xz1) * V42(fd.xz1) * PT_(1, 3)
            )
        )

        X_sum[21+fd.yx1] += (
            (V13(fd.yx3) * V24(fd.yx3) * P_(1, 2) + V13(fd.yx1) * V24(fd.yx1) * P_(2, 1)) + (
                V31(fd.yx3) * V42(fd.yx3) * PT_(1, 2) +
                V31(fd.yx1) * V42(fd.yx1) * PT_(2, 1)
            )
        )

        X_sum[21+fd.yz1] += (
            (V13(fd.yz3) * V24(fd.yz3) * P_(3, 2) + V13(fd.yz1) * V24(fd.yz1) * P_(2, 3)) + (
                V31(fd.yz3) * V42(fd.yz3) * PT_(3, 2) +
                V31(fd.yz1) * V42(fd.yz1) * PT_(2, 3)
            )
        )

        X_sum[21+fd.zx1] += (
            (V13(fd.zx3) * V24(fd.zx3) * P_(1, 3) + V13(fd.zx1) * V24(fd.zx1) * P_(3, 1)) + (
                V31(fd.zx3) * V42(fd.zx3) * PT_(1, 3) +
                V31(fd.zx1) * V42(fd.zx1) * PT_(3, 1)
            )
        )

        X_sum[21+fd.zy1] += (
            (V13(fd.zy3) * V24(fd.zy3) * P_(2, 3) + V13(fd.zy1) * V24(fd.zy1) * P_(3, 2)) + (
                V31(fd.zy3) * V42(fd.zy3) * PT_(2, 3) +
                V31(fd.zy1) * V42(fd.zy1) * PT_(3, 2)
            )
        )

        ### Yab2 = Vaa Vba2 + Vab2 Vbb + Vac2 Vbc2 + (w -- -w + t)

        X_sum[21+fd.xy2] += (
            (
                V13(fd.xx) * V24(fd.yx2) * P_(1, 1) +
                V13(fd.xy2) * V24(fd.yy) * P_(2, 2) +
                V13(fd.xz2) * V24(fd.yz2) * P_(3, 3)
            ) + (
                V31(fd.xx) * V42(fd.yx2) * PT_(1, 1) +
                V31(fd.xy2) * V42(fd.yy) * PT_(2, 2) +
                V31(fd.xz2) * V42(fd.yz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.xz2] += (
            (
                V13(fd.xx) * V24(fd.zx2) * P_(1, 1) +
                V13(fd.xz2) * V24(fd.zz) * P_(3, 3) +
                V13(fd.xy2) * V24(fd.zy2) * P_(2, 2)
            ) + (
                V31(fd.xx) * V42(fd.zx2) * PT_(1, 1) +
                V31(fd.xz2) * V42(fd.zz) * PT_(3, 3) +
                V31(fd.xy2) * V42(fd.zy2) * PT_(2, 2)
            )
        )

        X_sum[21+fd.yx2] += (
            (
                V13(fd.yy) * V24(fd.xy2) * P_(2, 2) +
                V13(fd.yx2) * V24(fd.xx) * P_(1, 1) +
                V13(fd.yz2) * V24(fd.xz2) * P_(3, 3)
            ) + (
                V31(fd.yy) * V42(fd.xy2) * PT_(2, 2) +
                V31(fd.yx2) * V42(fd.xx) * PT_(1, 1) +
                V31(fd.yz2) * V42(fd.xz2) * PT_(3, 3)
            )
        )

        X_sum[21+fd.yz2] += (
            (
                V13(fd.yy) * V24(fd.zy2) * P_(2, 2) +
                V13(fd.yz2) * V24(fd.zz) * P_(3, 3) +
                V13(fd.yx2) * V24(fd.zx2) * P_(1, 1)
            ) + (
                V31(fd.yy) * V42(fd.zy2) * PT_(2, 2) +
                V31(fd.yz2) * V42(fd.zz) * PT_(3, 3) +
                V31(fd.yx2) * V42(fd.zx2) * PT_(1, 1)
            )
        )

        X_sum[21+fd.zx2] += (
            (
                V13(fd.zz) * V24(fd.xz2) * P_(3, 3) +
                V13(fd.zx2) * V24(fd.xx) * P_(1, 1) +
                V13(fd.zy2) * V24(fd.xy2) * P_(2, 2)
            ) + (
                V31(fd.zz) * V42(fd.xz2) * PT_(3, 3) +
                V31(fd.zx2) * V42(fd.xx) * PT_(1, 1) +
                V31(fd.zy2) * V42(fd.xy2) * PT_(2, 2)
            )
        )

        X_sum[21+fd.zy2] += (
            (
                V13(fd.zz) * V24(fd.yz2) * P_(3, 3) +
                V13(fd.zy2) * V24(fd.yy) * P_(2, 2) +
                V13(fd.zx2) * V24(fd.yx2) * P_(1, 1)
            ) + (
                V31(fd.zz) * V42(fd.yz2) * PT_(3, 3) +
                V31(fd.zy2) * V42(fd.yy) * PT_(2, 2) +
                V31(fd.zx2) * V42(fd.yx2) * PT_(1, 1)
            )
        )

        ### Yab3 = Vab3 Vba1 + Vab1 Vba3 + (w -- -w + t)

        X_sum[21+fd.xy3] += (
            (V13(fd.xy3) * V24(fd.yx1) * P_(2, 1) + V13(fd.xy1) * V24(fd.yx3) * P_(1, 2)) + (
                V31(fd.xy3) * V42(fd.yx1) * PT_(2, 1) +
                V31(fd.xy1) * V42(fd.yx3) * PT_(1, 2)
            )
        )

        X_sum[21+fd.xz3] += (
            (V13(fd.xz3) * V24(fd.zx1) * P_(3, 1) + V13(fd.xz1) * V24(fd.zx3) * P_(1, 3)) + (
                V31(fd.xz3) * V42(fd.zx1) * PT_(3, 1) +
                V31(fd.xz1) * V42(fd.zx3) * PT_(1, 3)
            )
        )

        X_sum[21+fd.yx3] += (
            (V13(fd.yx3) * V24(fd.xy1) * P_(1, 2) + V13(fd.yx1) * V24(fd.xy3) * P_(2, 1)) + (
                V31(fd.yx3) * V42(fd.xy1) * PT_(1, 2) +
                V31(fd.yx1) * V42(fd.xy3) * PT_(2, 1)
            )
        )

        X_sum[21+fd.yz3] += (
            (V13(fd.yz3) * V24(fd.zy1) * P_(3, 2) + V13(fd.yz1) * V24(fd.zy3) * P_(2, 3)) + (
                V31(fd.yz3) * V42(fd.zy1) * PT_(3, 2) +
                V31(fd.yz1) * V42(fd.zy3) * PT_(2, 3)
            )
        )

        X_sum[21+fd.zx3] += (
            (V13(fd.zx3) * V24(fd.xz1) * P_(1, 3) + V13(fd.zx1) * V24(fd.xz3) * P_(3, 1)) + (
                V31(fd.zx3) * V42(fd.xz1) * PT_(1, 3) +
                V31(fd.zx1) * V42(fd.xz3) * PT_(3, 1)
            )
        )

        X_sum[21+fd.zy3] += (
            (V13(fd.zy3) * V24(fd.yz1) * P_(2, 3) + V13(fd.zy1) * V24(fd.yz3) * P_(3, 2)) + (
                V31(fd.zy3) * V42(fd.yz1) * PT_(2, 3) +
                V31(fd.zy1) * V42(fd.yz3) * PT_(3, 2)
            )
        )

        X[:, Rij, is, it, iu] .+= X_sum
    end
end

function getXBubble!(
    Workspace::OneLoopWorkspace,
    FlowParameter::Real,
    ComputeType::Type = Float64,
)
    Lam = FlowParameter
    Par = Workspace.Par
    (; N, lenIntw) = Par.NumericalParams
    (; NUnique) = Par.System
    iSigma = Workspace.State.iSigma
    DiSigma = Workspace.Deriv.iSigma

    for is = 1:N, it = 1:N
        spropX = zeros(3, 3, NUnique)
        spropY = zeros(3, 3, NUnique, NUnique)
        ns = is - 1
        nt = it - 1
        for nw = -lenIntw:lenIntw-1 # Matsubara sum

            set_spropX!(
                spropX,
                NUnique,
                iSigma,
                DiSigma,
                Lam,
                nw,
                nw + ns,
                ComputeType,
                Par.NumericalParams,
            )

            set_spropY!(
                spropY,
                NUnique,
                iSigma,
                DiSigma,
                Lam,
                nw,
                nw - nt,
                ComputeType,
                Par.NumericalParams,
            )
            for iu = 1:N
                nu = iu - 1
                if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                    continue
                end
                addY!(Workspace, is, it, iu, nw, spropY) # add to XTilde-type bubble functions

                ### If no u--t symmetry, then add all the bubbles
                ### If use u--t symmetry, then only add for nu smaller then nt (all other obtained by symmetry)
                # if(!Par.Options.use_symmetry || nu<=nt)

                addX!(Workspace, is, it, iu, nw, spropX)
                # end
            end
        end
    end
end


######################################################################
######### FLOW EQUATIONS ## FLOW EQUATIONS ## FLOW EQUATIONS #########
######################################################################

function get_iS(FlowParam::Real, iSigma::SigmaType, NumParams::NumericalParams)
    Lam = FlowParam
    T = NumParams.T

    @inline iSx(x, nw) = iS_(iSigma.x, x, Lam, nw, T) / 2
    @inline iSy(x, nw) = iS_(iSigma.y, x, Lam, nw, T) / 2
    @inline iSz(x, nw) = iS_(iSigma.z, x, Lam, nw, T) / 2

    return iSx, iSy, iSz

end


function get_iG_i(FlowParam::Real, iSigma_i::AbstractArray, NumParams::NumericalParams)
    Lam = FlowParam
    T = NumParams.T
    @inline iG_i(x, nw) = iG_(iSigma_i, x, Lam, nw, T)
    return iG_i
end

function get_iSKat(iSigma, DiSigma, FlowParam::Real, NumParams::NumericalParams)

    T = NumParams.T
    iSKat_x(x, nw) = iSKat_(iSigma.x, DiSigma.x, x, FlowParam, nw, T)
    iSKat_y(x, nw) = iSKat_(iSigma.y, DiSigma.y, x, FlowParam, nw, T)
    iSKat_z(x, nw) = iSKat_(iSigma.z, DiSigma.z, x, FlowParam, nw, T)

    return iSKat_x, iSKat_y, iSKat_z

end

function get_Theta(FlowParam::Real, _::NumericalParams)
    Lam = FlowParam
    Theta(w) = w^2 / (w^2 + Lam^2)
    return Theta
end

function get_f_int_factor(NumParams::NumericalParams)
    return NumParams.T
end

function get_get_w(NumParams::NumericalParams)
    return nw -> get_w(nw, NumParams.T)
end

function get_Dgamma_factor(NumParams::NumericalParams)
    return NumParams.T
end
function get_chi_factor(NumParams::NumericalParams)
    return NumParams.T
end

function get_props_factor(NumParams::NumericalParams)
    return NumParams.T
end


####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################


function flow_parameter_max_min(NumParams::NumericalParams)

    (; lambda_max, lambda_min) = NumParams
    return lambda_max, lambda_min
end



#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################


export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
