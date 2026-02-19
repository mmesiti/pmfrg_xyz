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
using LinearAlgebra
using SparseArrays

function addX!(Workspace, is::Integer, it::Integer, iu::Integer, nwpr::Integer, Props)
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
                Props[xk, xk, 1, 1] Props[xk, xk, 1, 2] Props[xk, xk, 1, 3]
                Props[xk, xk, 2, 1] Props[xk, xk, 2, 2] Props[xk, xk, 2, 3]
                Props[xk, xk, 3, 1] Props[xk, xk, 3, 2] Props[xk, xk, 3, 3]
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

function addY!(
    Workspace,
    is::Integer,
    it::Integer,
    iu::Integer,
    nwpr::Integer,
    Props;
    _l = 1.0,
)
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
            return Props[xi, xj, n, m]
        end

        function PT_(n::Int, m::Int)
            return Props[xj, xi, m, n]
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

function getKataninProp!(BubbleProp, NUnique, iSigma, DiSigma, T, Lam, nw1, nw2)
    iGx(x, nw) = iG_(iSigma.x, x, Lam, nw, T)
    iGy(x, nw) = iG_(iSigma.y, x, Lam, nw, T)
    iGz(x, nw) = iG_(iSigma.z, x, Lam, nw, T)

    iSKatx(x, nw) = iSKat_(iSigma.x, DiSigma.x, x, Lam, nw, T)
    iSKaty(x, nw) = iSKat_(iSigma.y, DiSigma.y, x, Lam, nw, T)
    iSKatz(x, nw) = iSKat_(iSigma.z, DiSigma.z, x, Lam, nw, T)

    for i = 1:NUnique, j = 1:NUnique
        BubbleProp[i, j, 1, 1] = iSKatx(i, nw1) * iGx(j, nw2) * T
        BubbleProp[i, j, 1, 2] = iSKatx(i, nw1) * iGy(j, nw2) * T
        BubbleProp[i, j, 1, 3] = iSKatx(i, nw1) * iGz(j, nw2) * T
        BubbleProp[i, j, 2, 1] = iSKaty(i, nw1) * iGx(j, nw2) * T
        BubbleProp[i, j, 2, 2] = iSKaty(i, nw1) * iGy(j, nw2) * T
        BubbleProp[i, j, 2, 3] = iSKaty(i, nw1) * iGz(j, nw2) * T
        BubbleProp[i, j, 3, 1] = iSKatz(i, nw1) * iGx(j, nw2) * T
        BubbleProp[i, j, 3, 2] = iSKatz(i, nw1) * iGy(j, nw2) * T
        BubbleProp[i, j, 3, 3] = iSKatz(i, nw1) * iGz(j, nw2) * T
    end

    ### Relative minus sign between paper & Nils' thesis
    return -BubbleProp
    # return SMatrix{NUnique, NUnique, 3, 3}(BubbleProp)
    ### SMatrix can only create 2d array (according to ChatGPT). Use SArray instead
end


function getXBubble!(Workspace, FlowParameter)
    Lam = FlowParameter
    Par = Workspace.Par
    (; T, N, lenIntw) = Par.NumericalParams
    (; NUnique) = Par.System
    iSigma = Workspace.State.iSigma
    DiSigma = Workspace.Deriv.iSigma

    for is = 1:N, it = 1:N
        BubbleProp = zeros(NUnique, NUnique, 3, 3)
        ns = is - 1
        nt = it - 1
        for nw = -lenIntw:lenIntw-1 # Matsubara sum

            spropX =
                getKataninProp!(BubbleProp, NUnique, iSigma, DiSigma, T, Lam, nw, nw + ns)
            spropY =
                getKataninProp!(BubbleProp, NUnique, iSigma, DiSigma, T, Lam, nw, nw - nt)
            for iu = 1:N
                nu = iu - 1
                if (ns + nt + nu) % 2 == 0# skip unphysical bosonic frequency combinations
                    continue
                end
                addY!(Workspace, is, it, iu, nw, spropY, _l = Lam) # add to XTilde-type bubble functions

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

function getDFint!(Workspace, FlowParam::Real)
    Lam = FlowParam
    (; State, Deriv, Par) = Workspace
    (; T, lenIntw_acc) = Par.NumericalParams
    NUnique = Par.System.NUnique

    iSigmax(x, nw) = iSigma_(State.iSigma.x, x, nw)
    iSigmay(x, nw) = iSigma_(State.iSigma.y, x, nw)
    iSigmaz(x, nw) = iSigma_(State.iSigma.z, x, nw)

    iGx(x, nw) = iG_(State.iSigma.x, x, Lam, nw, T)
    iGy(x, nw) = iG_(State.iSigma.y, x, Lam, nw, T)
    iGz(x, nw) = iG_(State.iSigma.z, x, Lam, nw, T)

    iSx(x, nw) = iS_(State.iSigma.x, x, Lam, nw, T)
    iSy(x, nw) = iS_(State.iSigma.y, x, Lam, nw, T)
    iSz(x, nw) = iS_(State.iSigma.z, x, Lam, nw, T)

    Theta(Lam, w) = w^2 / (w^2 + Lam^2)

    for x = 1:NUnique
        sumres = 0.0
        for nw = -lenIntw_acc:lenIntw_acc-1
            w = get_w(nw, T)
            sumres += iSx(x, nw) / iGy(x, nw) * Theta(Lam, w) * iSigmax(x, nw) / w
            sumres += iSy(x, nw) / iGy(x, nw) * Theta(Lam, w) * iSigmay(x, nw) / w
            sumres += iSz(x, nw) / iGz(x, nw) * Theta(Lam, w) * iSigmaz(x, nw) / w
        end
        Deriv.f_int[x] = -0.5 * T * sumres
    end
end

function get_Self_Energy!(Workspace, FlowParam)
    Lam = FlowParam
    Par = Workspace.Par
    @inline iSx(x, nw) =
        iS_(Workspace.State.iSigma.x, x, Lam, nw, Par.NumericalParams.T) / 2
    @inline iSy(x, nw) =
        iS_(Workspace.State.iSigma.y, x, Lam, nw, Par.NumericalParams.T) / 2
    @inline iSz(x, nw) =
        iS_(Workspace.State.iSigma.z, x, Lam, nw, Par.NumericalParams.T) / 2
    compute1PartBubble!(Workspace.Deriv.iSigma, Workspace.State.Gamma, [iSx, iSy, iSz], Par)
end

function addTo1PartBubble!(Dgamma::SigmaType, Gamma_::Function, Props, Par)

    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; siteSum, Nsum, OnsitePairs) = Par.System

    Threads.@threads for iw1 = 1:N
        nw1 = iw1 - 1
        for (x, Rx) in enumerate(OnsitePairs)
            for nw = -lenIntw_acc:lenIntw_acc-1
                jsum = zeros(3)
                wpw1 = nw1 + nw + 1
                wmw1 = nw - nw1
                flavTransform = (wmw1 * wpw1 < 0, false, false)
                for k_spl = 1:Nsum[Rx]
                    (; m, ki, xk) = siteSum[k_spl, Rx]
                    gam = @SVector [
                        Gamma_(n, ki, 0, -wmw1, -wpw1, flavTransform) for n = 1:21
                    ]

                    jsum[fd.xx] +=
                        (
                            gam[fd.xx] * Props[1](xk, nw) +
                            gam[fd.yx1] * Props[2](xk, nw) +
                            gam[fd.zx1] * Props[3](xk, nw)
                        ) * m
                    jsum[fd.yy] +=
                        (
                            gam[fd.xy1] * Props[1](xk, nw) +
                            gam[fd.yy] * Props[2](xk, nw) +
                            gam[fd.zy1] * Props[3](xk, nw)
                        ) * m
                    jsum[fd.zz] +=
                        (
                            gam[fd.xz1] * Props[1](xk, nw) +
                            gam[fd.yz1] * Props[2](xk, nw) +
                            gam[fd.zz] * Props[3](xk, nw)
                        ) * m
                end
                Dgamma.x[x, iw1] += -T * jsum[1]
                Dgamma.y[x, iw1] += -T * jsum[2]
                Dgamma.z[x, iw1] += -T * jsum[3]
            end
        end
    end
end

struct OneLoopWorkspace{T,ParType}
    State::StateType{T}
    Deriv::StateType{T}
    X::Array{T,5}
    Par::ParType
end

function OneLoopWorkspace(State, Deriv, X, Par)
    setZero!(Deriv)
    setZero!(X)

    return OneLoopWorkspace(StateType(State.x...), StateType(Deriv.x...), X, Par)
end

using JLD2
function getDeriv!(Deriv, State, setup, FlowParameter; saveArgs = true)

    (; X, Par) = setup # use pre-allocated X and XTilde to reduce garbage collector time

    Workspace = OneLoopWorkspace(State, Deriv, X, Par)

    getDFint!(Workspace, FlowParameter)
    get_Self_Energy!(Workspace, FlowParameter)
    getXBubble!(Workspace, FlowParameter)

    symmetrizeBubble!(Workspace.X, Par)
    addToVertexFromBubble!(Workspace.Deriv.Gamma, Workspace.X)
    symmetrizeVertex!(Workspace.Deriv.Gamma, Par)

    return
end

####################################################
######### SOLVE ## SOLVE ## SOLVE ## SOLVE #########
####################################################



function launchPMFRG!(State, setup, Deriv!::Function; method = DP5(), npoints = 600)

    Par = setup[end]
    (; lambda_max, lambda_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(lambda_max)
    tend = get_t_min(lambda_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    saved_values = SavedValues(eltype(State), Observables{eltype(State)})

    function save_func(State, t, integrator)
        chi_x = getChi_x(State, t_to_Lam(t), Par)
        chi_y = getChi_y(State, t_to_Lam(t), Par)
        chi_z = getChi_z(State, t_to_Lam(t), Par)

        return Observables(copy(chi_x), copy(chi_y), copy(chi_z))
    end

    ObsSaveat = gettMesh(lambda_min, lambda_max, npoints)
    saveCB = SavingCallback(
        save_func,
        saved_values,
        save_everystep = false,
        saveat = ObsSaveat,
        tdir = -1,
    )

    problem = ODEProblem(Deriv_subst!, State, (t0, tend), setup) # function, initial state, timespan, ??
    sol = solve(
        problem,
        method,
        reltol = accuracy,
        abstol = accuracy,
        save_everystep = false,
        callback = saveCB,
        dt = Lam_to_t(0.2 * lambda_max),
    )

    return sol, saved_values
end

function testPMFRG!(State, setup, Deriv!::Function; loadArgs = false)
    Par = setup[end]
    (; lambda_max, lambda_min, accuracy) = Par.NumericalParams

    t0 = Lam_to_t(lambda_max)
    tend = get_t_min(lambda_min)
    Deriv_subst! = generateSubstituteDeriv(Deriv!)

    der = copy(State)
    setZero!(der)

    Deriv_subst!(der, State, setup, t0, s = false)
end

SolveFRG(Par, isotropy; kwargs...) =
    launchPMFRG!(InitializeState(Par, isotropy), AllocateSetup(Par), getDeriv!; kwargs...)
TestFRG(Par, isotropy; kwargs...) =
    testPMFRG!(InitializeState(Par, isotropy), AllocateSetup(Par), getDeriv!; kwargs...)



#############################################################
######### OBSERVABLES ## OBSERVABLES ## OBSERVABLES #########
#############################################################

getChi_z(State::ArrayPartition, Lam::Real, Par) =
    getChi_z(State.x[2], State.x[3], State.x[5], Lam, Par)
getChi_x(State::ArrayPartition, Lam::Real, Par) =
    getChi_x(State.x[3], State.x[4], State.x[5], Lam, Par)
getChi_y(State::ArrayPartition, Lam::Real, Par) =
    getChi_y(State.x[4], State.x[2], State.x[5], Lam, Par)

function getChi_z(
    iSigmaX::AbstractArray,
    iSigmaY::AbstractArray,
    Gamma::AbstractArray,
    Lam::Real,
    Par,
)
    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGx(x, w) = iG_(iSigmaX, x, Lam, w, T)
    iGy(x, w) = iG_(iSigmaY, x, Lam, w, T)
    Vxy2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, fd.xy2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += T * iGx(xi, nK) * iGy(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGx(xi, nK)^2 * iGy(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += T^2 * GGGG * Vxy2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end

function getChi_x(
    iSigmaY::AbstractArray,
    iSigmaZ::AbstractArray,
    Gamma::AbstractArray,
    Lam::Real,
    Par,
)
    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGy(x, w) = iG_(iSigmaY, x, Lam, w, T)
    iGz(x, w) = iG_(iSigmaZ, x, Lam, w, T)
    Vyz2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, fd.yz2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += T * iGy(xi, nK) * iGz(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGy(xi, nK)^2 * iGz(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += T^2 * GGGG * Vyz2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end

function getChi_y(
    iSigmaZ::AbstractArray,
    iSigmaX::AbstractArray,
    Gamma::AbstractArray,
    Lam::Real,
    Par,
)
    (; T, N, lenIntw_acc) = Par.NumericalParams
    (; Npairs, invpairs, PairTypes, OnsitePairs) = Par.System

    iGz(x, w) = iG_(iSigmaZ, x, Lam, w, T)
    iGx(x, w) = iG_(iSigmaX, x, Lam, w, T)
    Vzx2(Rij, s, t, u, isFlavorTransform) =
        V_(Gamma, fd.zx2, s, t, u, isFlavorTransform, Rij, invpairs[Rij], N)

    Chi = zeros(_getFloatType(Par), Npairs)

    for Rij = 1:Npairs
        (; xi, xj) = PairTypes[Rij]
        for nK = -lenIntw_acc:lenIntw_acc-1
            if Rij in OnsitePairs
                Chi[Rij, 1] += T * iGz(xi, nK) * iGx(xi, nK)
            end
            for nK2 = -lenIntw_acc:lenIntw_acc-1
                npwpw2 = nK + nK2 + 1
                w2mw = nK2 - nK
                #use that Vc_0 is calculated from Vb
                GGGG = iGz(xi, nK)^2 * iGx(xj, nK2)^2
                flavTransform = (npwpw2 * w2mw > 0, false, false)
                Chi[Rij] += T^2 * GGGG * Vzx2(Rij, 0, npwpw2, -w2mw, flavTransform)
            end
        end
    end
    return (Chi)
end

export Params, SolveFRG, TestFRG, getChi_x, getChi_y, getChi_z

end
