
using SpinFRGLattices
using SpinFRGLattices.SquareLattice
import PMFRG_xyz:
    Params,
    AllocateSetup,
    InitializeState,
    OneLoopWorkspace,
    addX!,
    addY!,
    Xh_from_X,
    ThreadLocalBuffersT

function check_addXY_allocations()

    workspace, _ = create_synthetic_workspace_square(N = 10, lattice_size = 5)

    Par = workspace.Par
    (; NUnique) = Par.System


    buffs = ThreadLocalBuffersT(
        zeros(21),
        zeros(3, 3, NUnique),
        zeros(3, 3, NUnique, NUnique),
        zeros(3, 3),
        zeros(21),
        zeros(21),
        zeros(21),
        zeros(21),
    )


    X = workspace.X
    Gamma = workspace.State.Gamma
    System = Par.System
    N = Par.NumericalParams.N


    Xh = Xh_from_X(X)

    addX!(Xh, Gamma, System, N, 1, 1, 1, buffs.spropX, buffs)
    addY!(Xh, Gamma, System, N, 1, 1, 1, buffs.spropY, buffs)

    addXallocations = @allocations addX!(Xh, Gamma, System, N, 1, 1, 1, buffs.spropX, buffs)

    addYallocations = @allocations addY!(Xh, Gamma, System, N, 1, 1, 1, buffs.spropY, buffs)

    addXallocations, addYallocations
end


function create_synthetic_workspace_dimer(; N::Int = 8)
    System = SpinFRGLattices.getPolymer(2)
    par = Params(System, N = N, temp_max = 10.0, temp_min = 1.0)
    isotropy = zeros(System.Npairs, 3)

    for n = 1:System.Npairs
        isotropy[n, :] = [1.0, 0.5, 0.2]
    end

    State = InitializeState(par, isotropy)
    setup = AllocateSetup(par)

    X = setup[1]
    Par = setup[end]
    Deriv = copy(State)
    fill_with_zeros!(Deriv)

    workspace = OneLoopWorkspace(State, Deriv, X, Par)
    lam = 5.0

    return workspace, lam
end

function create_synthetic_workspace_square(; N::Int = 8, lattice_size::Int = 4)
    J1 = 1.0
    J2 = 0.5

    System = getSquareLattice(lattice_size, [J1, J2])
    par = Params(System, N = N, temp_max = 10.0, temp_min = 1.0)
    isotropy = zeros(System.Npairs, 3)

    for n = 1:System.Npairs
        isotropy[n, :] = [1.0, 0.5, 0.2]
    end

    State = InitializeState(par, isotropy)
    setup = AllocateSetup(par)

    (; X, Par) = setup
    Deriv = copy(State)
    fill_with_zeros!(Deriv)

    workspace = OneLoopWorkspace(State, Deriv, X, Par)
    lam = 5.0

    return workspace, lam
end

# level 3
function fill_with_zeros!(state)
    if hasfield(typeof(state), :f_int)
        fill!(state.f_int, 0.0)
    end
    if hasfield(typeof(state), :iSigma)
        if hasfield(typeof(state.iSigma), :x)
            fill!(state.iSigma.x, 0.0)
        end
        if hasfield(typeof(state.iSigma), :y)
            fill!(state.iSigma.y, 0.0)
        end
        if hasfield(typeof(state.iSigma), :z)
            fill!(state.iSigma.z, 0.0)
        end
    end
    if hasfield(typeof(state), :Gamma)
        fill!(state.Gamma, 0.0)
    end
end
