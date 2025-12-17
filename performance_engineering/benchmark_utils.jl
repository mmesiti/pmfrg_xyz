
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
    get_ThreadLocalBuffers,
    split_sitesum

function check_addXY_allocations(T)

    workspace, _ = create_synthetic_workspace_square(N = 10, lattice_size = 5)

    Par = workspace.Par
    (; NUnique, siteSum, Npairs, Nsum) = Par.System
    sitesum_split = split_sitesum(siteSum, 16, Npairs, Nsum)

    iuh_blocksize = 10 ÷ 2

    buffs = get_ThreadLocalBuffers(10, Par.System, iuh_blocksize, T, 1)[1]

    # Convert Gamma to ComputeType if needed
    Gamma =
        eltype(workspace.State.Gamma) == T ? workspace.State.Gamma :
        T.(workspace.State.Gamma)


    System = Par.System
    N = Par.NumericalParams.N


    Xh = zeros(T, iuh_blocksize, 42, Npairs, N, N)

    addX!(
        buffs.X_sum_addX,
        Gamma,
        System,
        N,
        1,
        1,
        1,
        buffs.spropX,
        buffs,
        sitesum_split,
        1,
        iuh_blocksize,
    )

    addY!(
        buffs.X_sum_addY,
        Gamma,
        System,
        N,
        1,
        1,
        1,
        buffs.spropY,
        buffs,
        1,
        iuh_blocksize,
    )


    addXallocations = @allocations addX!(
        buffs.X_sum_addX,
        Gamma,
        System,
        N,
        1,
        1,
        1,
        buffs.spropX,
        buffs,
        sitesum_split,
        1,
        iuh_blocksize,
    )

    addYallocations = @allocations addY!(
        buffs.X_sum_addY,
        Gamma,
        System,
        N,
        1,
        1,
        1,
        buffs.spropY,
        buffs,
        1,
        iuh_blocksize,
    )

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
