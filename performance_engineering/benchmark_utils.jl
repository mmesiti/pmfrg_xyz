
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
    get_ThreadLocalBuffers

using InteractiveUtils
function check_addXY_allocations(T)

    Workspace, _ = create_synthetic_workspace_square(N = 16, lattice_size = 5)

    Par = Workspace.Par
    N = Par.NumericalParams.N

    # Convert Gamma to ComputeType if needed
    Gamma =
        eltype(Workspace.State.Gamma) == T ? Workspace.State.Gamma :
        T.(Workspace.State.Gamma)


    iuh_blocksize = (T == Float64) ? 4 : 8

    Buffs = first(get_ThreadLocalBuffers(N, Par.System, iuh_blocksize, T, 1))


    is, it, nw = 1, 1, 1
    iuh_start = 1

    addX!(
        Buffs.X_sum_addX,
        Buffs.V12_addX,
        Buffs.V34_addX,
        Gamma,
        Workspace.Par.System,
        N,
        is,
        it,
        nw,
        Buffs.spropX,
        iuh_start,
    )

    addY!(
        Buffs.X_sum_addY,
        Gamma,
        Workspace.Par.System,
        N,
        is,
        it,
        nw,
        Buffs.spropY,
        Buffs,
        iuh_start,
    )


    addXallocations = @allocations addX!(
        Buffs.X_sum_addX,
        Buffs.V12_addX,
        Buffs.V34_addX,
        Gamma,
        Workspace.Par.System,
        N,
        is,
        it,
        nw,
        Buffs.spropX,
        iuh_start,
    )

    addYallocations = @allocations addY!(
        Buffs.X_sum_addY,
        Gamma,
        Workspace.Par.System,
        N,
        is,
        it,
        nw,
        Buffs.spropY,
        Buffs,
        iuh_start,
    )

    if addXallocations > 2
        open("addX-warntype.txt", "w") do f
            code_warntype(
                f,
                addX!,
                typeof.((
                    Buffs.X_sum_addX,
                    Buffs.V12_addX,
                    Buffs.V34_addX,
                    Gamma,
                    Workspace.Par.System,
                    N,
                    is,
                    it,
                    nw,
                    Buffs.spropX,
                    iuh_start,
                )),
            )
        end
    end




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
