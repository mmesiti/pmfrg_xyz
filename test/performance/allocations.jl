import PMFRG_xyz: addX!, addY!, ThreadLocalBuffersT
include("../../performance_engineering/benchmark_utils.jl")

function check_addXY_allocations()

    workspace, _ = create_synthetic_workspace_square(N=10, lattice_size=5)

    Par = workspace.Par
    (; NUnique, Npairs) = Par.System


    buffs = ThreadLocalBuffersT( zeros((21, Npairs)),
              zeros((21, Npairs)),
              zeros(21),
              zeros(3, 3, NUnique),
              zeros(3, 3, NUnique, NUnique),
              zeros(3,3),
              zeros(21),
              zeros(21),
              zeros(21),
              zeros(21))




    addX!(workspace,1,1,2,1,buffs.spropX,buffs)
    addY!(workspace,1,1,2,1,buffs.spropY,buffs)

    addXallocations = @allocations addX!(workspace,1,1,2,1,buffs.spropX,buffs)
    @test addXallocations <= 1

    addYallocations = @allocations addY!(workspace,1,1,2,1,buffs.spropY,buffs)
    @test addYallocations <= 1 
end


