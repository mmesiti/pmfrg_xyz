"""
Example configurations for profiling and benchmarking PMFRG_xyz.

Provides setup functions for different lattice examples (dimer, square lattice).
"""

module ExampleSetups
using SpinFRGLattices
using SpinFRGLattices.SquareLattice
import PMFRG_xyz: Params

"""
    setup_dimer_example(; N=8)

Set up the dimer anisotropy example system and parameters.

Returns: (system, params, isotropy_matrix)
"""
function dimer(; N = 8)
    system = SpinFRGLattices.getPolymer(2)
    par = Params(system, N = N, accuracy = 1e-10, temp_max = 10.0, temp_min = 1.0)

    isotropy = zeros(system.Npairs, 3)
    for n = 1:system.Npairs
        isotropy[n, :] = [1.0, 0.5, 0.2]
    end

    return par, isotropy
end

"""
    square_lattice(; lattice_size=6, N=8)

Set up the square lattice anisotropy example system and parameters.

Returns: (system, params, isotropy_matrix)
"""
function square_lattice(; lattice_size = 6, N = 8)
    J1 = 1.0
    J2 = 0.5

    system = getSquareLattice(lattice_size, [J1, J2])
    par = Params(system, N = N, accuracy = 1e-10, temp_max = 10.0, temp_min = 1.0)

    isotropy = zeros(system.Npairs, 3)
    for n = 1:system.Npairs
        isotropy[n, :] = [1.0, 0.5, 0.2]
    end

    return par, isotropy
end

"""
    square_lattice_large(; lattice_size=10, N=20)

Set up the square lattice anisotropy example system and parameters.

Returns: (system, params, isotropy_matrix)
"""
function square_lattice_large(; lattice_size = 10, N = 20)
    J1 = 1.0
    J2 = 0.5

    system = getSquareLattice(lattice_size, [J1, J2])
    par = Params(system, N = N, accuracy = 1e-10, temp_max = 10.0, temp_min = 1.0)

    isotropy = zeros(system.Npairs, 3)
    for n = 1:system.Npairs
        isotropy[n, :] = [1.0, 0.5, 0.2]
    end

    return par, isotropy
end


example_setups = Dict("dimer" => dimer, 
                      "square_lattice" => square_lattice,
                      "square_lattice_large" => square_lattice_large)

end


