# Examples 

This directory contains some examples 
that show how this code can be used.

Precise and up to date information 
can also be found in the scripts contained in 
`../github/workflows/test.yml`
(the continuous integration definition file)
which are run regularly.

## Julia Environment for examples
This directory contains also a `Project.toml` and a `Manifest.toml`
that define an environment 
that can be used to run these examples.
Refer to the documentation of [Pkg](https://docs.julialang.org/en/v1/stdlib/Pkg/) 
(the Julia Package manager)
for futther information.

## Example scripts

These minimal examples launch the `SolveFRG` function,
time it,
and save the results using the JLD2 library:

- `dimer-anisotropy.jl`
- `square-lattice-anisotropy.jl`

The examples here can be run from the terminal using the command

```bash
cd examples
julia --project=. dimer-anisotropy.jl  # shorter example
julia --project=. square-lattice-anisotropy.jl  # more complex, but still quick
```

Notice that `--project=.` activates the 
[environment](#julia-environment-for-examples) 
 in the current directory.

## Slurm example

The file `slurm-example-paderborn.sh`
can be used as an inspiration to run the simulation 
on a Slurm-managed cluster.

> [!CAUTION]
> Verify that all the `#SBATCH` directives in your slurm script 
> make sense and are appropriate for the problem you are launching
> and the cluster/partition you are running on.  
> Refer to the documentation of the HPC system you are using 
> to determine the right value of all the options.

Usage:

``` bash
cd examples 
sbatch <sbatch options> ./slurm-example-paderborn.sh dimer-anisotropy.jl 
# Or
sbatch <sbatch options> ./slurm-example-paderborn.sh square-lattice-anisotropy.jl
```

> [!IMPORTANT]
> There are many ways to set up Julia. 
> You might have to amend the commands and scripts shown here approprately,
> to be compatible with the way Julia is set up 
> on the machine/cluster you intend to launch the code.

