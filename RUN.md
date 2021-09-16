# Running NimbleSM

## Creating an input file
 
Examples of input files with the extension `.in` can be found 
in [the test directory](https://github.com/NimbleSM/NimbleSM/tree/develop/test).
The four subfolders, namely `contact`, `dynamics`, `quasistatic`, and `uq`, contain input files for four
different types of simulations.

For instance, consider the input file [wave_in_bar.in](https://github.com/NimbleSM/NimbleSM/blob/develop/test/dynamics/wave_in_bar/wave_in_bar.in)
```
genesis input file:               wave_in_bar.g
exodus output file:               wave_in_bar.e
final time:                       1.0e-5
number of load steps:             1000
output frequency:                 500
output fields:                    displacement velocity deformation_gradient ipt01_deformation_gradient ipt02_deformation_gradient ipt03_deformation_gradient ipt04_deformation_gradient ipt05_deformation_gradient ipt06_deformation_gradient ipt07_deformation_gradient ipt08_deformation_gradient stress ipt01_stress ipt02_stress ipt03_stress ipt04_stress ipt05_stress ipt06_stress ipt07_stress ipt08_stress
material parameters:              material_1 neohookean density 7.8 shear_modulus 1.5e12 bulk_modulus 1.0e12
element block:                 block_1 material_1
boundary condition:               initial_velocity nodelist_1 x 1000.0
boundary condition:               prescribed_velocity nodelist_2 x 0.0
boundary condition:               prescribed_velocity nodelist_2 y 0.0
boundary condition:               prescribed_velocity nodelist_2 z 0.0
```
* The `genesis input file` specifies the mesh. When running in parallel with `mpi`, a series of decomposed files with the same prefix will be searched by NimbleSM in the current directory.
* The `exodus output file` specifies the prefix for the Exodus output files.
* By default, this simulation will use the explicit Verlet integration scheme.

## Command-line options

A very basic command to run the code is
````
mpirun -np (# of ranks) $path_to_executable_NimbleSM [FLAG] $input_file [OPTION]
````
which runs with `mpirun` the executable, `NimbleSM`.

### FLAG list

* `--use_kokkos`
   * This flag activates a simulation that employs Kokkos objects. In particular, the `OpenMP` environment can be exploited.
   * Requirement: `NimbleSM` has to be linked to the [Kokkos library](https://github.com/Kokkos/kokkos).
* `--use_tpetra`
   * This flag activates a simulation that employs Tpetra objects. The simulation can use the explicit time integrator
     or a quasistatic time integration.
   * Requirement: `NimbleSM` has to be linked to the [Trilinos library](https://github.com/Trilinos/trilinos).
* `--use_uq`
   * This flag activates a simulation with the Uncertainty Quantification (UQ) module.
* `--use_vt`
   * This flag is only valid when building `NimbleSM` with the BVH library.

### OPTION list

[OPTION] is only accessible when compiling with the libraries BVH and [VT](https://github.com/DARMA-tasking/vt).
This feature will soon be released publicly.


##### [LICENSE](https://github.com/NimbleSM/NimbleSM/blob/develop/LICENSE)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.

