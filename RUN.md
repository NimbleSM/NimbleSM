# Running NimbleSM

## Creating an input file
 
Examples of input files can be found in [the test directory](https://github.com/NimbleSM/NimbleSM/tree/develop/test).
The four subfolders, namely `contact`, `dynamics`, `quasistatic`, and `uq`, contain input files for four
different types of simulations.

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

