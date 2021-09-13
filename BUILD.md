# Installing NimbleSM

## Configuring CMake

A very basic installation is done with:
````
cmake ${srcdir}
````
which builds and installs a default NimbleSM when you run `make`.

The full keyword listing is below.

## NimbleSM CMake Option Listing

* ArborX_DIR or AborX_ROOT: STRING
  * Location of ArborX install root.
  * Default: None

* bvh_DIR or bvh_ROOT: STRING
  * Location of BVH install root.
  * Default: None

* ENABLE_DEBUGGING: BOOL
  * Whether to turn on enhanced debugging info.
  * Default: OFF

* HAVE_ARBORX: BOOL
  * Whether to configure with [ArborX library](https://github.com/arborx/ArborX).
  * Default: OFF

* HAVE_BVH: BOOL
  * Whether to use BVH for parallel asynchronous contact.
  * Default: OFF

* HAVE_KOKKOS: BOOL
  * Whether to use [Kokkos library](https://github.com/kokkos/kokkos).
  * Default: OFF

* HAVE_TRILINOS: BOOL
  * Whether to use [Trilinos library](https://github.com/trilinos/Trilinos).
  * Default: OFF

* HAVE_UQ: BOOL
  * Whether to enable UQ sampling.
  * Default: OFF

* Kokkos_DIR or Kokkos_ROOT: STRING
  * Location of Kokkos install root.
  * Default: None

* NIMBLE_ENABLE_UNIT_TESTS
  * Whether to build Nimble with GTest testing.
  * Default: OFF

* TIME_CONTACT
  * Whether to add extra contact timers.
  * Default: OFF

* Trilinos_DIR or Trilinos_ROOT: STRING
  * Location of Trilinos install root.
  * Default: None

* USE_PURE_MPI: BOOL
  * Whether to configure with MPI library.
  * Default: OFF

* USE_TRILINOS_EXODUS: BOOL
  * Whether to use Exodus from Trilinos.
  * Default: OFF


## Examples of scripts or configuration steps

Examples of steps for building the code can be found in 
[the workflows directory](https://github.com/NimbleSM/NimbleSM/tree/develop/.github/workflows) 
(which builds NimbleSM for testing under different configurations)
or in [the scripts directory](https://github.com/NimbleSM/NimbleSM/tree/develop/scripts/cmake-example-scripts).


##### [LICENSE](https://github.com/NimbleSM/NimbleSM/blob/develop/LICENSE)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.
