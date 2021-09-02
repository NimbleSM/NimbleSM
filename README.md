# NimbleSM
NimbleSM is a Lagrangian finite-element code for solid mechanics.  Its primary application is the solution of mechanics problems on nonuniform, three-dimensional meshes using either explicit transient dynamic or implicit quasi-static time integration.  Additional features include multiscale modeling and the simulation of contact.  The NimbleSM code base is designed for performance portability across varying hardware architectures.

NimbleSM may be built with the following optional dependencies:

## ArborX

Performance portable algorithms for geometric search.
Required for contact simulations.

https://github.com/arborx/ArborX

## Exodus

Library for I/O.  It is highly recommended that you build NimbleSM with exodus support.  The alternative is an exceedingly slow text-based I/O system.

https://github.com/gsjaardema/seacas

## Kokkos

Performance portability programming ecosystem.  Utilized by NimbleSM for improved performance on GPU-based systems and other emerging hardware platforms.
Required for contact simulation with ArborX.

https://github.com/kokkos/kokkos

## MPI

Message passing library

https://www.open-mpi.org

## Trilinos

Algorithms for the solution of large-scale engineering and scientific problems.  Utilized by NimbleSM for parallel linear solvers.  
Required for parallel execution of implicit time integration (quasi-statics).

https://trilinos.org

## BVH & VT

This combination will be publicly released in the near future. 
It combines a novel geometric search library, BVH, with the [VT runtime tasking library](https://github.com/DARMA-tasking/vt) for C++, 
built for general scientific applications.


