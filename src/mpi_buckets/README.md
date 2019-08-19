# Outline of algorithm and functionality of MPI Buckets

## Goal

Each MPI rank contains a section of the simulation. 
We need to determine if any sections of the object are going to collide with one another.
In order to do that, we need each rank to obtain a list of the corresponding ranks that
Are likely to have pieces of their own objects nearby. 

## Methodology

Each MPI rank takes it's partition and identifies a list of grid cells 
that it's partition exists in. It provides those grid cells to a 
distributed hash table, which informs it of the ranks it has to communicate with
in order to resolve collisions. 

