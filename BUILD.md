# Installing NimbleSM

⚠️ WARNING:

The build instructions below are outdated and no longer valid for the latest versions of NimbleSM and its dependencies. A new, updated procedure has been tested and validated on Ubuntu 24.04. Please refer to the new instructions provided at the end of this document for a working example [here](#examples-of-scripts-or-configuration-steps)

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

* Kokkos_DIR or Kokkos_ROOT: STRING
  * Location of Kokkos install root.
  * Default: None

* NimbleSM_ENABLE_ARBORX: BOOL
  * Whether to configure with [ArborX library](https://github.com/arborx/ArborX).
  * Default: OFF

* NimbleSM_ENABLE_BVH: BOOL
  * Whether to use BVH for parallel asynchronous contact.
  * Default: OFF

* NimbleSM_ENABLE_DEBUG: BOOL
  * Whether to turn on enhanced debugging info.
  * Default: OFF

* NimbleSM_ENABLE_KOKKOS: BOOL
  * Whether to use [Kokkos library](https://github.com/kokkos/kokkos).
  * Default: OFF

* NimbleSM_ENABLE_MPI: BOOL
  * Whether to configure with MPI library.
  * Default: OFF

* NimbleSM_ENABLE_TRILINOS: BOOL
  * Whether to use [Trilinos library](https://github.com/trilinos/Trilinos).
  * Default: OFF

* NimbleSM_ENABLE_UNIT_TESTS
  * Whether to build Nimble with unit testing.
  * Default: OFF

* NimbleSM_ENABLE_UQ: BOOL
  * Whether to enable UQ sampling.
  * Default: OFF

* NimbleSM_TIME_Contact
  * Whether to add extra contact timers.
  * Default: OFF

* NimbleSM_USE_TRILINOS_EXODUS: BOOL
  * Whether to use Exodus from Trilinos.
  * Default: OFF

* Trilinos_DIR or Trilinos_ROOT: STRING
  * Location of Trilinos install root.
  * Default: None


## Examples of scripts or configuration steps

Examples of steps for building the code can be found in 
[the workflows directory](https://github.com/NimbleSM/NimbleSM/tree/develop/.github/workflows) 
(which builds NimbleSM for testing under different configurations)
or in [the scripts directory](https://github.com/NimbleSM/NimbleSM/tree/develop/scripts/cmake-example-scripts).

### This example provides instructions for setting up NimbleSM along with its dependencies.

**All these commands were used on Ubuntu 24.04.**

---


### Install Ninja Build System

```bat
sudo apt-get update
sudo apt-get install -y ninja-build build-essential cmake gcc-11 g++-11 mpich zlib1g-dev
```

## Step 1: Create Directory Structure and Set Environment Variables

Create the necessary directories for the project and set environment variables:

Set your working directory:
```bat
export WORKDIR=$HOME/dev/NGA/NimbleSM

export VT_SOURCE_DIR=$WORKDIR/vt
export VT_BUILD_DIR=$VT_SOURCE_DIR/build
export VT_INSTALL_DIR=$VT_SOURCE_DIR/install

export FMT_SOURCE_DIR=$WORKDIR/fmt
export FMT_BUILD_DIR=$FMT_SOURCE_DIR/build
export FMT_INSTALL_DIR=$FMT_SOURCE_DIR/install

export SPDLOG_SOURCE_DIR=$WORKDIR/spdlog
export SPDLOG_BUILD_DIR=$SPDLOG_SOURCE_DIR/build
export SPDLOG_INSTALL_DIR=$SPDLOG_SOURCE_DIR/install

export BVH_SOURCE_DIR=$WORKDIR/bvh
export BVH_BUILD_DIR=$BVH_SOURCE_DIR/build
export BVH_INSTALL_DIR=$BVH_SOURCE_DIR/install

export VTK_SOURCE_DIR=$WORKDIR/vtk
export VTK_BUILD_DIR=$VTK_SOURCE_DIR/build
export VTK_INSTALL_DIR=$VTK_SOURCE_DIR/install

export NIMBLESM_SOURCE_DIR=$WORKDIR/NimbleSM
export NIMBLESM_BUILD_DIR=$NIMBLESM_SOURCE_DIR/build
export NIMBLESM_INSTALL_DIR=$NIMBLESM_SOURCE_DIR/install
```
Create the directories:
```bat
mkdir -p $WORKDIR \
$VT_SOURCE_DIR \
$VT_BUILD_DIR \
$VT_INSTALL_DIR \
$FMT_SOURCE_DIR \
$FMT_BUILD_DIR \
$FMT_INSTALL_DIR \
$SPDLOG_SOURCE_DIR \
$SPDLOG_BUILD_DIR \
$SPDLOG_INSTALL_DIR \
$VTK_SOURCE_DIR \
$VTK_BUILD_DIR \
$VTK_INSTALL_DIR \
$BVH_SOURCE_DIR \
$BVH_BUILD_DIR \
$BVH_INSTALL_DIR \
$NIMBLESM_SOURCE_DIR \
$NIMBLESM_BUILD_DIR \
$NIMBLESM_INSTALL_DIR
```

## Step 2: Install Spack and Dependencies

```bat
cd $WORKDIR

git clone -c feature.manyFiles=true https://github.com/spack/spack.git

cd spack/bin

./spack install zlib

source $WORKDIR/spack/share/spack/setup-env.sh
```
Since we want to build all the dependencies with the same version of gcc, **ensure to only have** in ~/.spack/linux/compilers.yaml:

```
compilers:
- compiler:
    spec: gcc@=11.4.0
    paths:
      cc: /usr/bin/gcc-11
      cxx: /usr/bin/g++-11
      f77: /usr/bin/gfortran-11
      fc: /usr/bin/gfortran-11
    flags: {}
    operating_system: ubuntu24.04
    target: x86_64
    modules: []
    environment: {}
    extra_rpaths: []
```
If not do : `spack install gcc@11` and check again the file.

Install dependencies via Spack
```bat
spack install kokkos@4.4.00
spack install arborx@1.4.1
spack install seacas
```

## Step 3: Install external fmt 10.2.1

```bat
cd $FMT_SOURCE_DIR

export FMT_VERSION=10.2.1

wget https://github.com/fmtlib/fmt/archive/refs/tags/${FMT_VERSION}.tar.gz

tar -xzvf ${FMT_VERSION}.tar.gz

rm -rf ${FMT_VERSION}.tar.gz

cmake -S fmt-${FMT_VERSION}/ \
-B $FMT_BUILD_DIR \
-DCMAKE_INSTALL_PREFIX=$FMT_INSTALL_DIR \
-DCMAKE_C_COMPILER=gcc-11 \
-DCMAKE_CXX_COMPILER=g++-11

cd $FMT_BUILD_DIR

make -j$(nproc) install
```
## Step 4: Install VT

```bat
cd $VT_SOURCE_DIR

git clone git@github.com:DARMA-tasking/vt.git

export fmt_DIR=$FMT_BUILD_DIR/lib/cmake/fmt/

CC=gcc-11 \
CXX=g++-11 \
MPICC=$(spack location -i openmpi)/bin/mpicc \
MPICXX=$(spack location -i openmpi)/bin/mpicxx \
VT_EXTENDED_TESTS_ENABLED=0 \
VT_EXTERNAL_FMT=1 \
VT_TRACE_ENABLED=1 \
vt/ci/build_cpp.sh \
$VT_SOURCE_DIR/vt \
$VT_BUILD_DIR \
install
```

## Step 5: Install spdlog 1.13.0

```bat
cd $SPDLOG_SOURCE_DIR

export SPDLOG_VERSION=1.13.0

wget https://github.com/gabime/spdlog/archive/refs/tags/v${SPDLOG_VERSION}.tar.gz

tar -xzvf v${SPDLOG_VERSION}.tar.gz

rm v${SPDLOG_VERSION}.tar.gz

mv spdlog-${SPDLOG_VERSION} spdlog

cmake -S spdlog/ -B $SPDLOG_BUILD_DIR \
-DCMAKE_INSTALL_PREFIX=$SPDLOG_INSTALL_DIR \
-DSPDLOG_FMT_EXTERNAL=ON \
-Dfmt_DIR=$FMT_INSTALL_DIR/lib/cmake/fmt/ \
-DCMAKE_C_COMPILER=gcc-11 \
-DCMAKE_CXX_COMPILER=g++-11

cd $SPDLOG_BUILD_DIR

make -j$(nproc) install
```


## Step 6: Install VTK

```bat
cd $VTK_SOURCE_DIR

export VTK_VERSION=9.3.1

wget https://www.vtk.org/files/release/9.3/VTK-${VTK_VERSION}.tar.gz

tar -xzvf VTK-${VTK_VERSION}.tar.gz

cmake -S VTK-${VTK_VERSION} \
-B $VTK_BUILD_DIR \
-DCMAKE_INSTALL_PREFIX=$VTK_INSTALL_DIR \
-DCMAKE_C_COMPILER=gcc-11 \
-DCMAKE_CXX_COMPILER=g++-11

cd $VTK_BUILD_DIR

make -j$(nproc) install
```

## Step 7: Install distBVH

```bat
cd $BVH_SOURCE_DIR

git clone git@github.com:sandialabs/distBVH.git

export PATH=$(spack location -i openmpi)/bin:$PATH

cmake -S distBVH \
    -B $BVH_BUILD_DIR \
    -DCMAKE_INSTALL_PREFIX=$BVH_INSTALL_DIR \
    -DCMAKE_BUILD_TYPE=Debug \
    -DKokkos_ROOT=$(spack location -i kokkos@4.4.00)/lib/cmake \
    -Dvt_DIR=$VT_BUILD_DIR/vt/install/cmake \
    -DVTK_DIR=$VTK_INSTALL_DIR \
    -DBVH_DEBUG_LEVEL=5 \
    -Dspdlog_DIR=$SPDLOG_INSTALL_DIR/lib/cmake/spdlog \
    -Dfmt_DIR=$FMT_INSTALL_DIR/lib/cmake/fmt \
    -DCMAKE_C_COMPILER=gcc-11 \
    -DCMAKE_CXX_COMPILER=g++-11 \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DCMAKE_VERBOSE_MAKEFILE=OFF

cd $BVH_BUILD_DIR

make -j$(nproc) install
```


## Step 8: Build NimbleSM

Clone the NimbleSM repository:

```bat
cd $NIMBLESM_SOURCE_DIR

git clone git@github.com:NimbleSM/NimbleSM.git
```
Update CMakeLists.txt in $(NIMBLESM_SOURCE_DIR)/unit_tests line 20:

```python
  FetchContent_Declare(
    googletest
    # URL https://github.com/google/googletest/archive/refs/tags/release-1.15.2.tar.gz
    # URL_HASH MD5=ecd1fa65e7de707cd5c00bdac56022cd
    URL https://github.com/google/googletest/releases/download/v1.15.2/googletest-1.15.2.tar.gz       
  )
```

```bat
cmake -S NimbleSM \
    -B $NIMBLESM_BUILD_DIR \
    -DCMAKE_INSTALL_PREFIX=$NIMBLESM_INSTALL_DIR \
    -DKokkos_DIR=$(spack location -i kokkos@4.4.00)/lib/cmake \
    -Dspdlog_DIR=$SPDLOG_INSTALL_DIR/lib/cmake/spdlog \
    -Dbvh_DIR=$BVH_INSTALL_DIR/cmake \
    -DArborX_DIR=$(spack location -i arborx)/lib/cmake/ArborX \
    -DSEACASExodus_DIR=$(spack location -i seacas)/lib/cmake/SEACASExodus \
    -DNimbleSM_ENABLE_KOKKOS=ON \
    -DNimbleSM_ENABLE_BVH=ON \
    -DNimbleSM_ENABLE_ARBORX=ON \
    -DNimbleSM_ENABLE_MPI=ON \
    -DNimbleSM_ENABLE_TRILINOS=OFF \
    -DCMAKE_CXX_COMPILER=g++-11 \
    -DNimbleSM_ENABLE_UNIT_TESTS=ON

cd $NIMBLESM_BUILD_DIR

make -j$(nproc) install
```
You can check and run some tests, for example:
```bat
ctest

cd unit_tests

./NimbleSM_Unit

export test_folder=$NIMBLESM_SOURCE_DIR/NimbleSM/test/contact/sphere_plate_contact

cd $test_folder

mpirun -n 4 $NIMBLESM_BUILD_DIR/src/NimbleSM $test_folder/sphere_plate_contact.in
```

## Different troubleshooting :

unknown command mpi:
```bat
export PATH=$(spack location -i mpi)/bin:$PATH
```

exodiff not found:
```bat
export PATH=$(spack location -i seacas)/bin:$PATH
```


##### [LICENSE](https://github.com/NimbleSM/NimbleSM/blob/develop/LICENSE)

Under the terms of Contract DE-NA0003525 with NTESS,
the U.S. Government retains certain rights in this software.
