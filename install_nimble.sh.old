#!/bin/bash

echo "Usage: $0 <WORKDIR> <nproc>"
if [ "$#" -ne 2 ]; then
    WORKDIR=$HOME/NimbleSM
    nproc=6
    echo "Default values are set for the build : folder = $WORKDIR, nproc = $nproc"
else
    WORKDIR=$1
    nproc=$2
fi

sudo apt install git

echo "Testing SSH connection to GitHub..."

ssh -T git@github.com 2>&1 | grep "successfully authenticated" > /dev/null

if [ $? -eq 0 ]; then
    echo "Success: SSH key is correctly linked to GitHub."
else
    # Path to the .ssh directory
    SSH_DIR="$HOME/.ssh"

    # Check if the .ssh directory exists
    if [ -d "$SSH_DIR" ]; then
        echo ".ssh directory exists."
    else
        echo ".ssh directory does not exist, creating it..."
        mkdir -p "$SSH_DIR"
        chmod 700 "$SSH_DIR"
    fi

    # Check if an SSH key pair already exists
    if [ -f "$SSH_DIR/id_rsa" ] && [ -f "$SSH_DIR/id_rsa.pub" ]; then
        echo "An existing SSH key pair found:"
        ls "$SSH_DIR/id_rsa" "$SSH_DIR/id_rsa.pub"

        echo "To add this key to GitHub, follow these steps:"
        echo "1. Copy your public SSH key:"
        echo ""
        echo "   cat ~/.ssh/id_rsa.pub"
        echo ""
        echo "2. Go to GitHub and log in."
        echo "3. In the upper-right corner of any page, click your profile photo, then click 'Settings'."
        echo "4. In the 'Access' section of the sidebar, click 'SSH and GPG keys'."
        echo "5. Click 'New SSH key' or 'Add SSH key'."
        echo "6. Paste your public key into the 'Key' field."
        echo "7. Click 'Add SSH key'."
        echo "8. Confirm your GitHub password if prompted."
        
        # Continue with installations or other tasks
        echo "Proceeding with installations..."
        
        # Add your installation commands here
        
    else
        echo "No SSH key found, generating a new SSH key..."
        ssh-keygen -t rsa -b 4096 -f "$SSH_DIR/id_rsa" -N ""
        echo "SSH key successfully generated."
        
        # Display instructions for adding the newly generated key to GitHub
        echo "To add this newly generated key to GitHub, follow these steps:"
        echo "1. Copy your public SSH key:"
        echo ""
        echo "   cat ~/.ssh/id_rsa.pub"
        echo ""
        echo "2. Go to GitHub and log in."
        echo "3. In the upper-right corner of any page, click your profile photo, then click 'Settings'."
        echo "4. In the 'Access' section of the sidebar, click 'SSH and GPG keys'."
        echo "5. Click 'New SSH key' or 'Add SSH key'."
        echo "6. Paste your public key into the 'Key' field."
        echo "7. Click 'Add SSH key'."
        echo "8. Confirm your GitHub password if prompted."
    fi
    exit 1
fi


echo "Working Directory: $WORKDIR"
echo "Number of Cores: $nproc"

sudo apt-get update

sudo apt install gcc-11
sudo apt install g++-11
sudo apt install gfortran-11

sudo apt-get install -y ninja-build build-essential cmake gcc-11 g++-11 mpich zlib1g-dev




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

mkdir -p $WORKDIR \
$VT_SOURCE_DIR $VT_BUILD_DIR $VT_INSTALL_DIR \
$FMT_SOURCE_DIR $FMT_BUILD_DIR $FMT_INSTALL_DIR \
$SPDLOG_SOURCE_DIR $SPDLOG_BUILD_DIR $SPDLOG_INSTALL_DIR \
$VTK_SOURCE_DIR $VTK_BUILD_DIR $VTK_INSTALL_DIR \
$BVH_SOURCE_DIR $BVH_BUILD_DIR $BVH_INSTALL_DIR \
$NIMBLESM_SOURCE_DIR $NIMBLESM_BUILD_DIR $NIMBLESM_INSTALL_DIR

cd $WORKDIR

git clone -c feature.manyFiles=true https://github.com/spack/spack.git

source $WORKDIR/spack/share/spack/setup-env.sh

spack install zlib

rm ~/.spack/linux/compilers.yaml

echo "compilers:
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
    extra_rpaths: []" > ~/.spack/linux/compilers.yaml

spack install kokkos@4.4.00
spack install arborx@1.4.1
spack install seacas


cd $FMT_SOURCE_DIR
export FMT_VERSION=10.2.1
wget https://github.com/fmtlib/fmt/archive/refs/tags/${FMT_VERSION}.tar.gz
tar -xzvf ${FMT_VERSION}.tar.gz
rm -rf ${FMT_VERSION}.tar.gz
cmake -S fmt-${FMT_VERSION}/ -B $FMT_BUILD_DIR \
-DCMAKE_INSTALL_PREFIX=$FMT_INSTALL_DIR \
-DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11
cd $FMT_BUILD_DIR
make -j$(nproc) install


cd $VT_SOURCE_DIR
git clone git@github.com:DARMA-tasking/vt.git
export fmt_DIR=$FMT_INSTALL_DIR/lib/cmake/fmt/
CC=gcc-11 CXX=g++-11 \
MPICC=$(spack location -i openmpi)/bin/mpicc \
MPICXX=$(spack location -i openmpi)/bin/mpicxx \
VT_EXTENDED_TESTS_ENABLED=0 VT_EXTERNAL_FMT=1 VT_TRACE_ENABLED=1 \
vt/ci/build_cpp.sh $VT_SOURCE_DIR/vt $VT_BUILD_DIR install


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


cd $VTK_SOURCE_DIR
export VTK_VERSION=9.3.1
wget https://www.vtk.org/files/release/9.3/VTK-${VTK_VERSION}.tar.gz
tar -xzvf VTK-${VTK_VERSION}.tar.gz
cmake -S VTK-${VTK_VERSION} -B $VTK_BUILD_DIR \
-DCMAKE_INSTALL_PREFIX=$VTK_INSTALL_DIR \
-DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11
cd $VTK_BUILD_DIR
make -j$(nproc) install


cd $BVH_SOURCE_DIR
git clone git@github.com:sandialabs/distBVH.git
export PATH=$(spack location -i openmpi)/bin:$PATH
cmake -S distBVH -B $BVH_BUILD_DIR \
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


cd $NIMBLESM_SOURCE_DIR
git clone git@github.com:NimbleSM/NimbleSM.git
export PATH=$(spack location -i mpi)/bin:$PATH
export PATH=$(spack location -i seacas)/bin:$PATH
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