#!/bin/bash

install_vt(){
	set -ex

	export CC=gcc-11 
	export CXX=g++-11 
	export MPICC=$(spack location -i openmpi)/bin/mpicc
	export MPICXX=$(spack location -i openmpi)/bin/mpicxx
	export fmt_DIR=$FMT_INSTALL_DIR/lib/cmake/fmt/
	export VT_EXTENDED_TESTS_ENABLED=0
	export VT_EXTERNAL_FMT=1
	export VT_TRACE_ENABLED=1
	
	source_dir=${1}
	build_dir=${2}

	# Dependency versions, when fetched via git.
	checkpoint_rev=develop

	if test "${VT_DOXYGEN_ENABLED:-0}" -eq 1
	then
	    token=${3}
	else
	    target=${3:-install}
	fi

	if [ -z ${4} ]; then
	    dashj=""
	else
	    dashj="-j ${4}"
	fi

	if hash ccache &>/dev/null
	then
	    use_ccache=true
	fi

	if test "$use_ccache"
	then
	    { echo -e "===\n=== ccache statistics before build\n==="; } 2>/dev/null
	    ccache -s
	else
	    { echo -e "===\n=== ccache not found, compiling without it\n==="; } 2>/dev/null
	fi

	mkdir -p "${build_dir}"
	pushd "${build_dir}"

	# Match `nvcc_wrapper` and also a path ending with 'nvcc_wrapper'
	case $CXX in
	    *nvcc_wrapper)
		NVCC_WRAPPER_DEFAULT_COMPILER="$(which g++-"$(echo "${HOST_COMPILER}" | cut -d- -f2)")" \
		&& export NVCC_WRAPPER_DEFAULT_COMPILER;;
	esac

	if test "${VT_KOKKOS_ENABLED:-0}" -eq 1
	then
	  echo "The variable VT_KOKKOS_ENABLED is set."

	  if test -d "kokkos"
	  then
	    rm -Rf kokkos
	  fi

	  git clone -b master https://github.com/kokkos/kokkos.git
	  export KOKKOS_DIR=$PWD/kokkos
	  export KOKKOS_BUILD=${build_dir}/kokkos/build
	  export KOKKOS_INSTALL="$KOKKOS_BUILD/install"
	  mkdir -p "$KOKKOS_BUILD"
	  cd "$KOKKOS_BUILD"
	  cmake -G "${CMAKE_GENERATOR:-Ninja}" \
		-DCMAKE_INSTALL_PREFIX="$KOKKOS_INSTALL" \
		"$KOKKOS_DIR"
	  cmake --build . ${dashj} --target install
	  cd "${build_dir}"
	fi

	if test -d "checkpoint"
	then
	    rm -Rf checkpoint
	fi

	if test -d "vt-tv"
	then
	    rm -Rf vt-tv
	fi

	if test -d "${source_dir}/lib/checkpoint"
	then
	    { echo "Checkpoint already in lib... not downloading, building, and installing"; } 2>/dev/null
	else
	    if test "${VT_DOXYGEN_ENABLED:-0}" -eq 1
	    then
		cd "${source_dir}/lib"
		git clone -b "${checkpoint_rev}" --depth 1 https://github.com/DARMA-tasking/checkpoint.git
		cd -
	    else
		git clone -b "${checkpoint_rev}" --depth 1 https://github.com/DARMA-tasking/checkpoint.git
		export CHECKPOINT=$PWD/checkpoint
		export MAGISTRATE_BUILD=${build_dir}/checkpoint
		mkdir -p "$MAGISTRATE_BUILD"
		cd "$MAGISTRATE_BUILD"
		mkdir build
		cd build
		cmake -G "${CMAKE_GENERATOR:-Ninja}" \
		      -DCMAKE_INSTALL_PREFIX="$MAGISTRATE_BUILD/install" \
		      -DKokkos_ROOT="$KOKKOS_INSTALL" \
		      "$CHECKPOINT"
		cmake --build . ${dashj} --target install
	    fi
	fi

	if test "${VT_TV_ENABLED}" -eq 1
	then
	    if test -d "${source_dir}/lib/vt-tv"
	    then
		{ echo "vt-tv already in lib... not downloading"; } 2>/dev/null
	    else
		cd "${source_dir}/lib"
		vt_tv_rev="1.5.0"
		git clone -b "${vt_tv_rev}" --depth 1 https://github.com/DARMA-tasking/vt-tv.git
		cd -
	    fi
	fi

	if test "${VT_ZOLTAN_ENABLED:-0}" -eq 1
	then
	    export Zoltan_DIR=${ZOLTAN_DIR:-""}
	fi

	if test "${VT_CI_BUILD:-0}" -eq 1
	then
	    git config --global --add safe.directory "${source_dir}"
	fi

	export VT=${source_dir}
	export VT_BUILD=${build_dir}/vt
	mkdir -p "$VT_BUILD"
	cd "$VT_BUILD"
	rm -Rf ./*
	cmake -G "${CMAKE_GENERATOR:-Ninja}" \
	      -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
	      -Dvt_test_trace_runtime_enabled="${VT_TRACE_RUNTIME_ENABLED:-0}" \
	      -Dvt_lb_enabled="${VT_LB_ENABLED:-1}" \
	      -Dvt_trace_enabled="${VT_TRACE_ENABLED:-0}" \
	      -Dvt_trace_only="${VT_BUILD_TRACE_ONLY:-0}" \
	      -Dvt_doxygen_enabled="${VT_DOXYGEN_ENABLED:-0}" \
	      -Dvt_mimalloc_enabled="${VT_MIMALLOC_ENABLED:-0}" \
	      -Dvt_asan_enabled="${VT_ASAN_ENABLED:-0}" \
	      -Dvt_ubsan_enabled="${VT_UBSAN_ENABLED:-0}" \
	      -Dvt_werror_enabled="${VT_WERROR_ENABLED:-0}" \
	      -Dvt_pool_enabled="${VT_POOL_ENABLED:-1}" \
	      -Dvt_build_extended_tests="${VT_EXTENDED_TESTS_ENABLED:-1}" \
	      -Dvt_zoltan_enabled="${VT_ZOLTAN_ENABLED:-0}" \
	      -Dvt_tv_enabled="${VT_TV_ENABLED:-0}" \
	      -Dvt_production_build_enabled="${VT_PRODUCTION_BUILD_ENABLED:-0}" \
	      -Dvt_unity_build_enabled="${VT_UNITY_BUILD_ENABLED:-0}" \
	      -Dvt_diagnostics_enabled="${VT_DIAGNOSTICS_ENABLED:-1}" \
	      -Dvt_diagnostics_runtime_enabled="${VT_DIAGNOSTICS_RUNTIME_ENABLED:-0}" \
	      -Dvt_fcontext_enabled="${VT_FCONTEXT_ENABLED:-0}" \
	      -Dvt_fcontext_build_tests_examples="${VT_FCONTEXT_BUILD_TESTS_EXAMPLES:-0}" \
	      -Dvt_rdma_tests_enabled="${VT_RDMA_TESTS_ENABLED:-1}" \
	      -Dvt_code_coverage="${VT_CODE_COVERAGE:-0}" \
	      -DMI_INTERPOSE:BOOL=ON \
	      -DMI_OVERRIDE:BOOL=ON \
	      -Dvt_mpi_guards="${VT_MPI_GUARD_ENABLED:-0}" \
	      -DMPI_EXTRA_FLAGS="${MPI_EXTRA_FLAGS:-}" \
	      -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}" \
	      -DMPI_C_COMPILER="${MPICC:-mpicc}" \
	      -DMPI_CXX_COMPILER="${MPICXX:-mpicxx}" \
	      -DCMAKE_CXX_COMPILER="${CXX:-c++}" \
	      -DCMAKE_C_COMPILER="${CC:-cc}" \
	      -DCMAKE_EXE_LINKER_FLAGS="${CMAKE_EXE_LINKER_FLAGS:-}" \
	      -Dmagistrate_ROOT="$MAGISTRATE_BUILD/install" \
	      -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-}" \
	      -DCMAKE_INSTALL_PREFIX="$VT_BUILD/install" \
	      -Dvt_ci_build="${VT_CI_BUILD:-0}" \
	      -Dvt_ci_generate_lb_files="${VT_CI_TEST_LB_SCHEMA:-0}" \
	      -Dvt_debug_verbose="${VT_DEBUG_VERBOSE:-0}" \
	      -Dvt_tests_num_nodes="${VT_TESTS_NUM_NODES:-}" \
	      -Dvt_external_fmt="${VT_EXTERNAL_FMT:-0}" \
	      -Dfmt_DIR="${fmt_DIR}" \
	      -DLIBUNWIND_ROOT="${LIBUNWIND_ROOT:-/usr}" \
	      -Dvt_no_color_enabled="${VT_NO_COLOR_ENABLED:-0}" \
	      -DCMAKE_CXX_STANDARD="${CMAKE_CXX_STANDARD:-17}" \
	      -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS:-0}" \
	      "$VT"
	cmake_conf_ret=$?

	if test "${VT_DOXYGEN_ENABLED:-0}" -eq 1
	then
	    MCSS=$PWD/m.css
	    GHPAGE=$PWD/DARMA-tasking.github.io
	    git clone "https://${token}@github.com/DARMA-tasking/DARMA-tasking.github.io"
	    git clone https://github.com/mosra/m.css
	    cd m.css
	    git checkout master
	    cd ../

	    "$MCSS/documentation/doxygen.py" Doxyfile-mcss
	    cp -R docs "$GHPAGE"
	    cd "$GHPAGE"
	    git config --global user.email "jliffla@sandia.gov"
	    git config --global user.name "Jonathan Lifflander"
	    git add docs
	    git commit -m "Update docs (auto-build)"
	    git push origin master
	elif test "${VT_CI_BUILD:-0}" -eq 1
	then
	    # Generate output file with compilation warnings and errors

	    GENERATOR=$(cmake -L . | grep USED_CMAKE_GENERATOR:STRING | cut -d"=" -f2)
	    OUTPUT="$VT_BUILD"/compilation_errors_warnings.out
	    OUTPUT_TMP="$OUTPUT".tmp

	    # Because of the problem with new lines in Azure pipelines, all of them will be
	    # converted to this unique delimiter
	    DELIMITER="-=-=-=-"

	    WARNS_ERRS=""

	    # Unfortunately Ninja doesn't output compilation warnings and errors to stderr
	    # so it needs special treatment
	    if test "$GENERATOR" = "Ninja"
	    then
		# To easily tell if compilation of given file succeeded special progress bar is used
		# (controlled by variable NINJA_STATUS)
		export NINJA_STATUS="[ninja][%f/%t] "
		time cmake --build . ${dashj} --target "${target}" | tee "$OUTPUT_TMP"
		compilation_ret=${PIPESTATUS[0]}
		sed -i '/ninja: build stopped:/d' "$OUTPUT_TMP"

		# Now every line that doesn't start with [ninja][number]/[number] is an error or a warning
		WARNS_ERRS=$(grep -Ev '^(\[ninja\]\[[[:digit:]]+\/[[:digit:]]+\])|(--) .*$' "$OUTPUT_TMP" || true)
	    elif test "$GENERATOR" = "Unix Makefiles"
	    then
		# Gcc outputs warnings and errors to stderr, so there's not much to do
		time cmake --build . ${dashj} --target "${target}" 2> >(tee "$OUTPUT_TMP")
		compilation_ret=$?
		WARNS_ERRS=$(cat "$OUTPUT_TMP")
	    fi

	    # Convert new lines and redirect to an output file
	    WARNS_ERRS=${WARNS_ERRS//$'\n'/$DELIMITER}
	    echo "$WARNS_ERRS" > "$OUTPUT"
	else
	    time cmake --build . ${dashj} --target "${target}"
	    compilation_ret=$?
	fi

	if test "$use_ccache"
	then
	    { echo -e "===\n=== ccache statistics after build\n==="; } 2>/dev/null
	    ccache -s
	fi

	# Exit with error code if there was any
	if test "$cmake_conf_ret" -ne 0
	then
	    echo "There was an error during CMake configuration"
	    exit "$cmake_conf_ret"
	elif test "$compilation_ret" -ne 0
	then
	    echo "There was an error during compilation"
	    exit "$compilation_ret"
	fi
}

show_help() {
    echo "Usage: $0 <WORKDIR> <nproc> [options]"
    echo "For example ./install_nimble.sh 6 --all"
    echo "Options:"
    echo "  --all          Install all dependencies"
    echo "In case of the modification of this script by yourself to change a version of dependence, you can run the following options to rebuilt   "
    echo "  --fmt          Install FMT"
    echo "  --spdlog       Install spdlog"
    echo "  --vt           Install DARMA-vt"
    echo "  --vtk          Install VTK"
    echo "  --bvh          Install BVH"
    echo "  --nimble       Install NimbleSM"
    echo "  -h, --help     Show this help message and exit"
}

# Default WORKDIR and nproc
if [ "$#" -lt 1 ]; then
    show_help
    exit 1
else
    WORKDIR=$(pwd "$0")/dependencies
    nproc=$1
    shift 1
fi

# Parse optional dependencies
INSTALL_LIST=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --all)
            INSTALL_LIST=(fmt spdlog magistrate vt vtk bvh nimble)
            shift ;;
        --update) INSTALL_LIST+=("update"); shift ;;    
        --spack) INSTALL_LIST+=("spack"); shift ;;    
        --fmt) INSTALL_LIST+=("fmt"); shift ;;
        --spdlog) INSTALL_LIST+=("spdlog"); shift ;;
        --vt) INSTALL_LIST+=("vt"); shift ;;
        --vtk) INSTALL_LIST+=("vtk"); shift ;;
        --bvh) INSTALL_LIST+=("bvh"); shift ;;
        --nimble) INSTALL_LIST+=("nimble"); shift ;;
        -h|--help)
            show_help
            exit 0 ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1 ;;
    esac
done

sudo apt install git -y

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

#Create all the repository
echo "Working Directory: $WORKDIR"
echo "Number of Cores: $nproc"

export FMT_SOURCE_DIR=$WORKDIR/fmt
export FMT_BUILD_DIR=$FMT_SOURCE_DIR/build
export FMT_INSTALL_DIR=$FMT_SOURCE_DIR/install

export VT_SOURCE_DIR=$WORKDIR/vt
export VT_BUILD_DIR=$VT_SOURCE_DIR/build
export VT_INSTALL_DIR=$VT_BUILD_DIR/install

export SPDLOG_SOURCE_DIR=$WORKDIR/spdlog
export SPDLOG_BUILD_DIR=$SPDLOG_SOURCE_DIR/build
export SPDLOG_INSTALL_DIR=$SPDLOG_SOURCE_DIR/install

export BVH_SOURCE_DIR=$WORKDIR/bvh
export BVH_BUILD_DIR=$BVH_SOURCE_DIR/build
export BVH_INSTALL_DIR=$BVH_SOURCE_DIR/install

export VTK_SOURCE_DIR=$WORKDIR/vtk
export VTK_BUILD_DIR=$VTK_SOURCE_DIR/build
export VTK_INSTALL_DIR=$VTK_SOURCE_DIR/install

export NIMBLESM_SOURCE_DIR=$WORKDIR
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


if [[ " ${INSTALL_LIST[@]} " =~ " update " ]]; then
	sudo apt-get update -y

	sudo apt install -y gcc-11
	sudo apt install -y g++-11 
	sudo apt install -y gfortran-11 

	sudo apt-get install -y ninja-build build-essential cmake gcc-11 g++-11 mpich zlib1g-dev
	sudo apt install -y libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev
fi


if [[ " ${INSTALL_LIST[@]} " =~ " spack " ]]; then
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
fi

if [[ " ${INSTALL_LIST[@]} " =~ " fmt " ]]; then
	echo "==================================="
	echo "=         FMT INSTALLATION        ="
	echo "==================================="

	cd $FMT_SOURCE_DIR
	export FMT_VERSION=10.2.1
	wget https://github.com/fmtlib/fmt/archive/refs/tags/${FMT_VERSION}.tar.gz
	tar -xzf ${FMT_VERSION}.tar.gz
	rm -rf ${FMT_VERSION}.tar.gz
	cmake -S fmt-${FMT_VERSION}/ -B $FMT_BUILD_DIR \
	-DCMAKE_INSTALL_PREFIX=$FMT_INSTALL_DIR \
	-DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11
	cd $FMT_BUILD_DIR
	make -j$(nproc) install
fi


if [[ " ${INSTALL_LIST[@]} " =~ " vt " ]]; then
	echo "==================================="
	echo "=         VT INSTALLATION         ="
	echo "==================================="
	
	source $WORKDIR/spack/share/spack/setup-env.sh
	
	echo $VT_SOURCE_DIR
	
	cd $VT_SOURCE_DIR
	export VT_VERSION=1.5.0
	git clone git@github.com:DARMA-tasking/vt.git
	cd vt
	git checkout tags/${VT_VERSION}
	
	#export fmt_DIR=$FMT_INSTALL_DIR/lib/cmake/fmt/
	
	
	install_vt $VT_SOURCE_DIR/vt $VT_BUILD_DIR install $nproc
	
	
fi
  


if [[ " ${INSTALL_LIST[@]} " =~ " spdlog " ]]; then
	echo "==================================="
	echo "=       SPDLOG INSTALLATION       ="
	echo "==================================="
	
	cd $SPDLOG_SOURCE_DIR
	export SPDLOG_VERSION=1.13.0
	wget https://github.com/gabime/spdlog/archive/refs/tags/v${SPDLOG_VERSION}.tar.gz
	tar -xzf v${SPDLOG_VERSION}.tar.gz
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

fi

if [[ " ${INSTALL_LIST[@]} " =~ " vtk " ]]; then
	echo "==================================="
	echo "=         VTK INSTALLATION        ="
	echo "==================================="

	cd $VTK_SOURCE_DIR
	export VTK_VERSION=9.3.1
	wget https://www.vtk.org/files/release/9.3/VTK-${VTK_VERSION}.tar.gz
	tar -xzf VTK-${VTK_VERSION}.tar.gz
	rm VTK-${VTK_VERSION}.tar.gz
	cmake -S VTK-${VTK_VERSION} -B $VTK_BUILD_DIR \
	-DCMAKE_INSTALL_PREFIX=$VTK_INSTALL_DIR \
	-DCMAKE_C_COMPILER=gcc-11 -DCMAKE_CXX_COMPILER=g++-11
	cd $VTK_BUILD_DIR
	make -j$(nproc) install
fi

if [[ " ${INSTALL_LIST[@]} " =~ " bvh " ]]; then
	echo "==================================="
	echo "=         BVH INSTALLATION        ="
	echo "==================================="

	cd $BVH_SOURCE_DIR
	git clone git@github.com:sandialabs/distBVH.git
	export PATH=$(spack location -i openmpi)/bin:$PATH
	cmake -S distBVH -B $BVH_BUILD_DIR \
	-DCMAKE_INSTALL_PREFIX=$BVH_INSTALL_DIR \
	-DCMAKE_BUILD_TYPE=Debug \
	-DKokkos_ROOT=$(spack location -i kokkos@4.4.00)/lib/cmake \
	-Dvt_DIR=$VT_INSTALL_DIR/cmake \
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

fi


if [[ " ${INSTALL_LIST[@]} " =~ " nimble " ]]; then
	echo "==================================="
	echo "=      NIMBLESM INSTALLATION      ="
	echo "==================================="

	cd $NIMBLESM_SOURCE_DIR
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
fi



