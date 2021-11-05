#!/usr/bin/env bash

#
# //@HEADER
# // ************************************************************************
# //
# //                                NimbleSM
# //                             Copyright 2018
# //   National Technology & Engineering Solutions of Sandia, LLC (NTESS)
# //
# // Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
# // retains certain rights in this software.
# //
# // Redistribution and use in source and binary forms, with or without
# // modification, are permitted provided that the following conditions are
# // met:
# //
# // 1. Redistributions of source code must retain the above copyright
# // notice, this list of conditions and the following disclaimer.
# //
# // 2. Redistributions in binary form must reproduce the above copyright
# // notice, this list of conditions and the following disclaimer in the
# // documentation and/or other materials provided with the distribution.
# //
# // 3. Neither the name of the Corporation nor the names of the
# // contributors may be used to endorse or promote products derived from
# // this software without specific prior written permission.
# //
# // THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY EXPRESS OR IMPLIED
# // WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# // MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
# // NO EVENT SHALL NTESS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# // TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# // PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# // LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# // NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# // SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# //
# // Questions?  Contact David Littlewood (djlittl@sandia.gov)
# //
# // ************************************************************************
# //@HEADER
#

set -exo pipefail

if test $# -lt 1
then
    echo "usage: ./$0 <seacas-version>"
    exit 1
fi

mkdir -p /opt/seacas
pushd /opt/seacas

seacas_version="$1"
seacas_tar_name="v${seacas_version}.tar.gz"

echo "${seacas_version}"
echo "${seacas_tar_name}"

wget "http://github.com/gsjaardema/seacas/archive/${seacas_tar_name}"

tar xzf "${seacas_tar_name}"

pushd "seacas-${seacas_version}"
export ACCESS=$(pwd)

export CRAY=OFF
export MPI=ON
export JOBS=8
export STATIC=YES
export USE_PROXY=NO
export CC=gcc
export CXX=g++
export FORTRAN=OFF
export X11=OFF

./install-tpl.sh

mkdir build
pushd build
../cmake-config
make -j8 install

popd