#!/usr/bin/env bash

set -e
. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

pushd /opt/build/NimbleSM
success_flag=0
ctest --output-on-failure || success_flag=1
# We collect the test logs for exporting
echo "this is the Success flag: ${success_flag}"
mkdir -p /tmp/artifacts/
cp /opt/build/NimbleSM/Testing/Temporary/LastTest.log /tmp/artifacts/
cp /opt/build/NimbleSM/test/*/*/*.log /tmp/artifacts/
echo ${success_flag} > /tmp/artifacts/success_flag.txt
ls /tmp/artifacts
# Simply output the LastTest.log to screen
cat /opt/build/NimbleSM/Testing/Temporary/LastTest.log
popd
