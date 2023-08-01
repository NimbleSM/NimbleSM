#!/usr/bin/env bash

set -x
set -e
. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

pushd /opt/build/NimbleSM
ret_code=0
ctest --output-on-failure || ret_code=$?
# We collect the test logs for exporting
echo "ctest returned: $ret_code"
mkdir -p /tmp/artifacts/
cp /opt/build/NimbleSM/Testing/Temporary/LastTest.log /tmp/artifacts/
cp /opt/build/NimbleSM/test/*/*/*.log /tmp/artifacts/
echo ${ret_code} > /tmp/artifacts/success_flag.txt
ls /tmp/artifacts
popd

exit $ret_code
