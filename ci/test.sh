. /opt/spack/share/spack/setup-env.sh
spack env activate nimble

pushd /opt/build/NimbleSM
testlist=(
  "ctest --output-on-failure | tee ctest-output.log"
)

for this in "${testlist[@]}"; do
  $this || flag=1
done

popd
exit $flag
