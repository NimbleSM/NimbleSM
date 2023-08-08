FROM pierrpebay/nimblesm:pr-2 AS build_nimble-stage
# need to be repeated in the new stage
ARG NimbleSM_ENABLE_MPI
ARG NimbleSM_ENABLE_KOKKOS
ARG NimbleSM_ENABLE_TRILINOS
ARG NimbleSM_ENABLE_UQ
ARG NimbleSM_ENABLE_ARBORX
ARG NimbleSM_CONFIGURATION_NAME

# Add current source dir into the image
COPY . /opt/src/NimbleSM
RUN mkdir -p /opt/build/NimbleSM

# Build using the spack environment we created
RUN bash /opt/src/NimbleSM/ci/build.sh

FROM build_nimble-stage AS test-stage

RUN bash /opt/src/NimbleSM/ci/test.sh

FROM scratch AS export-stage
COPY --from=test-stage /tmp/artifacts /tmp/artifacts
