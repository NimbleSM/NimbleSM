FROM pierrpebay/nimblesm:latest as build_nimble-stage

# Add current source dir into the image
COPY . /opt/src/NimbleSM
RUN mkdir -p /opt/build/NimbleSM

# Build using the spack environment we created
RUN bash /opt/src/NimbleSM/ci/build-serial.sh

FROM build_nimble-stage as test-stage

RUN bash /opt/src/NimbleSM/ci/test.sh

FROM scratch as export-stage
COPY --from=test-stage /tmp/artifacts /tmp/artifacts
