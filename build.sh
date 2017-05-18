#!/bin/bash

pushd `dirname $0`

mkdir -p build && \
cd build && \
cmake .. && \
make

popd
