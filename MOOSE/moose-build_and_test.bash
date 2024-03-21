#!/bin/env bash

# Only allow this script to be sourced, not executed
(return 0 2>/dev/null) && sourced=1 || sourced=0

if [ $sourced -eq 0 ]; then
    echo "moose-build_and_test.bash must be sourced, not executed"
fi

# Activate conda env - depends on conda init having been run already
if [ ! command -v conda &> /dev/null ]; then
    echo "conda isn't on PATH"
    exit 1
fi
conda activate moose

# Clone repo
git clone https://github.com/idaholab/moose.git repo -b master


# Build and run tests
cd ./repo/test
make -j $np
./run_tests -j 8
cd -