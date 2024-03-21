#!/bin/env bash

if [ ! command -v conda &> /dev/null ]; then
    echo "conda isn't on PATH"
    exit 1
fi

conda update --all --yes
conda config --add channels https://conda.software.inl.gov/public
conda create -n moose moose-dev=2024.03.15