#!/bin/env bash

poisson_dir="examples/poisson"
poisson_script="https://docs.fenicsproject.org/dolfinx/main/python/_downloads/b94ac7be61dc3726ca331afd20f195d2/demo_poisson.py"

mkdir -p "$poisson_dir"
wget "$poisson_script" -O "$poisson_dir/demo_poisson.py"