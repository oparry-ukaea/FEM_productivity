#!/bin/env bash
miniforge_dir=$SOFTWARE/miniforge

curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$miniforge_dir"
PATH="$miniforge_dir/bin:$PATH"

conda init --all