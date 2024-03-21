#!/bin/env bash

rm -rf env firedrake-install
curl -O https://raw.githubusercontent.com/firedrakeproject/firedrake/master/scripts/firedrake-install
python3 firedrake-install --venv-name env --install irksome