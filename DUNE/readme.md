[Home](../readme)
# DUNE

## Overview

Description from the [DUNE project page](https://www.dune-project.org/):

> "DUNE, the Distributed and Unified Numerics Environment is a modular toolbox for solving partial differential equations (PDEs) with grid-based methods. It supports the easy implementation of methods like Finite Elements (FE), Finite Volumes (FV), and also Finite Differences (FD)."

- UFL based
  
  
## Review 

Using dune-fem module 
Suggested workflow - develop using python-bindings, then transfer to C++ "for efficiency" (they argue that the interface is similar enough to make this trivial), but suggests that the Python bindings can't be used to run production code.

### Gripes
- Suggested (apt-based) installation method didn't work, at least not easily - couldn't find the dune-pdelab package, which seemed to be essential
- dune.fem tutorials require the latest development version, didn't work with `pip install dune.fem`

### Installation (Ubuntu)


**apt**
- abandoned: needs sudo and, although core modules installed ok, 

**pip**
- core modules easy to install, no sudo needed
- additional (dune-fem, dune-fem-dg etc.) modules easy to install
- tutorials work out-of-the-box once you have the correct version


Dependencies
- SciPy
- Eigen (optional)
- PAPI (optional)
- PETSc (optional)
- SIONlib (optional)
- SuiteSparse (optional)



## Links

- [Installation](https://www.dune-project.org/doc/installation/)
- [Core modules on GitLab](https://gitlab.dune-project.org/core/)