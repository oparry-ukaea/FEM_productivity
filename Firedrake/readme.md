[Home](../readme)
# Firedrake

## Overview
Description from the [Firedrake project page](https://www.firedrakeproject.org/):

> "Firedrake is an automated system for the solution of partial differential equations using the finite element method (FEM). Firedrake uses sophisticated code generation to provide mathematicians, scientists, and engineers with a very high productivity way to create sophisticated high performance simulations."


## Review 

Tight integration with PETSc - effectively write callback functions for it in C, translated from user Python
No postproc tools at first glance - just generates vtk and relies on ParaView.
Surprisingly difficult to infer from the docs that its doing code gen
Parallelism all handled under-the-hood; 

### Installation

System install of several dependencies...

- PETSc

## Links

- [GitHub page](https://github.com/firedrakeproject/firedrake)
- [Irkesome](https://github.com/firedrakeproject/Irksome): Time-stepping package used to generate Runge Kutta methods for Firedrake
- [firedrakeproject/ufl](https://github.com/firedrakeproject/ufl): A fork of [FEniCS/UFL](https://github.com/fenics/ufl) - used to provide a DSL for building Firedrake solvers
- [Installation](https://www.firedrakeproject.org/download.html)