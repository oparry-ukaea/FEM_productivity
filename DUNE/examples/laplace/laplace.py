import numpy as np

from dune.grid import structuredGrid

gridView = structuredGrid([0, 0], [1, 1], [20, 20])

from dune.fem.space import lagrange

space = lagrange(gridView, order=1)
u_h = space.interpolate(0, name="u_h")

from ufl import (
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    dx,
    grad,
    inner,
    dot,
    sin,
    cos,
    pi,
)

x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)

f = (8 * pi**2 + 1) * cos(2 * pi * x[0]) * cos(2 * pi * x[1])
a = (inner(grad(u), grad(v)) + u * v) * dx
l = f * v * dx


from dune.fem import assemble

mat, rhs = assemble(a == l)

from scipy.sparse.linalg import spsolve as solver

A = mat.as_numpy
b = rhs.as_numpy
y = u_h.as_numpy
y[:] = solver(A, b)

u_h.plot()


from dune.fem import integrate

exact = cos(2 * pi * x[0]) * cos(2 * pi * x[1])
e_h = u_h - exact
squaredErrors = integrate([e_h**2, inner(grad(e_h), grad(e_h))])
print("L^2 and H^1 errors:", [np.sqrt(e) for e in squaredErrors])
