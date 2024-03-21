# Adapted from ET implementation; his original preamble included below
#
# Attempt at time-dependent solver of 1D SOL equations, with upwind flux
# based on:
# https://www.firedrakeproject.org/demos/DG_advection.py.html
# irksome implementation follows
# https://www.firedrakeproject.org/Irksome/demos/demo_cahnhilliard.py.html
# and
# https://www.firedrakeproject.org/Irksome/demos/demo_monodomain_FHN.py.html

import firedrake as fd
import irksome as irk
import math
import os.path


def check_dict_empty(d, key_desc, is_fatal=False):
    prefix = "ERROR: Invalid" if is_fatal else "WARNING: Ignoring"
    [print(f"***{prefix} {key_desc} '{k}'***\n") for k in d.keys()]
    if is_fatal:
        raise RuntimeError("check_dict_empty failed with is_fatal=True")


def run_1D_SOL(model_params={}, settings={}, **invalid_kws):
    check_dict_empty(invalid_kws, "keyword")

    # Set defaults for model params that weren't supplied
    nstar = fd.Constant(model_params.pop("nstar", 1.0))
    Temp = fd.Constant(model_params.pop("Temp", 1.0))
    check_dict_empty(model_params, "model parameter")
    # Fix Mach number of 1
    mach_num = fd.Constant(1.0)

    # Set defaults for settings that weren't supplied
    # ic_settings = settings.pop("ICs",)
    meshres = settings.pop("meshres", 200)
    nsteps = settings.pop("nsteps", 250)
    tfinal = settings.pop("tfinal", 10.0)
    check_dict_empty(model_params, "setting")

    # Mesh
    mesh = fd.IntervalMesh(meshres, -1, 1)

    # Function spaces
    V1 = fd.FunctionSpace(mesh, "DG", 1)
    # IMPORTANT: velocity space needs to be continuous
    V2 = fd.VectorFunctionSpace(mesh, "CG", 1)
    V = V1 * V2

    # Set timestep
    t = fd.Constant(0.0)
    dt = fd.Constant(tfinal / nsteps)

    # parameters for irksome
    butcher_tableau = irk.GaussLegendre(2)

    (x,) = fd.SpatialCoordinate(mesh)
    nu = fd.Function(V)
    n, u = fd.split(nu)
    v1, v2 = fd.TestFunctions(V)

    # sonic outflow equilibrium init data
    # nu.sub(0).interpolate((nstar/fd.sqrt(Temp))*(1+fd.sqrt(1-x*x)))
    # nu.sub(1).interpolate((fd.sqrt(Temp)/(x))*(1-fd.sqrt(1-x*x)))

    # Density ICs - Gaussian blob
    mach_num = fd.Constant(model_params.pop("mach_num", 1.0))
    width = 0.1
    nu.sub(0).interpolate(
        (
            1.0
            + 0.2
            * (1 / fd.sqrt(2 * math.pi * width**2))
            * fd.exp(-(x**2) / (2 * width**2))
        )
    )
    # Velocity ICs v=x
    nu.sub(1).interpolate(fd.as_vector([0.0 + 1.0 * x]))

    # Density source - const nstar everywhere
    nstarFunc = fd.Function(V)
    nstarFunc.sub(0).interpolate(nstar + 0.0 * x)

    # Construct a quantity that =u.n for faces with u.n +ve, 0 otherwise
    norm = fd.FacetNormal(mesh)
    u_n = 0.5 * (fd.dot(u, norm) + abs(fd.dot(u, norm)))

    # outflow BC imposed weakly in here
    F = (
        -((irk.Dt(n) * v1) * fd.dx + (n * fd.dot(irk.Dt(u), v2)) * fd.dx)
        + (n * fd.dot(u, fd.grad(v1)) + v1 * nstar) * fd.dx
        + (
            nstar * fd.dot(u, v2)
            + n * u[0] * fd.grad(fd.dot(u, v2))[0]
            + n * u[0] * fd.dot(u, fd.grad(v2[0]))
            + Temp * n * fd.grad(v2[0])[0]
        )
        * fd.dx
        - (v1("+") - v1("-")) * (u_n("+") * n("+") - u_n("-") * n("-")) * fd.dS
        + (u("+")[0] * v2("+")[0] - u("-")[0] * v2("-")[0])
        * (n("+") * u_n("+") - n("-") * u_n("-"))
        * fd.dS
        - fd.conditional(fd.dot(u, norm) > 0, v1 * fd.dot(u, norm) * n, 0.0) * fd.ds
    )
    # Solver params taken from Cahn-Hilliard example cited above
    solver_params = {
        "snes_monitor": None,
        "snes_max_it": 100,
        "snes_linesearch_type": "l2",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "mat_type": "aij",
        "pc_factor_mat_solver_type": "mumps",
    }

    # Dirichlet BCs are needed for boundary velocity; only works for mach number = 1 for now
    u_bc_left = fd.DirichletBC(V.sub(1), fd.as_vector([-mach_num * fd.sqrt(Temp)]), 1)
    u_bc_right = fd.DirichletBC(V.sub(1), fd.as_vector([mach_num * fd.sqrt(Temp)]), 2)

    stepper = irk.TimeStepper(
        F,
        butcher_tableau,
        t,
        dt,
        nu,
        solver_parameters=solver_params,
        bcs=[u_bc_left, u_bc_right],
    )

    # Set up output
    outdir = os.path.dirname(__file__)
    outfile = fd.output.VTKFile(os.path.join(outdir, "SOL_1D_DG_upwind.pvd"))

    # Calc solution vectors
    n_sln = nstar / fd.sqrt(Temp) * (fd.Constant(1) + fd.sqrt(1 - x * x))
    u_sln = fd.sqrt(Temp) * (fd.Constant(1) - fd.sqrt(1 - x * x)) / x

    # Run
    while float(t) < float(tfinal):
        if (float(t) + float(dt)) >= tfinal:
            dt.assign(tfinal - float(t))

        print(f"t = {t}")

        outfile.write(nu.sub(0), nu.sub(1))
        stepper.advance()

        t.assign(float(t) + float(dt))

        # Print fractional errors
        print(f"    n_err = {fd.norm(n-n_sln)/fd.norm(n_sln)}")
        print(f"    u_err = {fd.norm(u[0]-u_sln)/fd.norm(u_sln)}")

    print("done.\n")


run_1D_SOL()
