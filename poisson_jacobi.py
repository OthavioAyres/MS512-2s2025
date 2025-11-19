#!/usr/bin/env python3
import math


def gbc(m, h):
    """
    Build boundary condition array g[j][i], analogous to subroutine gbc in Fortran.
    Indices: i -> x, j -> y. Python indices go from 0 to m-1.
    """
    # Initialize g with zeros
    g = [[0.0 for _ in range(m)] for _ in range(m)]

    for j in range(m):
        y = h * j
        for i in range(m):
            x = h * i

            # The order of these 'if's matches the Fortran code.
            if i == 0:
                g[j][i] = 0.0
            if i == m - 1:
                g[j][i] = y
            if j == 0:
                g[j][i] = (x - 1.0) * math.sin(x)
            if j == m - 1:
                g[j][i] = x * (2.0 - x)

    return g


def jacobi_step(u, b, h):
    """
    Perform one Jacobi iteration.
    Returns uk (previous u) so we can use it in convcheck, like in the Fortran code.
    u is updated in-place to the new values.
    """
    m = len(u)
    # Copy current u to uk (previous iteration)
    uk = [row[:] for row in u]

    # Update only interior points (1 .. m-2)
    for j in range(1, m - 1):
        for i in range(1, m - 1):
            u[j][i] = (
                uk[j][i - 1]   # left
                + uk[j][i + 1]  # right
                + uk[j - 1][i]  # down
                + uk[j + 1][i]  # up
                + (h ** 2) * b[j][i]
            ) / 4.0

    return uk


def convcheck(u, uk, eps):
    """
    Check convergence using the relative L2 norm of the update:
    erro = sqrt( sum (u - uk)^2 / sum u^2 ).
    Returns (erro, converged_bool).
    """
    m = len(u)
    eu = 0.0
    ed = 0.0

    for j in range(m):
        for i in range(m):
            val = u[j][i]
            diff = val - uk[j][i]
            eu += val * val
            ed += diff * diff

    # Just in case eu is zero (shouldn't happen here because of nonzero boundary),
    # guard against division by zero:
    if eu == 0.0:
        erro = 0.0
    else:
        erro = math.sqrt(ed / eu)

    converged = (erro < eps)
    return erro, converged


def write_solution(filename, u, h):
    """
    Write solution in the same format as the Fortran code (x, y, u(x,y)).
    This lets you use gnuplot with g-sol.plt (which expects 'fort.100').
    """
    m = len(u)
    with open(filename, "w") as f:
        for j in range(m):
            y = h * j
            for i in range(m):
                x = h * i
                f.write(f"{x:13.6f}  {y:13.6f}  {u[j][i]:15.8e}\n")


def poisson_jacobi(m=21, nitmax=1_000_000, eps=1e-8, write_output=True):
    """
    Main driver: solve Laplace's equation on [0,1]x[0,1] with the same
    boundary conditions and parameters as the Fortran program.
    """
    h = 1.0 / (m - 1)

    # Initialize arrays: u (solution), b (source term, here zero everywhere)
    u = [[0.0 for _ in range(m)] for _ in range(m)]
    b = [[0.0 for _ in range(m)] for _ in range(m)]

    # Apply boundary conditions
    g = gbc(m, h)
    # Copy g into u
    for j in range(m):
        for i in range(m):
            u[j][i] = g[j][i]

    # Jacobi iterations
    for k in range(1, nitmax + 1):
        uk = jacobi_step(u, b, h)
        erro, converged = convcheck(u, uk, eps)

        if converged:
            print(f"Converged in {k:8d} iterations, Error = {erro:15.10e}")
            break
    else:
        # If we did not break from the loop (no convergence within nitmax)
        print(f"Did not converge within {nitmax} iterations. Last error = {erro:15.10e}")

    # Optionally write solution for plotting with gnuplot
    if write_output:
        write_solution("fort.100", u, h)
        print("Solution written to 'fort.100' (use gnuplot 'g-sol.plt' to visualize).")

    return u


if __name__ == "__main__":
    # Parameters chosen to match the Fortran code
    poisson_jacobi(m=21, nitmax=1_000_000, eps=1e-8, write_output=True)