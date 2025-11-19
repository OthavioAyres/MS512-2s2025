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

    Returns a tuple (u, iterations, final_error), where:
      - u is the solution array,
      - iterations is the number of Jacobi steps performed (up to convergence
        or nitmax), and
      - final_error is the last value of the relative L2 error from convcheck.
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

    iterations = 0
    erro = 0.0

    # Jacobi iterations
    for k in range(1, nitmax + 1):
        uk = jacobi_step(u, b, h)
        erro, converged = convcheck(u, uk, eps)

        iterations = k

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

    return u, iterations, erro


def run_example_8_2_8():
    """
    Reproduce Example 8.2.8: run the Jacobi solver for h = 1/10, 1/20, 1/40
    and print a table with (h, matrix dimension, iterations to convergence).
    """
    # Grid sizes chosen so that h = 1 / (m - 1) is 1/10, 1/20, 1/40
    m_values = [11, 21, 41,81]
    nitmax = 1_000_000
    eps = 1e-8

    results = []

    for m in m_values:
        h = 1.0 / (m - 1)
        # We suppress file output here, since we are only interested
        # in the iteration counts for the convergence study.
        u, iterations, erro = poisson_jacobi(
            m=m,
            nitmax=nitmax,
            eps=eps,
            write_output=False,
        )

        # Number of interior unknowns in each direction is (m-2),
        # so the total number of unknowns is (m-2)^2.
        matrix_dim = (m - 2) ** 2
        h_label = f"1/{m - 1}"

        results.append((h_label, matrix_dim, iterations))
    print()
    print("Example 8.2.8 - Jacobi method convergence study")
    header_h = "h"
    header_dim = "Matrix dimension"
    header_it = "Iterations to convergence"
    print(f"{header_h:>8} {header_dim:>18} {header_it:>26}")

    for h_label, matrix_dim, iterations in results:
        print(f"{h_label:>8} {matrix_dim:18d} {iterations:26d}")


if __name__ == "__main__":
    # Running this module directly reproduces Example 8.2.8:
    # it executes the Jacobi method for h = 1/10, 1/20, 1/40, 1/80 and
    # prints a table with the matrix dimension and iterations
    # to convergence for each grid.
    run_example_8_2_8()


