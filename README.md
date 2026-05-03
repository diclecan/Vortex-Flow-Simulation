# Taylor-Green Vortex Simulation

This project simulates the two-dimensional Taylor-Green vortex using the vorticity-streamfunction formulation of the incompressible Navier-Stokes equations.

The Taylor-Green vortex is a classical test case in numerical fluid dynamics because an analytical solution is known. Therefore, the numerical solution can be compared directly with the analytical reference solution in order to evaluate accuracy, convergence behaviour, energy decay and enstrophy decay.

## Project Overview

The simulation solves the vorticity transport equation

\[
\frac{\partial \omega}{\partial t} + \mathbf{u} \cdot \nabla \omega = \nu \Delta \omega
\]

together with the streamfunction-vorticity relation

\[
\Delta \psi = -\omega
\]

where:

- `ω` is the vorticity,
- `ψ` is the streamfunction,
- `u` and `v` are the velocity components,
- `ν` is the kinematic viscosity.

The velocity field is computed from the streamfunction as

\[
u = \frac{\partial \psi}{\partial y}, \qquad
v = -\frac{\partial \psi}{\partial x}.
\]

## Numerical Method

The implementation uses:

- finite differences on a periodic Cartesian grid,
- central difference operators for first and second derivatives,
- sparse matrices for efficient storage,
- LU factorization for repeatedly solving the Poisson equation,
- time integration of the vorticity transport equation,
- comparison with the analytical Taylor-Green vortex solution.

The computational domain is

\[
[0, 2\pi] \times [0, 2\pi]
\]

with periodic boundary conditions.

## Repository Structure

```text
.
├── Sim_TG_Wirbel.m          # Main MATLAB script
├── Simulation.m             # Simulation class with numerical operators and solver methods
├── assets/
│   ├── Lagrange_Partikel_vergleich_ana.gif
│   └── Stromfunktion_3D_vergleich_ana.gif
├── report/
│   └── Wirbelstromprojekt_Bericht.pdf
└── README.md
