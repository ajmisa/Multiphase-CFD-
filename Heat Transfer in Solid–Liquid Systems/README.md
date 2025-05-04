# Cooling of a Sphere in a Well-Stirred Liquid

This project models the unsteady cooling of a solid spherical particle immersed in a well-stirred liquid. The system assumes negligible external heat transfer resistance and focuses on the internal diffusion of heat within the sphere. The goal is to solve the dimensionless transient heat conduction equation and compare the numerical results with the analytical solution from classical transport theory.

## Objective

- Formulate the heat conduction model for a sphere with symmetric boundary conditions
- Non-dimensionalize the governing equations and identify relevant parameters (e.g. Biot and Fourier numbers)
- Implement a numerical solution using finite difference or finite volume techniques
- Compare numerical results with the analytical solution from *Transport Phenomena* by Bird, Stewart, and Lightfoot

## Physical Setup

- A sphere of radius \( R \) and initial temperature \( T_1 \) is immersed in a large, well-stirred liquid at \( T_0 \)
- External heat transfer resistance is neglected
- Heat is transferred only by conduction inside the sphere
- The temperature of the liquid changes as a function of heat released from the sphere

## Mathematical Model

- One-dimensional unsteady heat conduction in spherical coordinates
- Boundary conditions:
  - Symmetry at the center of the sphere
  - Convective condition at the surface simplified to a continuity condition with the liquid
- Numerical model implemented using an implicit scheme

## Output and Verification

- The dimensionless temperature of the sphere and liquid are plotted over time
- Results are validated against the analytical solution for comparison
- Grid convergence and time step sensitivity are briefly assessed

