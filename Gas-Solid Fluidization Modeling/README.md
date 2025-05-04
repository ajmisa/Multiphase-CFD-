# Comparison of Euler–Euler Models in OpenFOAM

This project investigates gas–solid flow in a fluidized bed using two Euler–Euler modeling approaches available in OpenFOAM: the standard Euler–Euler model and the Kinetic Theory of Granular Flow (KTGF)-enhanced model. The objective is to compare their accuracy, computational efficiency, and agreement with experimental results from the literature.

## Objective

- Implement and simulate gas–solid fluidized beds using:
  - Standard Euler–Euler model (with empirical drag laws)
  - KTGF-based Euler–Euler model (accounting for granular temperature)
- Validate simulations using experimental data from Taghipour et al. (2005)
- Perform sensitivity analysis on time step and spatial resolution
- Evaluate model performance in terms of accuracy and computational cost

## Physical and Numerical Setup

- System: Bubbling fluidized bed filled with glass beads
- Validation case based on Taghipour et al. (2005)
- Mesh: 2D rectangular domain with grid refinement from dx = 0.02 m to dx = 0.005 m
- Solver: `twoPhaseEulerFoam` (OpenFOAM)
- Drag Models: Gidaspow and Wen-Yu (Ergun for dense, Wen-Yu for dilute regimes)
- Boundary conditions and turbulence models defined in `momentumTransportProperties` and `turbulenceProperties.particles`

## Key Results

- **Validation**: Both models reproduce the time-averaged void fraction profiles reasonably well.
- **Temporal Analysis**: Smaller time steps improve agreement with experimental data. The KTGF model is more sensitive to time step size.
- **Spatial Analysis**: Grid refinement improves accuracy but increases computational cost. The Standard EE model performs robustly across all grid sizes.
- **Model Comparison**:
  - Standard EE model shows better stability and slightly better agreement with experiments
  - KTGF offers higher resolution of dynamic structures but requires more careful tuning and computational effort

## Conclusion

- The Standard Euler–Euler model provided more consistent results with lower computational cost
- KTGF-based modeling may improve physical fidelity in dense regimes but requires finer grids and time steps
- Further improvements could be made by exploring turbulence modeling (e.g., LES) or alternative drag coefficients

