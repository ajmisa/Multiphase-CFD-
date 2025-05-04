# Numerical Simulation of Droplet Impact on Inclined Surfaces

This project investigates the dynamics of a water droplet impacting and spreading on a hydrophilic inclined surface using the Volume of Fluid (VOF) method in OpenFOAM. The simulations are benchmarked against experimental data and focus on comparing two contact angle boundary condition models: constant and dynamic.

## Objective

- Simulate droplet impact on a 64.25° inclined hydrophilic surface
- Implement and compare constant and dynamic contact angle models
- Analyze the effect of tangential velocity on spreading dynamics
- Validate numerical predictions against experimental results from Sookran (2018)

## Methodology

- **Geometry**: A 20 mm × 20 mm × 21.5 mm inclined domain created via `blockMesh`
- **VOF Solver**: Tracks the liquid–air interface throughout the impact and spreading
- **Contact Angle Models**:
  - `constantAlphaContactAngle`: static contact angle (88°)
  - `dynamicAlphaContactAngle`: varying with contact line velocity (88° advancing, 32° receding)
- **Physical Properties**:
  - Water droplet diameter: 2.5 mm
  - Impact velocity: 1.35 m/s normal, 2.9 m/s tangential
  - Surface tension: 72.8 mN/m

## Key Results

- **Constant Contact Angle**:
  - Captures basic spreading behavior but underpredicts length
  - Overestimates splashing and shows blank regions in the droplet body
- **Dynamic Contact Angle**:
  - Closer agreement with experiment in spreading phase
  - Smoother droplet shape and improved wetting behavior
- **Limitations**:
  - Both models affected by coarse grid resolution
  - Need for grid and timestep refinement to fully resolve thin liquid layers

## Conclusion

- The dynamic contact angle model better represents droplet dynamics, but further resolution improvements are needed.
- Visual agreement with experiment supports the use of dynamic wetting models in inclined impact simulations.
- Future work should focus on mesh refinement and time step sensitivity studies to ensure model accuracy.

