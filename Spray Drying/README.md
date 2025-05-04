# Effect of Turbulence Closures on Spray Jet Modeling

This project investigates the impact of turbulence model selection on the simulation of spray jet behavior in a scaled spray dryer using the Euler–Lagrangian framework in OpenFOAM. The study compares two widely used turbulence closures—**k-epsilon** and **k-omega**—in their ability to predict flow dynamics, jet spreading, and velocity profiles.

## Objective

- Reproduce the experimental jet configuration from Pawar (2014)
- Simulate jet development using OpenFOAM's `sprayFoam` solver
- Compare simulation results using the k-epsilon and k-omega turbulence models
- Evaluate model performance based on:
  - Velocity field comparison with experimental data
  - Error norm analysis (L1, L2, Linf)
  - Visual symmetry and spreading patterns

## Simulation Setup

- Geometry: Rectangular spray chamber (0.2 × 0.2 × 0.3 m)
- Inlet nozzle: Diameter 0.025 m
- Fluid: Air (N₂) as continuous phase, water as dispersed phase
- Spray injection: Rosin–Rammler distribution for droplet sizes (0.2–0.7 mm)
- Reynolds number: 2.07 × 10⁵
- Mass flow rate: 4.07 kg/s
- Time duration: 0.01 s (transient phase)
- Solver: `sprayFoam` (OpenFOAM)

## Turbulence Models Compared

| Model      | Near-Wall Accuracy | Jet Penetration | Computational Cost | Recommendation               |
|------------|--------------------|------------------|--------------------|-------------------------------|
| k-epsilon  | Low                | Overpredicts     | Lower              | Suitable for far-field flows |
| k-omega    | High               | Accurate         | Moderate           | Better for shear/jet regions |

## Key Results

- **Velocity Profiles**: k-omega matched experimental peak velocity and spreading more closely.
- **Error Analysis**:
  - Both models showed first-order convergence with timestep refinement.
  - k-omega consistently produced lower L1, L2, and Linf error norms.
- **Visual Symmetry**: Both models showed symmetric jet behavior, but k-omega produced sharper and more accurate jet core decay.
- **Conclusion**: k-omega is preferred for near-wall and shear-layer accuracy, while k-epsilon performs adequately in bulk regions.

## Limitations and Recommendations

- Finer grid refinement and longer simulation times are recommended.
- Advanced turbulence (e.g. LES) and droplet interaction models (e.g. TAB or KHRT) may improve accuracy.
- Additional validation against broader experimental datasets is needed.

