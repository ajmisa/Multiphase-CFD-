# Multiphase Computational Fluid Dynamics (MCFD) Projects

This repository contains a series of simulation-based projects focused on multiphase flow modeling, conducted as part of the *Micro Flow Chemistry and Multiphase CFD* course at Eindhoven University of Technology.

Each project addresses a different aspect of CFD modeling in gas–liquid or liquid–solid–gas systems using OpenFOAM or Python-based solvers. Topics range from laminar multiphase convection-reaction models to reactor-scale flow simulations and interface tracking techniques.

## Projects Included

### 1. Multicomponent Convection–Reaction Modeling
Simulates reactive transport in laminar, steady-state conditions using a backward Euler scheme and Newton’s method. General reaction kinetics are implemented with temperature and species coupling.

### 2. Interface Tracking in Multiphase Systems
Applies the Volume of Fluid (VOF) method for simulating dynamic liquid–gas interfaces. The project focuses on mesh refinement, time step sensitivity, and interface accuracy for transient flows.

### 3. Multiphase Reactor Simulation
CFD-based modeling of gas–liquid–solid flow in a lab-scale reactor using Eulerian two-fluid models or simplified hybrid approaches. Boundary conditions and multiphase transport models are tuned for physically meaningful results.

### 4. Validation and Benchmarking
Selected case studies are validated against analytical or literature benchmarks, including mass transfer coefficients, velocity profiles, and reaction conversions.

## Structure

Each folder includes:
- A `README.md` explaining the specific model or simulation
- Source code or solver files
- Output visualizations or postprocessing scripts (Python, ParaView, etc.)
- Final reports summarizing objectives, methodology, results, and conclusions

## Course Information

Micro Flow Chemistry and Multiphase CFD  
TU/e – Eindhoven University of Technology  
Academic Year 2024–2025

## Author

Adam Jordani Misa  
MSc Chemical Engineering – TU/e  
Email: aj.misa@outlook.com
