import numpy as np
import matplotlib.pyplot as plt

# numerical solution functions and setup
def initialize_solid_liquid(gridsteps):
    """Initialize temperature for solid and liquid."""
    T_s = np.zeros(gridsteps)  # Solid initially zero
    return T_s

def calculate_fo(alpha, delta_t, delta_xi):
    """Calculate the Fourier number."""
    return alpha * delta_t / delta_xi**2

def update_temperature_solid(T_s_old, T_l, Fo, gridsteps):
    """Update solid temperature using the implicit scheme."""
    T_s_new = np.zeros_like(T_s_old)
    delta_xi = 1 / (gridsteps - 1)  # Calculate spatial step size

    for i in range(gridsteps):
        if i == 0:
            # Boundary condition at ξ = 0 (symmetry condition)
            T_s_new[i] = T_s_old[i + 1]
        elif i == gridsteps - 1:
            # Boundary condition at ξ = 1 (liquid-solid interaction)
            T_s_new[i] = T_l
        else:
            # Compute the position ξ_i, ξ_{i+1/2}, and ξ_{i-1/2}
            xi = i * delta_xi
            xi_plus_half = xi + delta_xi / 2
            xi_minus_half = xi - delta_xi / 2

            # Correctly calculate the term involving ξ^2
            xi2_term = Fo * ((xi_plus_half**2) * T_s_old[i + 1]
                             - (xi_plus_half**2 + xi_minus_half**2) * T_s_old[i]
                             + (xi_minus_half**2) * T_s_old[i - 1])
            
            # Update T_s for the current grid point
            T_s_new[i] = T_s_old[i] + (1 / (xi**2 if xi != 0 else 1)) * xi2_term

    return T_s_new

# Analytical solution function using Fourier series
def temperature_distribution_fourier(xi, tau, num_terms=50):
    """Calculate the dimensionless temperature distribution using Fourier series."""
    theta = np.zeros_like(xi)
    for n in range(1, num_terms + 1):
        lambdan = (n * np.pi)  # nth root approximation
        term = (2 * (-1)**n / (n * np.pi * (xi + 1e-10))) * np.sin(n * np.pi * xi) * np.exp(-n**2 * np.pi**2 * tau)
        theta += term
    return 1 + theta

# Simulation parameters
gridsteps = 51
alpha = 1  # Thermal diffusivity
tau_max_values = [0.01, 0.04, 0.1, 0.2, 0.4]  # Different tau_max values
timesteps = 10000  # Number of timesteps

# Initialize plot
plt.figure(figsize=(12, 8))

for tau_max in tau_max_values:
    delta_tau = tau_max / timesteps
    delta_xi = 1 / (gridsteps - 1)
    Fo = calculate_fo(alpha, delta_tau, delta_xi)
    
    # Initialize temperatures
    T_s = initialize_solid_liquid(gridsteps)
    solid_temps = [T_s.copy()]

    # Time evolution (Numerical)
    for _ in range(timesteps):
        T_s = update_temperature_solid(T_s, 1.0, Fo, gridsteps)
        solid_temps.append(T_s.copy())

    # Analytical solution preparation
    xi = np.linspace(0, 1, gridsteps)
    tau = tau_max
    analytical_temp = temperature_distribution_fourier(xi, tau)

    # Plot final profile for this tau_max
    plt.plot(xi, solid_temps[-1], '--', label=f'Numerical τ_max = {tau_max}')
    plt.plot(xi[1:], analytical_temp[1:], label=f'Analytical τ_max = {tau_max}')  # Exclude xi=0 for analytical

# Finalize plot
plt.xlabel(r'$\xi$')
plt.ylabel(r'$\theta_{\text{solid}}$')
plt.title('Solid Temperature Profile for Different $\tau_{\text{max}}$')
plt.legend()
plt.grid()
plt.show()