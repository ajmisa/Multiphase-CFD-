import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
timesteps = 1100
tau_max = 0.2
gridsteps = 51
alpha = 1  # Example thermal diffusivity
delta_tau = tau_max / timesteps  # Time step size
delta_xi = 1 / (gridsteps - 1)  # Spatial step size

# Calculate the Fourier number
Fo = alpha * delta_tau / delta_xi**2

# Values of B to consider
B_values = [0.5, 1, 2, 4, 10, 1000]

# Initialize plot
plt.figure(figsize=(10, 6))

for B in B_values:
    # Initialize temperatures for solid and liquid
    T_s = np.zeros(gridsteps)  # Solid initially zero
    T_l = 1.0  # Liquid initially one (normalized)

    # Time evolution of temperatures
    theta_l_all = [T_l]  # Track liquid temperature

    # Implicit matrix setup
    A = np.zeros((gridsteps, gridsteps))
    for i in range(1, gridsteps - 1):
        xi = i * delta_xi
        xi_plus_half = xi + delta_xi / 2
        xi_minus_half = xi - delta_xi / 2
        A[i, i - 1] = -Fo * xi_minus_half**2 / xi**2
        A[i, i] = 1 + Fo * (xi_plus_half**2 + xi_minus_half**2) / xi**2
        A[i, i + 1] = -Fo * xi_plus_half**2 / xi**2

    A[0, 0] = 1
    A[gridsteps - 1, gridsteps - 1] = 1

    for _ in range(timesteps):
        # Setup the right-hand side (RHS) of the implicit equation
        b = T_s.copy()
        b[0] = T_s[1]  # Boundary condition at ξ = 0 (symmetry)
        b[gridsteps - 1] = T_l  # Boundary condition at ξ = 1 (interaction)

        # Solve the linear system A * T_s_new = b
        T_s_new = np.linalg.solve(A, b)
        T_s = T_s_new.copy()  # Update the old array with new values

        # Update liquid temperature
        interaction_term = (-3 * delta_tau / (B * 2 * delta_xi)) * (
            3 * T_s[-1] - 4 * T_s[-2] + T_s[-3])
        T_l += interaction_term
        theta_l_all.append(T_l)  # Track liquid temperature

    # Calculate the rescaled values for (1 + B) * (1 - theta_liquid)
    rescaled_theta_l = [(1 + B) * (1 - theta_l) for theta_l in theta_l_all]

    # Plot for the current B
    tau = np.linspace(0, timesteps * delta_tau, timesteps + 1)
    plt.plot(tau, rescaled_theta_l, label=f'B = {B}')

# Finalize plot
plt.xlabel(r'$\tau$')
plt.ylabel(r'$(1 + B) \cdot (1 - \theta_{\mathrm{liquid}})$')
plt.title('Rescaled Liquid Temperature Evolution vs. $\tau$ for Different B Values')
plt.legend()
plt.grid()
plt.xlim(-0.01, 0.2)  # Set x-axis limits
plt.ylim(0, 1.0)  # Set y-axis limits
plt.show()