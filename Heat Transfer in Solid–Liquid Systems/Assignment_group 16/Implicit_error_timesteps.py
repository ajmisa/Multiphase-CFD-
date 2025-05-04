import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root


# Root and analytical solution functions
def root_equation(b_k, B):
    return np.tan(b_k) - 3 * b_k / (3 + B * b_k**2)

def find_roots(B):
    initial_guesses = np.linspace(0, 1000, 10000)
    roots = []
    for guess in initial_guesses:
        sol = root(root_equation, guess, args=(B,))
        if sol.success and np.abs(root_equation(sol.x[0], B)) < 1e-5:
            if len(roots) == 0 or all(np.abs(sol.x[0] - r) > 1e-2 for r in roots):
                roots.append(sol.x[0])
    return np.array(roots)

# Parameters
tau_max = 0.2
gridsteps = 51
alpha = 1
delta_xi = 1 / (gridsteps - 1)
B = 2  # Fixed B for this analysis

# Analytical solution
timesteps_analytical = 2000  # High-resolution analytical timesteps
delta_tau_analytical = tau_max / timesteps_analytical
tau_analytical = np.linspace(0, tau_max, timesteps_analytical + 1)
b_k_values = find_roots(B)
Theta_f = np.zeros_like(tau_analytical)
for k in range(1, len(b_k_values)):
    b_k = b_k_values[k]
    Theta_f += 6 * B * np.exp(-b_k**2 * tau_analytical) / (9 * (1 + B) + B**2 * b_k**2)
Theta_f += B / (1 + B)
analytical_solution = (1 + B) * (1 - Theta_f)

# Different timesteps for numerical solution
timesteps_list = [100, 200, 500, 1000, 2000, 5000, 10000]

# Error metrics storage
error_metrics = {"Timesteps": [], "RMSE": [], "MAE": [], "L-infinity": []}

# Loop over different timesteps
for timesteps in timesteps_list:
    delta_tau = tau_max / timesteps
    tau = np.linspace(0, tau_max, timesteps + 1)
    Fo = alpha * delta_tau / delta_xi**2

    # Initialize temperatures
    T_s = np.zeros(gridsteps)
    T_l = 1.0
    theta_l_all = [T_l]

    # Setup implicit matrix
    A = np.zeros((gridsteps, gridsteps))
    for i in range(1, gridsteps - 1):
        xi = i * delta_xi
        xi_plus_half = xi + delta_xi / 2
        xi_minus_half = xi - delta_xi / 2
        A[i, i - 1] = -Fo * xi_minus_half**2 / xi**2
        A[i, i] = 1 + Fo * (xi_plus_half**2 + xi_minus_half**2) / xi**2
        A[i, i + 1] = -Fo * xi_plus_half**2 / xi**2

    A[0, 0] = 1  # Boundary condition at ξ = 0
    A[gridsteps - 1, gridsteps - 1] = 1  # Boundary condition at ξ = 1

    # Time evolution
    for _ in range(timesteps):
        b = T_s.copy()
        b[0] = T_s[1]  # Boundary condition at ξ = 0
        b[gridsteps - 1] = T_l  # Boundary condition at ξ = 1

        T_s_new = np.linalg.solve(A, b)
        T_s = T_s_new.copy()

        # Update liquid temperature
        interaction_term = (-3 * delta_tau / (B * 2 * delta_xi)) * (
            3 * T_s[-1] - 4 * T_s[-2] + T_s[-3])
        T_l += interaction_term
        theta_l_all.append(T_l)

    numerical_solution = [(1 + B) * (1 - theta_l) for theta_l in theta_l_all]

    # Interpolate analytical solution to match numerical tau points
    analytical_interp = np.interp(tau, tau_analytical, analytical_solution)

    # Calculate errors
    absolute_error = np.abs(np.array(analytical_interp) - np.array(numerical_solution))
    rmse = np.sqrt(np.mean(absolute_error**2))
    mae = np.mean(absolute_error)
    l_infinity = np.max(absolute_error)

    # Store errors
    error_metrics["Timesteps"].append(timesteps)
    error_metrics["RMSE"].append(rmse)
    error_metrics["MAE"].append(mae)
    error_metrics["L-infinity"].append(l_infinity)

# Plot error metrics vs. timesteps
plt.figure(figsize=(10, 6))
plt.plot(error_metrics["Timesteps"], error_metrics["RMSE"], marker='o', label="RMSE")
plt.plot(error_metrics["Timesteps"], error_metrics["MAE"], marker='s', label="MAE")
plt.plot(error_metrics["Timesteps"], error_metrics["L-infinity"], marker='^', label="L-infinity")
plt.xlabel("Timesteps")
plt.ylabel("Error")
plt.title("Error Metrics vs. Numerical Timesteps (Implicit Method)")
plt.legend()
plt.xscale('log')  # Log scale for timesteps
plt.yscale('log')  # Log scale for error
plt.grid()
plt.show()
