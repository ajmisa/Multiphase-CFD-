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
timesteps = 1100
tau_max = 0.2
delta_tau = tau_max / timesteps
tau = np.linspace(0, tau_max, timesteps + 1)
B_values = [0.5, 1, 2, 4, 10]
gridsteps = 51
alpha = 1
delta_xi = 1 / (gridsteps - 1)

# Initialize error storage
error_metrics = {B: {} for B in B_values}

for B in B_values:
    # Analytical solution
    b_k_values = find_roots(B)
    Theta_f = np.zeros_like(tau)
    for k in range(1, len(b_k_values)):
        b_k = b_k_values[k]
        Theta_f += 6 * B * np.exp(-b_k**2 * tau) / (9 * (1 + B) + B**2 * b_k**2)
    Theta_f += B / (1 + B)
    analytical = (1 + B) * (1 - Theta_f)

    # Numerical solution
    Fo = alpha * delta_tau / delta_xi**2
    T_s = np.zeros(gridsteps)
    T_l = 1.0
    theta_l_all = [T_l]

    for _ in range(int(tau_max / delta_tau)):
        T_s_new = np.zeros_like(T_s)
        for i in range(gridsteps):
            if i == 0:
                T_s_new[i] = T_s[i+1]
            elif i == gridsteps-1:
                T_s_new[i] = T_l
            else:
                xi = i * delta_xi
                xi_plus_half = xi + delta_xi / 2
                xi_minus_half = xi - delta_xi / 2
                xi2_term = Fo * ((xi_plus_half**2) * T_s[i + 1]
                                 - (xi_plus_half**2 + xi_minus_half**2) * T_s[i]
                                 + (xi_minus_half**2) * T_s[i - 1])
                T_s_new[i] = T_s[i] + (1 / (xi**2 if xi != 0 else 1)) * xi2_term
        T_s = T_s_new.copy()
        interaction_term = (-3 * delta_tau / (B * 2 * delta_xi)) * (
            3 * T_s[-1] - 4 * T_s[-2] + T_s[-3])
        T_l += interaction_term
        theta_l_all.append(T_l)
    numerical = [(1 + B) * (1 - theta_l) for theta_l in theta_l_all]

    # Calculate errors
    absolute_error = np.abs(np.array(analytical) - np.array(numerical))
    relative_error = np.abs(absolute_error / np.array(analytical))
    rmse = np.sqrt(np.mean(absolute_error**2))
    mae = np.mean(absolute_error)
    l_infinity = np.max(absolute_error)  # Maximum absolute error

    # Store errors
    error_metrics[B] = {"RMSE": rmse, "MAE": mae, "L-infinity": l_infinity}

# Plot RMSE, MAE, and L-infinity for different B values
B_vals = list(error_metrics.keys())
rmse_vals = [error_metrics[B]["RMSE"] for B in B_vals]
mae_vals = [error_metrics[B]["MAE"] for B in B_vals]
l_infinity_vals = [error_metrics[B]["L-infinity"] for B in B_vals]

plt.figure(figsize=(10, 6))
plt.plot(B_vals, rmse_vals, marker='o', label="RMSE")
plt.plot(B_vals, mae_vals, marker='s', label="MAE")
plt.plot(B_vals, l_infinity_vals, marker='^', label="L-infinity")
plt.xlabel(r"$B$")
plt.ylabel("Error")
plt.title("Error Metrics for Different $B$ Values")
plt.xscale('log')  # Optional: Log scale for better visualization
plt.yscale('log')
plt.legend()
plt.grid()
plt.show()
