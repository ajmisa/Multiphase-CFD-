import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

# Define the root equation for analytical solution
def root_equation(b_k, B):
    return np.tan(b_k) - 3 * b_k / (3 + B * b_k**2)

# Find roots for the analytical solution
def find_roots(B):
    initial_guesses = np.linspace(0, 1000, 10000)
    roots = []
    for guess in initial_guesses:
        sol = root(root_equation, guess, args=(B,))
        if sol.success and np.abs(root_equation(sol.x[0], B)) < 1e-5:
            if len(roots) == 0 or all(np.abs(sol.x[0] - r) > 1e-2 for r in roots):
                roots.append(sol.x[0])
    return np.array(roots)

# Analytical and numerical parameters
timesteps = 2000
tau_max = 0.2
delta_tau = tau_max / timesteps
tau = np.linspace(0, tau_max, timesteps + 1)
B_values = [0.5, 1, 2, 4, 10]
gridsteps = 51
alpha = 1
delta_xi = 1 / (gridsteps - 1)

# Initialize plot
plt.figure(figsize=(12, 8))

for B in B_values:
    # Analytical Solution
    b_k_values = find_roots(B)
    Theta_f = np.zeros_like(tau)
    for k in range(1, len(b_k_values)):
        b_k = b_k_values[k]
        Theta_f += 6 * B * np.exp(-b_k**2 * tau) / (9 * (1 + B) + B**2 * b_k**2)
    Theta_f += B / (1 + B)
    analytical = (1 + B) * (1 - Theta_f)

    # Numerical Solution
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
    
    # Plot comparison
    plt.plot(tau, analytical, label=f'Analytical B={B}', linestyle='--')
    plt.plot(tau, numerical, label=f'Numerical B={B}', linestyle='-')

# Finalize plot
plt.xlabel(r'$\tau$')
plt.ylabel(r'$(1 + B) \cdot (1 - \theta)$')
plt.title('Comparison of Analytical and Numerical Solutions')
plt.legend()
plt.grid()
plt.xlim(0, 0.2)
plt.ylim(0, 1.0)
plt.show()
