import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

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

# Dimensionless Parameters
timesteps = 2000
tau_max = 0.2  # Set maximum dimensionless time to 0.2
delta_tau = tau_max / timesteps
tau = np.linspace(0, tau_max, timesteps + 1)
B_values = [0.5, 1, 2, 4, 10, 1000]

plt.figure(figsize=(10, 6))

for B in B_values:
    # Find roots for the given B
    b_k_values = find_roots(B)
    # Compute Theta_f using the derived formula
    Theta_f = np.zeros_like(tau)
    for k in range(1, len(b_k_values)):
        b_k = b_k_values[k]
        Theta_f += 6 * B * np.exp(-b_k**2 * tau) / (9 * (1 + B) + B**2 * b_k**2)
    Theta_f += B / (1 + B)
    # Plot the results
    plt.plot(tau, (1 + B) * (1 - Theta_f), label=f'B = {B}')

# Finalize plot
plt.xlabel(r'$\frac{\alpha_s t}{R^2}$')
plt.ylabel(r'$(1 + B) \left( 1 - \frac{T - T_0}{T_i - T_0} \right)$')
plt.title('Dimensionless Temperature Evolution for Different B Values')
plt.legend()
plt.grid()
plt.xlim(-0.01, 0.2)  # Set x-axis limits
plt.ylim(0, 1.0)  # Set y-axis limits
plt.show()