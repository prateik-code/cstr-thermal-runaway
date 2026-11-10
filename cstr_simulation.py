import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


Volume = 100.0             # Reactor liquid volume (Liters)
flow_rate = 10.0              # Volumetric flow rate in and out (L/min)
C_H2O2 = 5.0            # Feed concentration of H2O2 (mol/L)
T_feed = 298.15       # Temperature of incoming feed fluid (Kelvin)

# Kinetic parameters (Arrhenius relationship)
A_factor = 2.0e11     # Pre-exponential frequency factor (1/min)
Ea_dividedby_R = 9000.0    # Activation Energy divided by Universal Gas Constant (Kelvin)

# Thermodynamic properties of the aqueous system
delta_H = -98200.0    # Enthalpy of reaction (J/mol)
rho_c = 4184.0       # Volumetric heat capacity: density * Cp (J / (L * K))

# Cooling jacket specifications
UA = 3000.0           # Overall heat transfer coefficient * Surface area (J / (min * K))
T_coolant = 290.0     # Temperature of cooling fluid in jacket (Kelvin)


def cstr_odes(t, state):
    """
    Computes instantaneous derivatives [dCA/dt, dT/dt] given current system state.
    state[0]: Concentration of H2O2 inside reactor (mol/L)
    state[1]: Temperature inside reactor (Kelvin)
    """
    C_A = state[0]
    T = state[1]

    # Reaction rate constant calculated based on current temperature
    k = A_factor * np.exp(-Ea_dividedby_R / T)

    # 1. Mass Balance: Accumulation = Inflow - Outflow - Reaction Consumption
    dCA_dt = (flow_rate / Volume) * (C_H2O2 - C_A) - (k * C_A)

    # 2. Energy Balance: Accumulation = Inflow - Outflow + Generation - Cooling
    heat_generated = -delta_H * k * C_A                    # Joules generated per Liter-min
    cooling_removed = UA * (T - T_coolant) / Volume        # Joules removed per Liter-min
    
    dT_dt = (flow_rate / Volume) * (T_feed - T) + (heat_generated - cooling_removed) / rho_c

    return [dCA_dt, dT_dt]


time_span = (0, 60)                     # Simulation duration (0 to 60 minutes)
initial_state = [0.0, 298.15]           # Pure water at room temperature in the tank at the start
time_points = np.linspace(0, 60, 600)   # 600 evaluation timestamps(evenly distributed)

# Solved initial value problem(used Runge-Kutta method (RK45))
solution = solve_ivp(
    fun=cstr_odes,
    t_span=time_span,
    y0=initial_state,
    t_eval=time_points,
    method='RK45'
)

# Extract trajectory arrays
t = solution.t
C_A_profile = solution.y[0]
T_profile = solution.y[1]


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

# H2O2 Concentration Dynamics
ax1.plot(t, C_A_profile, color='#1f77b4', linewidth=2.2, label='[H2O2] in Tank')
# ax1.axhline(C_H2O2, color='gray', linestyle='--', alpha=0.6, label='Inflow Concentration (C_H2O2)')
ax1.set_ylabel('Concentration (mol/L)', fontsize=11, fontweight='bold')
ax1.set_title('Non-Isothermal CSTR Dynamic Response: H2O2 Decomposition', fontsize=13, fontweight='bold')
ax1.grid(True, linestyle=':', alpha=0.7)
ax1.legend(loc='upper right')

# Reactor Thermal Dynamics
ax2.plot(t, T_profile, color='#d62728', linewidth=2.2, label='Reactor Core Temp (T)')
ax2.axhline(T_coolant, color='blue', linestyle=':', alpha=0.7, label='Coolant Jacket Temp (T_cool)')
ax2.axhline(373.15, color='black', linestyle='--', linewidth=1.5, label='Aqueous Boiling Limit (373.15 K)')
ax2.set_xlabel('Time (minutes)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Temperature (K)', fontsize=11, fontweight='bold')
ax2.grid(True, linestyle=':', alpha=0.7)
ax2.legend(loc='upper right')

plt.tight_layout()
plt.savefig("cstr_dynamics_plot.png", dpi=300, bbox_inches="tight")
plt.show()