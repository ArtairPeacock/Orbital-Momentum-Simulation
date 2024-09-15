import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11 # Gravitational constant
MS = 1.989e30 # Mass of the Sun in kg
ME = 5.972e24 # Mass of the Sun in kg
r = 1.0 # Initial distance of Earth to the Sun in AU (1 AU)
v_K = np.sqrt(G * MS / r) # Keplerian velocity in m/s

# Function to calculate acceleration
def AccelerationGravity(r_vec):
	r_mag = np.linalg.norm(r_vec)
	return (-G * MS / r_mag**3) * r_vec
	
# Function implementing Runge-Kutta 2nd order method for orbit simulation
def SimulateOrbit_rk2(dt):
	num_steps = int(total_time / dt)
	time = np.linspace(0, total_time, num_steps + 1)
	
	r_vec = np.array([r, 0, 0]) # Initial position vector in AU
	v_vec = np.array([0, v_K, 0]) # Initial velocity vector in m/s
	
	angular_momentum_z = [] # List to store angular momentum
	
	for t in time:
		# Calculate angular momentum
		L = np.cross(r_vec, MS * v_vec)
		angular_momentum_z.append(L[2]) # Append z-component of angular momentum
		
		k1v = AccelerationGravity(r_vec)
		k1r = v_vec
		
		k2r = v_vec + k1v * dt
		k2v = AccelerationGravity(r_vec + k1r * dt)
		
		r_vec += 0.5 * (k1r + k2r) * dt
		v_vec += 0.5 * (k1v + k2v) * dt
		
	return np.array(angular_momentum_z)
	
# Function implementing Euler method for orbit simulation (with corrected angular
momentum calculation)
def SimulateOrbit_euler(dt):
	num_steps = int(total_time / dt)
	time = np.linspace(0, total_time, num_steps + 1)
	
	r_vec = np.array([r, 0, 0]) # Initial position vector in AU
	v_vec = np.array([0, v_K, 0]) # Initial velocity vector in m/s
	
	angular_momentum_z = [] # List to store angular momentum
	
	for t in time:
		# Calculate angular momentum
		L = np.cross(r_vec, ME * v_vec)
		angular_momentum_z.append(L[2]) # Append z-component of angular momentum
		
		# Update position and velocity
		a_vec = AccelerationGravity(r_vec)
		v_vec += a_vec * dt
		r_vec += v_vec * dt
		
	return np.array(angular_momentum_z)
	
# Simulation parameters
total_time = 10 # Total time for one Earth year in years
time_steps = [0.1, 0.01, 0.001] # Different time steps

fig, axs = plt.subplots(1, len(time_steps), figsize=(18, 6))

for i, dt in enumerate(time_steps):
	angular_momentum_rk2 = SimulateOrbit_rk2(dt)
	angular_momentum_euler = SimulateOrbit_euler(dt)
	
	time = np.linspace(0, total_time, int(total_time / dt) + 1)
	
	axs[i].plot(time, angular_momentum_rk2, label='Runge-Kutta 2nd Order',color='blue')
	axs[i].plot(time, angular_momentum_euler, label='Euler Method', color='red', linestyle='--')
	axs[i].set_title('Time Evolution of Angular Momentum')
	axs[i].set_xlabel('Time (years)')
	axs[i].set_ylabel('Z-component of Angular Momentum (AU/year)')
	axs[i].legend()
	axs[i].grid(True)
	
plt.tight_layout()
plt.show()
