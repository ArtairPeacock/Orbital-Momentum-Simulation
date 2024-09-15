import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11 # Gravitational constant
MS = 1.9891e30 # Mass of the Sun in kg
ME = 5.972e24 # Mass of the Sun in kg
AU = 1 # 1 Astronomical Unit
r = AU # Initial distance of Earth to the Sun in AU (1 AU)
v_K = np.sqrt(G * MS / r) # Keplerian velocity in m/s

# Time parameters
total_time_years = 10 # Total time for 10 Earth years
time_steps = [0.1, 0.01, 0.001] # Different time steps in years

# Function to calculate acceleration
def AccelerationGravity(r_vec):
	r_mag = np.linalg.norm(r_vec)
	return (-G * MS / r_mag**3) * r_vec
	
# Function implementing Euler method for orbit simulation
def SimulateOrbit(dt):
	num_steps = int(total_time_years / dt)
	time = np.linspace(0, total_time_years, num_steps + 1)
	
	r_vec = np.array([r, 0, 0], dtype=float) # Initial position vector
	v_vec = np.array([0, v_K, 0], dtype=float) # Initial velocity vector
	
	angular_momentum_z = []
	
	for t in time:
		a_vec = AccelerationGravity(r_vec)
		v_vec += a_vec * dt
		r_vec += v_vec * dt
		
		# Calculate linear momentum
		p = ME * v_vec
		
		# Calculate angular momentum L = r x p
		L = np.cross(r_vec, p)
		angular_momentum_z.append(L[2]) # Append z-component of angular momentum
		
	return time, angular_momentum_z
	
# Plotting the z-component of angular momentum for different time steps
plt.figure(figsize=(10, 6))
for dt in time_steps:
	time, angular_momentum_z = SimulateOrbit(dt)
	plt.plot(time, angular_momentum_z, label=f"Time Step: {dt} years")

plt.title('Time Evolution of z-component of Angular Momentum of Earth')
plt.xlabel('Time (years)')
plt.ylabel('Z-component of Angular Momentum (AU/year)')
plt.legend()
plt.grid(True)
plt.show()
