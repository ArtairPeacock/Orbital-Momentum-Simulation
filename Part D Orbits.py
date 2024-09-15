import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11 # gravitational constant
MS = 1.989e30 # mass of the Sun in kg
ME = 5.972e24 # mass of the Earth in kg

def f(t, r):
	# r[0] = x, r[1] = y, r[2] = vx, r[3] = vy
	dx = r[2]
	dy = r[3]
	r_mag = np.sqrt(r[0]**2 + r[1]**2)
	F = -G * MS * ME / r_mag**3
	dvx = F * r[0] / ME
	dvy = F * r[1] / ME
	return np.array([dx, dy, dvx, dvy])

def RungeKuttaSecondOrder(t, r, h):
	k1 = h * f(t, r)
	k2 = h * f(t + h, r + k1)
	return r + 0.5 * (k1 + k2)
	
# Initial conditions
t = 0.0
dt_values = [0.1, 0.01, 0.001] # time steps in years

# Lists to store the orbit positions for each time step
orbit_positions = []

# Simulation for each time step
for dt in dt_values:
	t = 0.0
	r = np.array([1.4710e11, 0.0, 0.0, 3.0287e4]) # initial position and velocity
	of Earth [x, y, vx, vy]
	x_positions = []
	y_positions = []
		
	num_steps = int(365 / dt) # simulate for 1 year with the given time step
	
	for i in range(num_steps):
		x_positions.append(r[0] / 1.496e11) # convert x position to AU
		y_positions.append(r[1] / 1.496e11) # convert y position to AU
		r = RungeKuttaSecondOrder(t, r, dt * 24 * 3600) # convert dt from years to
		seconds
		t += dt
		orbit_positions.append((x_positions, y_positions))
		
# Plotting the orbits for different time steps
plt.figure(figsize=(8, 8))
for i, (x_positions, y_positions) in enumerate(orbit_positions):
	plt.plot(x_positions, y_positions, label=f'Time Step: {dt_values[i]} years')
	
plt.scatter(0, 0, color='yellow', marker='o', label='Sun')
plt.title('Orbit of Earth around the Sun')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()
