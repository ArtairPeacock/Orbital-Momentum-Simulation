import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.6743e-11 # Gravitational constant
M = 1.9891e30 # Mass of the sun
GM = 4 * (np.pi ** 2)
AU = 1 # 1 astronomical unit

# Initial conditions
r0 = AU # Initial distance of Earth to the sun
v0 = np.sqrt(GM / r0) # Initial velocity of Earth
t0 = 0 # Initial time
tf = 10 # Final time (in years)
dt = [0.1, 0.01, 0.001] # Time steps

# List to hold the values for the number of steps
numsteps = [int(tf / step) for step in dt]

# Arrays to hold the positions
x = np.zeros((len(dt), max(numsteps)))
y = np.zeros((len(dt), max(numsteps)))

# Arrays to hold the velocities
v_x = np.zeros((len(dt), max(numsteps)))
v_y = np.zeros((len(dt), max(numsteps)))

# Euler method
for i in range(len(dt)):
	x[i][0] = AU
	y[i][0] = 0
	v_x[i][0] = 0
	v_y[i][0] = v0

	for j in range(1, numsteps[i]):
		r = np.sqrt(x[i][j - 1] ** 2 + y[i][j - 1] ** 2)
		v_x[i][j] = v_x[i][j - 1] - (GM * x[i][j - 1] / r ** 3) * dt[i]
		v_y[i][j] = v_y[i][j - 1] - (GM * y[i][j - 1] / r ** 3) * dt[i]
		x[i][j] = x[i][j - 1] + v_x[i][j] * dt[i]
		y[i][j] = y[i][j - 1] + v_y[i][j] * dt[i]

	# Plotting
	plt.scatter(x[i][:numsteps[i]], y[i][:numsteps[i]], label=f"dt={dt[i]}",
alpha=0.7)

plt.title("Orbit of Earth using Euler's Method")
plt.xlabel("x (AU)")
plt.ylabel("y (AU)")
plt.legend()
plt.grid(True)
plt.show()
