import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11 # gravitational constant in m^3/kg/s^2
M_sun = 1.989e30 # mass of the sun in kg
M_earth = 5.972e24 # mass of Earth in kg
M_jupiter = 1.898e27 # mass of Jupiter in kg
AU = 1 # 1 astronomical unit in meters
r_earth = 1.0 # initial distance of Earth from the sun in AU
r_jupiter = 5.2 # distance of Jupiter from the sun in AU

# Runge-Kutta 2nd order method
def RungeKutta(h, steps, theta_earth, theta_jupiter, r_earth, r_jupiter, v_earth, v_jupiter):
	earth_x = []
	earth_y = []
	jupiter_x = []
	jupiter_y = []
	
	for i in range(steps):
		# Calculating positions in AU
		x_earth = r_earth * np.cos(theta_earth)
		y_earth = r_earth * np.sin(theta_earth)
		x_jupiter = r_jupiter * np.cos(theta_jupiter)
		y_jupiter = r_jupiter * np.sin(theta_jupiter)
		
		earth_x.append(x_earth)
		earth_y.append(y_earth)
		jupiter_x.append(x_jupiter)
		jupiter_y.append(y_jupiter)
		
		# Calculating distances in AU
		rEJ = np.sqrt((x_earth - x_jupiter) ** 2 + (y_earth - y_jupiter) ** 2)
		
		# Calculating accelerations
		acc_earth_sun = G * M_sun / (r_earth ** 2)
		acc_jupiter_sun = G * M_sun / (r_jupiter ** 2)
		acc_earth_jupiter = G * M_jupiter / (rEJ ** 2)
		acc_jupiter_earth = acc_earth_jupiter
		
		# Total accelerations in x and y directions
		acc_x_earth = -acc_earth_sun * np.cos(theta_earth) + acc_earth_jupiter * np.cos(theta_jupiter)
		acc_y_earth = -acc_earth_sun * np.sin(theta_earth) + acc_earth_jupiter * np.sin(theta_jupiter)
		acc_x_jupiter = -acc_jupiter_sun * np.cos(theta_jupiter) + acc_jupiter_earth * np.cos(theta_earth)
		acc_y_jupiter = -acc_jupiter_sun * np.sin(theta_jupiter) + acc_jupiter_earth * np.sin(theta_earth)

		# Runge-Kutta method
		k1_theta_earth = h * v_earth
		k1_v_earth = h * acc_x_earth
		k1_theta_jupiter = h * v_jupiter
		k1_v_jupiter = h * acc_x_jupiter
			
		k2_theta_earth = h * (v_earth + 0.5 * k1_v_earth)
		k2_v_earth = h * (acc_x_earth + 0.5 * k1_v_earth)
		k2_theta_jupiter = h * (v_jupiter + 0.5 * k1_v_jupiter)
		k2_v_jupiter = h * (acc_x_jupiter + 0.5 * k1_v_jupiter)
		
		theta_earth += k2_theta_earth
		v_earth += k2_v_earth
		theta_jupiter += k2_theta_jupiter
		v_jupiter += k2_v_jupiter
		
	return earth_x, earth_y, jupiter_x, jupiter_y

# Run simulation
time_steps = [0.1, 0.01, 0.001]
	 
fig, axs = plt.subplots(1, len(time_steps), figsize=(16, 6))

for i, h in enumerate(time_steps):
	steps = int(10 / h) # Calculate the number of steps for each time step
	v_earth = np.sqrt(G * M_sun / (r_earth * AU)) / 1000
	v_jupiter = np.sqrt(G * M_sun / (r_jupiter * AU)) / 1000
	theta_earth = 0.0
	theta_jupiter = 2 * np.pi
	
	earth_x_vals, earth_y_vals, jupiter_x_vals, jupiter_y_vals = RungeKutta(h, steps, theta_earth, theta_jupiter, r_earth, r_jupiter, v_earth, v_jupiter)
	axs[i].plot(earth_x_vals, earth_y_vals, 'o', markersize=1, label=f'Earth h={h})')
	axs[i].plot(jupiter_x_vals, jupiter_y_vals, 'o', markersize=1, label=f'Jupiter (h={h})')
	
	axs[i].set_xlabel('X (AU)')
	axs[i].set_ylabel('Y (AU)')
	axs[i].set_title(f'Orbits of Earth and Jupiter (h={h})')
	axs[i].legend()
	axs[i].grid()
	axs[i].set_aspect('equal', adjustable='box')
	
plt.tight_layout()
plt.show()
