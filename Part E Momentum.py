import numpy as np
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11 # gravitational constant in m^3/kg/s^2
M_sun = 1.989e30 # mass of the sun in kg
M_earth = 5.972e24 # mass of Earth in kg
M_jupiter = 1.898e27 # mass of Jupiter in kg
year = 365 * 24 * 3600 # 1 year in seconds
r_earth = 1.0 # initial distance of Earth from the sun in AU
r_jupiter = 5.2 # distance of Jupiter from the sun in AU

# Initial conditions
v_earth = np.sqrt(G * M_sun / r_earth) # Keplerian velocity for Earth in m/s
v_jupiter = np.sqrt(G * M_sun / r_jupiter) # Keplerian velocity for Jupiter in m/s
theta_earth = 0.0 # initial angle for Earth
theta_jupiter = np.pi # initial angle for Jupiter
t_max = 10 # maximum time in years

def RungeKutta(h, t_max, theta_earth, theta_jupiter, v_earth, v_jupiter):
	L_earth = [] # Angular momentum for Earth
	L_jupiter = [] # Angular momentum for Jupiter
	time_steps = [] # Time steps
	
	steps = int(t_max / h) # number of steps
	
	for i in range(steps):
		# Calculating Earth's and Jupiter's positions in AU
		x_earth = r_earth * np.cos(theta_earth)
		y_earth = r_earth * np.sin(theta_earth)
		x_jupiter = r_jupiter * np.cos(theta_jupiter)
		y_jupiter = r_jupiter * np.sin(theta_jupiter)
		
		# Calculating distances between planets
		r_ej = np.sqrt((x_earth - x_jupiter)**2 + (y_earth - y_jupiter)**2)
		
		# Calculating accelerations
		a_earth = G * M_jupiter / (r_ej**3) * r_earth
		a_jupiter = G * M_earth / (r_ej**3) * r_jupiter
		
		# Update velocities
		v_earth += h * a_earth
		v_jupiter += h * a_jupiter
		
		# Update angles using velocities and time step
		theta_earth += (h * v_earth) / r_earth
		theta_jupiter += (h * v_jupiter) / r_jupiter
		
		# Calculating angular momentum for Earth and Jupiter
		L_earth.append(M_earth * r_earth * v_earth)
		L_jupiter.append(M_jupiter * r_jupiter * v_jupiter)
		
		time_steps.append(i * h)
		
	return time_steps, L_earth, L_jupiter
	
# Different time steps
time_step_options = [0.1, 0.01, 0.001] # Define different time steps
fig, axs = plt.subplots(2, figsize=(10, 10))

for h in time_step_options:
	time, L_earth_vals, L_jupiter_vals = RungeKutta(h, t_max, theta_earth, theta_jupiter, v_earth, v_jupiter)
	axs[0].plot(np.array(time) / year, L_earth_vals, label=f'Earth (h={h} years)')
	axs[1].plot(np.array(time) / year, L_jupiter_vals, label=f'Jupiter (h={h} years)')
	
axs[0].set_xlabel('Time (years)')
axs[0].set_ylabel('Angular Momentum (AU/year)')
axs[0].set_title('Earth - Z-component of Angular Momentum')
axs[0].legend()
axs[0].grid()
axs[1].set_xlabel('Time (years)')
axs[1].set_ylabel('Angular Momentum (AU/year)')
axs[1].set_title('Jupiter - Z-component of Angular Momentum')
axs[1].legend()
axs[1].grid()

plt.tight_layout()
plt.show()
