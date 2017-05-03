import numpy as np

class PlanetaryObject():
    """
    A simple class used to store pertinant information about the plantary object
    """
    def __init__(self, date, L, e, SMA, i, peri, asc, r, v, anom, fp, mu):
        self.date = date   # Event Date
        self.L = L         # Longitude
        self.e = e         # Eccentricity
        self.SMA = SMA     # SMA
        self.i = i         # Inclination
        self.peri = peri   # Longitude of Perihelion
        self.asc = asc     # Longitude of Ascending Node
        self.r = r         # Radius
        self.v = v         # Velocity
        self.anom = anom   # True Anomaly
        self.fp = fp       # Flight Path Angle
        self.mu = mu       # Gravitation parameter


mu_sun = 132712439935.5
        
def eccentricity(r_1, r_2, theta_1, theta_2):
    """
    Calculates the eccentricity of the transfer ellipse.  This is calculated through 
    the following equation:
    
    .. math::
        \frac {r_2 - r_1} {r_1 * \cos{\theta_1} - r_2 * \cos{\theta_2}}
    
    :param r_1: radius of the departing planetary object
    :param r_2: radius of the arriving planetary object
    :param theta_1: True anomaly of the departing planetary object in degrees
    :param theta_2: True anomaly of the arriving planetary object in degrees
    """
    return (r_2 - r_1) / ((r_1 * np.cos(np.radians(theta_1))) - (r_2 * np.cos(np.radians(theta_2))))

def periapsis_radius(r, e, theta):
    """
    Calculates the periapsis radius of the transfer ellipse.  This is calculated 
    using the following equation:
    
    .. math::
        \frac {r_1 * [1 + e \cos{\theta]}} {1 + e}
    
    :param r: radius of the departing planetary object
    :param e: eccentricity of the transfer ellipse
    """
    return (r * (1 + e * np.cos(np.radians(theta)))) / (1 + e)

def semimajor_axis(r=None, r_a=None, r_p=None, mu=None, V=None, e=None):
    """
    Calculates the semi-major axis of the transfer ellipse.  This is calculated 
    using one of the following equations:
    
    .. math::
        \frac {r_a + r_p} {2}
        
        \frac {\mu r} {2 \mu - V^2 r}
        
        \frac {r_p} {1 - e}
        
        \frac {r_a} {1 + e}
    
    :param r: general radius of the elliptical orbit
    :param r_a: Radius of apoapsis
    :param r_p: Radius of periapsis
    :param mu: gravitation parameter
    :param V: Velocity of the orbiting object
    :param e: Eccentricity of the elliptical orbit
    """
    if r_a != None and r_p != None:
        return (r_a + r_p) / 2
    if mu != None and r !=None and V != None:
        return (mu * r) / (2 * mu - V ** 2 * r)
    if r_p != None and e != None:
        return r_p / (1 - e)
    if r_a != None and e != None:
        return r_a / (1 + e)
    
    # If we reach this point, then the passed in arguments doesn't match
    #    any equations we have defined.  Raise an Error
    raise TypeError("Invalid arguments!")
    

def time_since_periapsis(e, n, theta=None, E=None):
    """
    Calculates the time since the periapsis.  This is calculated using the
    following equation:
    
    .. math::
        \frac {E - e \sin{E}} {n}
        
    If E, isn't defined, it will be calculated using the param theta and 
    the following equation:
    
    ..math:: 
        
        \cos {E} = \frac {e + \cos{\theta}} {1 + e \cos{\theta}}
    
    :param e: eccentricity of the transfer ellipse
    :param n: mean motion
    :param theta: degrees to periapsis
    :param E: eccentric anomaly in radians
    """
    if theta == None and E == None:
        raise TypeError("theta or E MUST be defined")
    if theta != None and E != None:
        raise TypeError("theta OR E must be defined.  Not both")
        
    if E == None:
        cos_E = (e + np.cos(np.radians(theta))) / (1 + e * np.cos(np.radians(theta)))
        E = np.arccos(cos_E)
        
    return (E - e * np.sin(E)) / n

def mean_motion(mu, a):
    """
    Calculates the mean motion of an elliptical orbit.  This is calculated 
    using the following equation:
    
    .. math::
        \sqrt{\frac{\mu} {a^3}}
    
    :param mu: gravitation parameter (Mass * Gravitation constant)
    :param a: semimajor axis
    """
    
    return np.sqrt(mu / a ** 3)

def velocity(mu, r, a):
    """
    Calculates the Velocity (V) of an object based on the elliptical orbit.  
    This is calculated using the following equation:
    
    .. math::
        \sqrt{\frac{2 * \mu} {r} - \frac{\mu} {a}}
        
    :param mu: gravitation parameter (Mass * Gravition constant)
    :param a: semimajor axis
    """
    return np.sqrt(2 * mu / r - mu / a)

def flight_path_angle(e, theta):
    """
    Calculates the Flight Path Angle (γ).  This is calculated using
    the following equation:
        
    .. math::
        \tan{γ} = {\frac{e * \sin{\theta}}{1 + 3 * \cos{\theta}}
        
    :param e: eccentricity of the elliptical orbit
    :param theta: 
    """
    tan_y = (e * np.sin(np.radians(theta))) / (1 + e * np.cos(np.radians(theta)))
    return np.arctan(tan_y)

def inclination(Omega, L_s, L_t, i):
    a = np.radians(Omega + 180 - L_s)
    b = np.radians(L_t - (180 + Omega))
    alpha = np.radians(180 - i)
    cos_c = np.cos(a) * np.cos(b) + np.sin(a) * np.sin(b) * np.cos(alpha)
    c = np.arccos(cos_c)
    sin_i_t = (np.sin(alpha) * np.sin(b)) / np.sin(c)
    return np.arcsin(sin_i_t)


# We'll split this problem up into 3 different sections:
# 
# 1. Transfer Ellipse
# 2. Departure Trajectory
# 3. Arrival Trajectory
# 
# This approach is called the "patched conics approach".  

# ## Designing the Transfer Ellipse
# 
# In order to create the transfer ellipse between Earth and Mars, we're going to do the following:
# 
# 1. Determine the required # of days based off departure and arrival date
# 2. Calculate the angular difference between Earth and Mars
# 3. Design the Eccentricity, Periapsis Axis and Semi-Major Axis using an iterative approach
# 4. Finally, calculating the respective velocities and flight path angles for each planetary object
# 
# To find the required number of days to travel to Mars, we simply subtract the launch and arrival dates:

# In[4]:

def transfer_ellipse(start_planet, end_planet, return_trials=False):
    time_of_flight = end_planet.date - start_planet.date
    time_of_flight = time_of_flight.days
    
    longs = []
    tofs = []
    
    line_of_apisides = 180  # trial start
    tof = 9999999999 # large number to get us started
    
    while  tof / 3600 / 24 > time_of_flight:
        true_anom = line_of_apisides + (end_planet.L - start_planet.L)
        longs.append((line_of_apisides, true_anom))
        e = eccentricity(start_planet.r, end_planet.r, line_of_apisides, true_anom)
        r_p = periapsis_radius(start_planet.r, e, line_of_apisides)
        a = semimajor_axis(r_p=r_p, e=e)
        
        n = mean_motion(mu_sun, a)
        
        peri_to_start = time_since_periapsis(e, n, theta=line_of_apisides)
        end_to_peri = time_since_periapsis(e, n, theta=true_anom)
        
        tof = peri_to_start - end_to_peri
        tofs.append(tof / 3600 / 24)
        line_of_apisides += 1
        
    # Calculate the Relative Velocities
    V_start = velocity(mu_sun, start_planet.r, a)
    V_end = velocity(mu_sun, end_planet.r, a)
    
    y_start = flight_path_angle(e, line_of_apisides)
    y_end = flight_path_angle(e, true_anom)
    
    r_dict = {
        'line_of_apisides': line_of_apisides - 1,  # subtract the 1 we added during the loop
        'true_anom': true_anom,
        'eccentricity': e,
        'SMA': a,
        'time_of_flight': tof,
        'V_start': V_start,
        'V_end': V_end,
        'y_start': np.degrees(y_start),
        'y_end': np.degrees(y_end)
    }
        
        
    if return_trials:
        r_dict.update({'runs':{'longs': longs, 'tofs':tofs}})
        
    return r_dict


# In[5]:

time_of_flight = mars.date - earth.date
time_of_flight = time_of_flight.days

longs = []
tofs = []

line_of_apisides = 180  # trial start
tof = 9999999999 # large number to get us started

bottom_angle = 90
top_angle = 270

for line_of_apisides in range(180, 271):
    true_anom = line_of_apisides + (mars.L - earth.L)
    longs.append((line_of_apisides, true_anom))
    e = eccentricity(earth.r, mars.r, line_of_apisides, true_anom)
    r_p = periapsis_radius(earth.r, e, line_of_apisides)
    a = semimajor_axis(r_p=r_p, e=e)

    n = mean_motion(mu_sun, a)

    peri_to_start = time_since_periapsis(e, n, theta=line_of_apisides)
    end_to_peri = time_since_periapsis(e, n, theta=true_anom)

    tof = peri_to_start - end_to_peri
    tofs.append(tof / 3600 / 24)


# In[ ]:




# In[ ]:

tf = transfer_ellipse(earth, mars, return_trials=True)
tf


# ## The Departure Trajectory
# 
# The departure trajectory entails the following:
# 
# * Determining the plane change
# * Calculating $V_{HE}$ and $C3$
# 
# $V_{HE}$ is the vector required for the spacecraft to leave the Earth's gravitational pull, and enter into the calculated transfer ellipse above.

# We'll use the function defined above to find the inclination of the transfer plane $i_t$:

# In[ ]:

i_t = inclination(mars.asc, earth.L, mars.L, mars.i)
print("i_t = {:.2f}°".format(np.degrees(i_t)))


# Using $i_t$, we can now determine the $V_{HE}$ and $C3$ of the departing trajectory:

# In[ ]:

cos_alpha_2 = np.cos(i_t) * np.cos(earth.fp + abs(y_ne))
alpha_2 = np.arccos(cos_alpha_2)
C3 = earth.v ** 2 + V_ne ** 2 - 2 * earth.v * V_ne * np.cos(alpha_2)
V_he = np.sqrt(C3)
print("C3 = {:.2f} km^2/s^2; V_he = {:.2f} km/s".format(C3, V_he))


# ## The Arrival Trajectory
# 
# Similar to what was done for the departure trajectory, we'll calculate the plane change and then determine $V_\infty$.  

# In[ ]:



