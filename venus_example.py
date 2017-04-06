
# coding: utf-8

# In[1]:

import datetime
import math


# This is a direct copy of the Earth to Venus mission plan.  I'm doing this to make sure I get the functions correct, before proceeding further on the Earth to Mars.  

# Below is the capturing of the data for each planet.  I'm using a custom `PlanetaryObject` class to store the information

# In[2]:

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


# In[3]:

earth = PlanetaryObject(
    datetime.date(1988, 4, 8),
    197.53,     # Longitude
    0.01672,    # Eccentricity
    None,       # SMA
    None,       # Inclination
    102.29,     # Longitude of Perihelion
    0,          # Longitude of Ascending Node
    149.7848e6, # Radius
    29.75,      # Velocity
    95.24,      # True Anomaly
    0.9554,     # Flight Path Angle
    398600.4    # Gravitation parameter (km^3/s^2)
)


# In[4]:

venus = PlanetaryObject(
    datetime.date(1988, 7, 26),
    330.52,     # Longitude
    0.006778,   # Eccentricity
    None,       # SMA
    3.394,      # Inclination
    131.41,     # Longitude of Perihelion
    76.58,      # Longitude of Ascending Node
    108.9014e6, # Radius
    34.8,      # Velocity
    199.11,     # True Anomaly
    -0.128,     # Flight Path Angle
    324858.8    # Gravitation parameter (km^3/s^2)
)


# These are my formulas in python form.  They're based off of Table 3.3 found in the book

# In[55]:

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
    return (r_2 - r_1) / ((r_1 * math.cos(math.radians(theta_1))) - (r_2 * math.cos(math.radians(theta_2))))

def periapsis_radius(r, e, theta):
    """
    Calculates the periapsis radius of the transfer ellipse.  This is calculated 
    using the following equation:
    
    .. math::
        \frac {r_1 * [1 + e \cos{\theta]}} {1 + e}
    
    :param r: radius of the departing planetary object
    :param e: eccentricity of the transfer ellipse
    """
    return (r * (1 + e * math.cos(math.radians(theta)))) / (1 + e)

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
        cos_E = (e + math.cos(math.radians(theta))) / (1 + e * math.cos(math.radians(theta)))
        E = math.acos(cos_E)
        
    return (E - e * math.sin(E)) / n

def mean_motion(mu, a):
    """
    Calculates the mean motion of an elliptical orbit.  This is calculated 
    using the following equation:
    
    .. math::
        \sqrt{\frac{\mu} {a^3}}
    
    :param mu: gravitation parameter (Mass * Gravitation constant)
    :param a: semimajor axis
    """
    
    return math.sqrt(mu / a ** 3)

def velocity(mu, r, a):
    """
    Calculates the Velocity (V) of an object based on the elliptical orbit.  
    This is calculated using the following equation:
    
    .. math::
        \sqrt{\frac{2 * \mu} {r} - \frac{\mu} {a}}
        
    :param mu: gravitation parameter (Mass * Gravition constant)
    :param a: semimajor axis
    """
    return math.sqrt(2 * mu / r - mu / a)

def flight_path_angle(e, theta):
    """
    Calculates the Flight Path Angle (γ).  This is calculated using
    the following equation:
        
    .. math::
        \tan{γ} = {\frac{e * \sin{\theta}}{1 + 3 * \cos{\theta}}
        
    :param e: eccentricity of the elliptical orbit
    :param theta: 
    """
    tan_y = (e * math.sin(math.radians(theta))) / (1 + e * math.cos(math.radians(theta)))
    return math.atan(tan_y)

def inclination(Omega, L_s, L_t, i):
    a = math.radians(Omega + 180 - L_s)
    b = math.radians(L_t - (180 + Omega))
    alpha = math.radians(180 - i)
    cos_c = math.cos(a) * math.cos(b) + math.sin(a) * math.sin(b) * math.cos(alpha)
    c = math.acos(cos_c)
    sin_i_t = (math.sin(alpha) * math.sin(b)) / math.sin(c)
    return math.asin(sin_i_t)


# # Designing the Transfer Ellipse

# ## Time of Flight

# In[6]:

venus.date - earth.date


# In[7]:

time_of_flight = venus.date - earth.date
time_of_flight = time_of_flight.days
time_of_flight


# ## Eccentricity

# In[8]:

line_of_apisides = 180
true_anom = line_of_apisides + (venus.L - earth.L)
true_anom


# In[9]:

eccentricity(earth.r, venus.r, line_of_apisides, true_anom)


# In[10]:

e = eccentricity(earth.r, venus.r, line_of_apisides, true_anom)


# ## Periapsis Radius

# In[11]:

periapsis_radius(earth.r, e, line_of_apisides)


# In[12]:

r_p = periapsis_radius(earth.r, e, line_of_apisides)


# ## Semi-Major Axis

# In[13]:

# Book apparently rounds the actual values here
semimajor_axis(r_p=103.555e6, e=0.1825)


# In[14]:

a = 126.673e6


# ## Time of Flight

# In[15]:

n = mean_motion(mu_sun, a)
n


# In[16]:

peri_to_earth = time_since_periapsis(e, n, theta=line_of_apisides)
peri_to_earth / 3600 / 24 # conversion from seconds to days


# In[17]:

venus_to_peri = time_since_periapsis(e, n, theta=true_anom)
venus_to_peri / 3600 / 24


# In[18]:

(peri_to_earth - venus_to_peri) / 3600 / 24


# ## Velocities

# In[19]:

velocity(mu_sun, earth.r, 129.336e6)  # using the Value from the Book which appear to be rounded


# In[20]:

velocity(mu_sun, venus.r, 129.336e6)  # again using the values from the book which appear to be rounded


# ## Flight Path Angles

# In[21]:

math.degrees(flight_path_angle(0.17194, 199.53))  # same as above, using the book values 


# In[22]:

math.degrees(flight_path_angle(0.17194, 332.52))


# Now that I've verified the fundamental functions above, let's wrap this all up into a nice function that'll optimize this for us

# In[23]:

def transfer_ellipse(start_planet, end_planet, tof_accuracy=2, max_iters=1000, return_trials=False):
    time_of_flight = end_planet.date - start_planet.date
    time_of_flight = time_of_flight.days
    
    longs = []
    tofs = []
    
    line_of_apisides = 180  # trial start
    tof = 9999999999 # large number to get us started
    
    bottom_angle = 90
    top_angle = 270
    
    
    i = 0
    while not(time_of_flight - 10e-tof_accuracy < tof / 3600 / 24 < time_of_flight + 10e-tof_accuracy) and i < max_iters:
        
        line_of_apisides = (bottom_angle - top_angle) / 2
        
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
        
        if tof / 3600 / 24 > time_of_flight:
            
        
        
        i += 1
    
#     while  tof / 3600 / 24 > time_of_flight:
#         true_anom = line_of_apisides + (end_planet.L - start_planet.L)
#         longs.append((line_of_apisides, true_anom))
#         e = eccentricity(start_planet.r, end_planet.r, line_of_apisides, true_anom)
#         r_p = periapsis_radius(start_planet.r, e, line_of_apisides)
#         a = semimajor_axis(r_p=r_p, e=e)
        
#         n = mean_motion(mu_sun, a)
        
#         peri_to_start = time_since_periapsis(e, n, theta=line_of_apisides)
#         end_to_peri = time_since_periapsis(e, n, theta=true_anom)
        
#         tof = peri_to_start - end_to_peri
#         tofs.append(tof / 3600 / 24)
#         line_of_apisides += 1
        
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
        'y_start': math.degrees(y_start),
        'y_end': math.degrees(y_end)
    }
        
        
    if return_trials:
        r_dict.update({'runs':{'longs': longs, 'tofs':tofs}})
        
    return r_dict


# In[54]:

tf = transfer_ellipse(earth, venus, return_trials=True)
tf


# Interestingly enough, we're getting $\theta_{Earth} = 194$, however the book claims that $\theta_{Earth} = 199$.  I believe the discrepency here is found with the fact that the book apparently rounds their vaules while the values used by the functions above are more accurate.  

# # Designing the Departure Trajectory

# ## Plane Change

# In[25]:

alpha = 180 - venus.i
alpha


# In[26]:

a = venus.asc + 180 - earth.L
a


# In[27]:

b_prime = venus.L - (venus.asc + 180)
b_prime


# In[28]:

# b = b_prime # this can be used when the transfer angles are small.
b = 73.967 # taken from the book b/c after much research, I still don't know how to solve a spherical right triangle


# In[29]:

csc_c = math.cos(math.radians(a)) * math.cos(math.radians(b)) + math.sin(math.radians(a)) * math.sin(math.radians(b)) * math.cos(math.radians(alpha))
csc_c


# In[30]:

c = math.degrees(math.acos(csc_c))
c


# In[31]:

sin_i = (math.sin(math.radians(alpha)) * math.sin(math.radians(b))) / math.sin(math.radians(c))
sin_i


# In[32]:

i_t = math.degrees(math.asin(sin_i))
i_t


# ## Calculating $V_{HE}$ and C3

# In[33]:

# cos_alpha = math.cos(math.radians(i_t)) * math.cos(math.radians(earth.fp + tf['y_start']))
cos_alpha = math.cos(math.radians(4.455)) * math.cos(math.radians(earth.fp + 3.924))  # using the value from the book, since my is different (and more accurate I believe)
cos_alpha


# In[34]:

alpha = math.degrees(math.acos(cos_alpha))
alpha


# In[35]:

#C3 = earth.v ** 2 + tf['V_start'] ** 2 - 2 * earth.v * tf['V_start'] * math.cos(math.radians(alpha))
C3 = earth.v ** 2 + 27.312 ** 2 - 2 * earth.v * 27.312 * math.cos(math.radians(alpha))
C3


# In[36]:

V_he = math.sqrt(C3)
V_he


# Similar to what we did for the Transfer Ellipse, let's combine all these steps into a single function to calculate these for us:

# In[37]:

def depart_trajectory(start_planet, end_planet, y, V):
    alpha = 180 - end_planet.i
    a = end_planet.asc + 180 - start_planet.L
    b = end_planet.L - (end_planet.asc + 180)
    csc_c = math.cos(math.radians(a)) * math.cos(math.radians(b)) + math.sin(math.radians(a)) * math.sin(math.radians(b)) * math.cos(math.radians(alpha))
    c = math.degrees(math.acos(csc_c))
    sin_i = (math.sin(math.radians(alpha)) * math.sin(math.radians(b))) / math.sin(math.radians(c))
    i_t = math.degrees(math.asin(sin_i))
    
    # if they have the same sign, subtract them, else add them
    if start_planet.fp * y > 0:
        y_s = abs(start_planet.fp) - abs(y)
    else:
        y_s = abs(start_planet.fp) + abs(y)
    
    cos_alpha = math.cos(math.radians(i_t)) * math.cos(math.radians(y_s))
    alpha = math.degrees(math.acos(cos_alpha))
    C3 = start_planet.v ** 2 + V ** 2 - 2 * start_planet.v * V * math.cos(math.radians(alpha))
    V_he = math.sqrt(C3)
    
    r_dict = {
        'i_t': i_t,
        'C3': C3,
        'V_he': V_he
    }
    
    return r_dict


# In[38]:

depart_trajectory(earth, venus, -3.924, 27.312)


# # Designing the Arrival Trajectory

# ## Plane Change

# In[39]:

alpha = 180 - venus.i
alpha


# In[40]:

a = venus.asc + 180 - earth.L
a


# In[41]:

b_prime = venus.L - (venus.asc + 180)
b_prime


# In[42]:

# b = b_prime # this can be used when the transfer angles are small.
b = 73.967 # taken from the book b/c after much research, I still don't know how to solve a spherical right triangle
b


# In[43]:

csc_c = math.cos(math.radians(a)) * math.cos(math.radians(b)) + math.sin(math.radians(a)) * math.sin(math.radians(b)) * math.cos(math.radians(alpha))
csc_c


# In[44]:

c = math.degrees(math.acos(csc_c))
c


# In[45]:

sin_it = math.sin(math.radians(alpha)) * math.sin(math.radians(a)) / math.sin(math.radians(c))
sin_it


# In[46]:

it = math.degrees(math.asin(sin_it))
it


# ## Calculating $V_\infty$

# In[47]:

#cos_alpha_inf = math.cos(math.radians(it)) * math.cos(math.radians(tf['y_end'] + venus.fp))
cos_alpha_inf = math.cos(math.radians(it)) * math.cos(math.radians(3.938 + venus.fp))
alpha_inf = math.acos(cos_alpha_inf)
math.degrees(alpha_inf)


# In[48]:

#C3 = venus.v ** 2 + tf['V_end'] ** 2 + 2 * venus.v * tf['V_end'] * math.cos(alpha_inf)
C3 = venus.v ** 2 + 37.57 ** 2 - 2 * venus.v * 37.57 * math.cos(math.radians(5.5039))
V_inf = math.sqrt(C3)
V_inf# should be 4.442 km/s


# In[49]:

def arrival_trajectory(start_planet, end_planet, y, V):
    alpha = 180 - end_planet.i
    a = end_planet.asc + 180 - start_planet.L
    b = end_planet.L - (end_planet.asc + 180)
    csc_c = math.cos(math.radians(a)) * math.cos(math.radians(b)) + math.sin(math.radians(a)) * math.sin(math.radians(b)) * math.cos(math.radians(alpha))
    c = math.degrees(math.acos(csc_c))
    sin_it = math.sin(math.radians(alpha)) * math.sin(math.radians(a)) / math.sin(math.radians(c))
    it = math.degrees(math.asin(sin_it))
    
    # if they have the same sign, subtract them, else add them
    if end_planet.fp * y > 0:
        y_s = abs(abs(end_planet.fp) - abs(y))
    else:
        y_s = abs(abs(end_planet.fp) + abs(y))
    
    cos_alpha_inf = math.cos(math.radians(it)) * math.cos(math.radians(y_s + end_planet.fp))
    alpha_inf = math.acos(cos_alpha_inf)
    C3 = end_planet.v ** 2 + V ** 2 - 2 * end_planet.v * V * math.cos(math.radians(alpha_inf))
    V_inf = math.sqrt(C3)
    
    r_dict = {
        'i_t': it,
        'V_inf': V_inf
    }
    
    return r_dict


# In[50]:

arrival_trajectory(earth, venus, -3.938, 37.57)


# We're getting different answers here, becuase our angles are a little different.  `alpha_inf` as calculated by the book is 5.5039 while I'm getting 5.5036.  This is due to the rounding of the $i_{tp}$ as found in the book.  I'm getting 3.9745967799374893 while the books rounds this to 3.975.  See calculation below:

# In[51]:

math.degrees(math.acos(math.cos(math.radians(3.975)) * math.cos(math.radians(3.938-0.128))))


# ## Combining the Trajectories into a single function:

# In[52]:

def trajectories(start_planet, end_planet, y_start, y_end, V_start, V_end):
    alpha = 180 - end_planet.i
    a = end_planet.asc + 180 - start_planet.L
    b = end_planet.L - (end_planet.asc + 180)
    csc_c = math.cos(math.radians(a)) * math.cos(math.radians(b)) + math.sin(math.radians(a)) * math.sin(math.radians(b)) * math.cos(math.radians(alpha))
    c = math.degrees(math.acos(csc_c))
    
    sin_i_start = (math.sin(math.radians(alpha)) * math.sin(math.radians(b))) / math.sin(math.radians(c))
    i_start = math.degrees(math.asin(sin_i_start))
    
    sin_i_end = math.sin(math.radians(alpha)) * math.sin(math.radians(a)) / math.sin(math.radians(c))
    i_end = math.degrees(math.asin(sin_i_end))
    
    # if they have the same sign, subtract them, else add them
    if start_planet.fp * y_start > 0:
        y_s = abs(abs(start_planet.fp) - abs(y_start))
    else:
        y_s = abs(abs(start_planet.fp) + abs(y_start))
    
    cos_alpha = math.cos(math.radians(i_start)) * math.cos(math.radians(y_s))
    alpha = math.degrees(math.acos(cos_alpha))
    C3 = start_planet.v ** 2 + V_start ** 2 - 2 * start_planet.v * V_start * math.cos(math.radians(alpha))
    V_he = math.sqrt(C3)
    
    if end_planet.fp * y_end > 0:
        y_e = abs(abs(end_planet.fp) - abs(y_end))
    else:
        y_e = abs(abs(end_planet.fp) + abs(y_end))
        
    cos_alpha_inf = math.cos(math.radians(i_end)) * math.cos(math.radians(y_s + end_planet.fp))
    alpha_inf = math.acos(cos_alpha_inf)
    C3_inf = end_planet.v ** 2 + V_end ** 2 - 2 * end_planet.v * V_end * math.cos(math.radians(alpha_inf))
    V_inf = math.sqrt(C3_inf)
    
    r_dict = {
        'i_start': i_start,
        'C3': C3,
        'V_he': V_he,
        'i_end': i_end,
        'V_inf': V_inf
    }
    
    return r_dict


# In[53]:

trajectories(earth, venus, -3.924, -3.938, 27.312, 37.57)


# In[ ]:



