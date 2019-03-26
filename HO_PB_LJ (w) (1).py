import numpy as np
import numpy.random as npr

# distance vector between two particles
def distance_vec(pos_vec1, pos_vec2):
    d_vec = pos_vec1 - pos_vec2
    
    if boundarycond == True:
        if d_vec[0] < -L_x/2:
            d_vec[0] = d_vec[0] + L_x
        elif d_vec[0] > L_x/2:
            d_vec[0] = d_vec[0] - L_x
        
        if d_vec[1] < -L_y/2:
            d_vec[1] = d_vec[1] + L_y
        elif d_vec[1] > L_y/2:
            d_vec[1] = d_vec[1] - L_y
            
    return d_vec

# norm of the above distance vector
def distance(dist_vec):
    return np.linalg.norm(dist_vec)

# Force between two Lennard-Jones particles
def lennardjonesforce(dist_vec, timestep):
    g = distance(dist_vec)**2
    calculateU_lj(dist_vec, timestep)
    return 24*epsilon*(2*sigma**12/g**7-sigma**6/g**4)*dist_vec  # BE CAREFUL OF THE SIGN HERE! IT MUST BE +!

def calculateU_lj(dist_vec, timestep):
    energy_lj = 4 * epsilon * ((sigma / distance(dist_vec))**12 - ((sigma / distance(dist_vec))**6)) - 4 * epsilon * ((sigma / r_cutoff)**12 - (sigma / r_cutoff)**6)
    U_lj[timestep] = U_lj[timestep] + energy_lj

def calculateU_h(dist, timestep):
    energy_h = 0.5 * k * dist**2
    U_h[timestep] = U_h[timestep] + energy_h
                      
def calculateK(vel_vec1, vel_vec2, timestep):
    energy_k = 0.5 * m[0] * np.linalg.norm(vel_vec1)**2 + 0.5 * m[1] * np.linalg.norm(vel_vec2)**2       
    K[timestep] = K[timestep] + energy_k
        
# script for writing the .xyz file for visualization in VMD
def write_xyz(data, filename) -> None:
    with open(filename, 'w') as f:
        for j in range(data.shape[1]):
            print(f'{data.shape[0]}\n', file=f)
            for i in range(data.shape[0]):
                print(f'p{i}   {data[i, j, 0]}   {data[i, j, 1]}   0', file=f)

########################################### Simulation parameters #####################################################

N = 2                                        # Number of particles in the HO
M = 3                                        # Number of HOs
num_step = 10000                             # Number of timesteps
dt = 0.0001                                   # Timestep size (accuracy)
dim = 2                                      # Number of dimensions

boundarycond = False                          # Set boundary conditions on (True) or off (false)
harmonicF = True                             # Set harmonic force on (True) or off (false)
lennardjonF = True                           # Set Lennard-Jones force on (True) or off (false)

sigma = 1
epsilon = 1
r_cutoff = 2.5*sigma

######################################### Periodic boundary conditions ################################################

min_x = -10
max_x = 10
min_y = -10
max_y = 10

L_x = max_x - min_x
L_y = max_y - min_y

############################################ Spring properties ########################################################

bond_length = 5                              # Bond length (equilibrial position)
k = 30                                       # Spring constant

########################################### Particle properties #######################################################

m = np.array([1,1])                          # Particle masses

####################################### Randomized initial conditions #################################################

x1 = npr.uniform(-10, 10,(M, dim))              # Particle 1, initial position
x2 = npr.uniform(-10, 10,(M, dim))              # Particle 2, initial position
x1_timestep = x1
x2_timestep = x2

v1 = npr.uniform(-10,10,(M, dim))              # Particle 1, initial velocity
v2 = npr.uniform(-10,10,(M, dim))              # Particle 2, initial velocity

U_h = np.zeros((num_step))
K = np.zeros((num_step))
U_lj = np.zeros((num_step))
E = np.zeros((num_step))

####################################### Determined initial conditions #################################################

#x1 = np.array([-5,3])
#x2 = np.array([3,-1])

#v1 = np.array([3,2])
#v2 = np.array([-1,5])

particles = np.zeros((N*M, num_step, dim))

for j in range(num_step):
    
    for i in range(M):

        # Leapfrog steps 1 and 2
        F1 = 0
        F2 = 0

        if harmonicF == True:
            d_vec = distance_vec(x1_timestep[i, :], x2_timestep[i, :])        # not normalized
            d = distance(d_vec)
            F1 = F1 - k * d_vec * (1 - (bond_length / d))
            F2 = F2 + k * d_vec * (1 - (bond_length / d))
            calculateU_h(d - bond_length, j)

        if lennardjonF == True:
            
            for p in range(M):

                d_vec1 = distance_vec(x1_timestep[i, :], x2_timestep[p, :]) 
                F1 = F1 + lennardjonesforce(d_vec1, j)

                if i != p:
                    d_vec2 = distance_vec(x1_timestep[i, :], x1_timestep[p, :]) 
                    F1 = F1 + lennardjonesforce(d_vec2, j)

                d_vec3 = distance_vec(x2_timestep[i, :], x1_timestep[p, :]) 
                F2 = F2 + lennardjonesforce(d_vec3, j)

                if i != p:
                    d_vec4 = distance_vec(x2_timestep[i, :], x2_timestep[p, :]) 
                    F2 = F2 + lennardjonesforce(d_vec4, j)

        # Leapfrog step 3
        v1[i, :] = v1[i, :] + (dt / m[0]) * F1
        v2[i, :] = v2[i, :] + (dt / m[1]) * F2
        
        calculateK(v1[i, :], v2[i, :], j)

        # Leapfrog step 4
        x1[i, :] = x1[i, :] + v1[i, :] * dt
        x2[i, :] = x2[i, :] + v2[i, :] * dt
        
        if boundarycond == True:
          
            #Periodic boundaries in x-axis
            if x1[i, 0] < min_x:
                x1[i, 0] = x1[i, 0] + L_x
            elif x1[i, 0] > max_x:
                x1[i, 0] = x1[i, 0] - L_x
            
            if x2[i, 0] < min_x:
                x2[i, 0] = x2[i, 0] + L_x
            elif x2[i, 0] > max_x:
                x2[i, 0] = x2[i, 0] - L_x
            
            #Periodic boundaries in y-axis
            if x1[i, 1] < min_y:
                x1[i, 1] = x1[i, 1] + L_y
            elif x1[i, 1] > max_y:
                x1[i, 1] = x1[i, 1] - L_y
            
            if x2[i, 1] < min_y:
                x2[i, 1] = x2[i, 1] + L_y
            elif x2[i, 1] > max_y:
                x2[i, 1] = x2[i, 1] - L_y
            
            # Saving positions/timestep
            #particles[2 * i, j, :] = x1[i, :]
            #particles[2 * i + 1, j, :] = x2[i, :]        
        

        
        # Saving positions/timestep, no boundary conditions
        particles[2 * i, j, :] = x1[i, :]
        particles[2 * i + 1, j, :] = x2[i, :]

    # Update part
    x1_timestep = x1
    x2_timestep = x2
    
    

# Visualizing trajectory
write_xyz(particles, str(M) + "xHO_2D_LJ.xyz")

E = U_lj + U_h + K

import matplotlib.pyplot as plt
f = plt.figure(figsize=(20, 20))

plt.plot(U_lj)
plt.plot(U_h)
plt.plot(K)
plt.plot(E)

':)'
