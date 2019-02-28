import numpy as np
import matplotlib.pyplot as plt
#from numpy import linalg as linis

import random
import math
import sys

#from RandomClass import Random
#from collections import Counter


#real parameters
k = 1.0
r_0 = 1.0
mass_p1 = 1.0
mass_p2 = 1.0
t = 0.0
dt = 0.001
r = 1.0
E_kin_p1 = 0.0
E_kin_p2 = 0.0
E_kin_tot = 0.0
E_pot = 0.0
E_tot = 0.0

#real dimension
particlesNumber = 2
position_p1 = np.array((0,0,0))
position_p2 = np.array((0,0,0))
velocity_p1 = np.array((0,0,0))
velocity_p2 = np.array((0,0,0))
accelleration_p1 = np.array((0,0,0))
accelleration_p2 = np.array((0,0,0))
lin_mom_p1 = np.array((0,0,0))
lin_mom_p2 = np.array((0,0,0))
ang_mom_p1 = np.array((0,0,0))
ang_mom_p2 = np.array((0,0,0))
center_mass = np.array((0,0,0))
total_lin_mom = np.array((0,0,0))
total_ang_mom = np.array((0,0,0))

#input
print("Please give positions for particles")
position_p1 = np.array((1,1,1))
position_p2 = np.array((-1,-1,-1))

velocity_p1 = np.array((1,1,1)) #random
velocity_p2 = np.array((-1,-1,-1)) #random

#Initializa
mass_p1 = 1.0
mass_p2 = 1.0
r_0 = 1.0
eq_r = 1.0
k = 1.0
t = 0.0
dt = 0.0001

r = np.dot(position_p1, position_p2)
center = np.add(position_p1,0.5*np.subtract(position_p2,position_p1))

#file preperation

#Calculate momenta
lin_mom_p1 = velocity_p1*mass_p1
lin_mom_p2 = velocity_p2*mass_p2
total_lin_mom = lin_mom_p1 + lin_mom_p2

ang_mom_p1 = (position_p1 - center)*velocity_p1
ang_mom_p2 =(position_p2 - center)*velocity_p2
total_ang_mom = ang_mom_p1 + ang_mom_p2
#print momenta

#Calculate energies
E_kin_p1 = 0.5*mass_p1*velocity_p1**2
E_kin_p2 = 0.5*mass_p2*velocity_p2**2
E_kin_tot = E_kin_p1 + E_kin_p2
E_pot = (r-eq_r)**2
E_tot = E_kin_tot + E_pot
#print energies

#Time evolution
i = 0
n = 11000
t = t + dt
while i < n:
    
    #file.write("%s \n\n" % particlesNumber)  # First step to xyz file
    #Calculate acceleration
    acceleration_p1 = -k*2*(r - eq_r)*((position_p1 - position_p2) / r) / mass_p1
    acceleration_p2 = -k*2*(r - eq_r)*((position_p2 - position_p1) / r) / mass_p2

    #Calculate new velocity
    velocity_p1 = velocity_p1 + acceleration_p1*dt
    velocity_p2 = velocity_p2 + acceleration_p2*dt

    #Calculate new position
    position_p1 = position_p1 + velocity_p1*dt
    position_p2 = position_p2 + velocity_p2*dt

    #Calculate distance
    r = np.dot(position_p1, position_p2)
    center = np.add(position_p1,0.5*np.subtract(position_p2,position_p1))

    #Calculate momenta
    lin_mom_p1 = velocity_p1*mass_p1
    lin_mom_p2 = velocity_p2*mass_p2
    total_lin_mom = lin_mom_p1 + lin_mom_p2

    ang_mom_p1 = (position_p1 - center)*velocity_p1
    ang_mom_p2 =(position_p2 - center)*velocity_p2
    total_ang_mom = ang_mom_p1 + ang_mom_p2
    #print momenta

    #Calculate energies
    E_kin_p1 = 0.5*mass_p1*np.linalg.norm(velocity_p1)**2
    E_kin_p2 = 0.5*mass_p2*np.linalg.norm(velocity_p2)**2
    E_kin_tot = E_kin_p1 + E_kin_p2
    E_pot = (r-eq_r)**2 #*U0
    E_tot = E_kin_tot + E_pot
    #print (E_kin_p1, E_kin_p2, E_kin_tot, E_pot, E_kin_tot)
    print(position_p1, position_p2)
    i += 1

velocity_p1 = numpy.zeros((0,0,0)) #random
velocity_p2 = numpy.zeros((0,0,0)) #random

#Initializa
mass_p1 = 1.0
mass_p2 = 1.0
r_0 = 1.0
eq_r = 1.0
k = 1.0
t = 0.0
dt = 0.01

r = numpy.dot(position_p1, position_p2)
center = numpy.add(position_p1,0.5*numpy.subtract(position_p2,position_p1))

#file preperation

#Calculate momenta
lin_mom_p1 = velocity_p1*mass_p1
lin_mom_p2 = velocity_p2*mass_p2
total_lin_mom = lin_mom_p1 + lin_mom_p2

ang_mom_p1 = (position_p1 - center)*velocity_p1
ang_mom_p2 =(position_p2 - center)*velocity_p2
total_ang_mom = ang_mom_p1 + ang_mom_p2
#print momenta

#Calculate energies
E_kin_p1 = 0.5*mass_p1*velocity_p1**2
E_kin_p2 = 0.5*mass_p2*velocity_p2**2
E_kin_tot = E_kin_p1 + E_kin_p2
E_pot = (r-eq_r)**2
E_tot = E_kin_tot + E_pot
#print energies


#Time evolution
i = 0
n = 11000
t = t + dt
while i < n:
    
    #Calculate acceleration
    acceleration_p1 = -k*2*(r - eq_r)*((position_p1 - position_p2) / r) / mass_p1
    acceleration_p2 = -k*2*(r - eq_r)*((position_p2 - position_p1) / r) / mass_p2

    #Calculate new velocity
    velocity_p1 = velocity_p1 + acceleration_p1*dt
    velocity_p2 = velocity_p2 + acceleration_p2*dt

    #Calculate new position
    position_p1 = position_p1 + velocity_p1*dt
    position_p2 = position_p2 + velocity_p2*dt

    #Calculate distance
    r = (position_p1[:,:,None,:]*position_p2[...,None]).sum(1)
    center = numpy.add(position_p1,0.5*numpy.subtract(position_p2,position_p1))

    #Calculate momenta
    lin_mom_p1 = velocity_p1*mass_p1
    lin_mom_p2 = velocity_p2*mass_p2
    total_lin_mom = lin_mom_p1 + lin_mom_p2

    ang_mom_p1 = (position_p1 - center)*velocity_p1
    ang_mom_p2 =(position_p2 - center)*velocity_p2
    total_ang_mom = ang_mom_p1 + ang_mom_p2
    #print momenta

    #Calculate energies
    E_kin_p1 = 0.5*mass_p1*velocity_p1**2
    E_kin_p2 = 0.5*mass_p2*velocity_p2**2
    E_kin_tot = E_kin_p1 + E_kin_p2
    E_pot = (r-eq_r)**2
    E_tot = E_kin_tot + E_pot
    #print energies


# In[ ]:




