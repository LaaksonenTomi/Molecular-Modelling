# READ ME - INSTRUCTIONS

# 1. To start the simulation you only need to run this script.
# 2. To configure the simulations ONLY EDIT section labelled as PARAMETERS.
# 3. To change the output file format change the Truth Values in the OUTPUT OPTIONS section. Only one file can be written at a time.
# 4. At the very bottom in the DEFINE SCRIPT section you may remove a hashtag to perform a script

import numpy as np
import random
import time

#_________MODULES_________

starttime = time.time()

# PARAMETERS
# Give initial values here for forces and particles below

L = 25 # Box side length
rcl = 1 # Latice Cut Off
mass = 1  # masses of particles
r_0 = 1  # equilibrium bond length for harmonic force
U = 100  # some parameter for harmonic force
sigma = 1  # Lennard jones parameter
epsilon = 1  # Lennard jones parameter
LJ_cutoff =  2.5*sigma
periodic_boundaries = True  # using periodic boundaries or not?
xwidth = L  # boundary sizes
ywidth = L
zwidth = L
limits = np.array([xwidth, ywidth, 2 ** 63])  # 2**63 or any large number representing infinite dimension!

# Give simulation parameters below
dt = 0.0001  # timestep
steps = 1000  # number of steps
n = 60      # number of particles should be even
saveInterval = 10 # the interval of saved simulation steps
ke = 32     # this is the 'temperature' that is tied to the thermostat

# END_PARAMETERS

# OUTPUT_OPTIONS
vmd = True
csv = False
if vmd:
    file = open("vmd_list.xyz", "w")
elif csv:
    file = open("csv.csv", "w")

#END_OUTPUT_OPTIONS

# DO NOT CHANGE
# initialise arrays for values with n times 3dimensional coordinates, type is a 64bit float
pos = np.zeros((n, 3), dtype=np.float64)
vel = np.zeros((n, 3), dtype=np.float64)
forces = np.zeros((n, 3), dtype=np.float64)


def generate_random_with_enough_spacing(n_generate, positions, velocities, mindist=1.1 * sigma, maxdist=2.5 * sigma):
    made = 0
    while made < n_generate:
        random_position = np.random.rand(3)
        random_position[0] *= xwidth
        random_position[1] *= ywidth
        random_position[2] *= zwidth  # place particles inside boundaries
        # test if particle too far or close to other particles
        ok = False
        for i9 in range(n_generate):
            distance = np.linalg.norm(random_position - positions[i9])
            if distance < mindist:
                # print("Particle {0} too close to some particle @ {1}".format(made, random_position))
                break
            elif maxdist > distance > mindist:
                ok = True
            if i9 == n_generate - 1 and ok:  # no collision after checks and withint good range to atleast one particle
                positions[made] = random_position  # place particle to the position
                velocities[made] = np.array([random.uniform(-1,1), random.uniform(-1,1), random.uniform(-1,1)])
                made += 1


def harmonic_force(pos1, pos2):
    vector_distance = distancevector(pos1, pos2)
    scalar_distance = np.linalg.norm(vector_distance)
    scalar_force = 2 * U * (scalar_distance - r_0)
    vector_force = scalar_force * vector_distance / scalar_distance
    return vector_force


def lennard_jones_force(pos1, pos2):

    vector_distance = distancevector(pos2, pos1)
    scalar_distance = np.linalg.norm(vector_distance)
    if scalar_distance > LJ_cutoff:
        return np.zeros(3)  # return a zerovector
    r = scalar_distance
    scalar_force = 4 * epsilon * (((sigma / r) ** 12) - ((sigma / r) ** 6))
    vector_force = vector_distance / scalar_distance * scalar_force
    return vector_force


def output_vmd():
    file.write("{0} \n\n".format(n))
    # print(pos)
    for particle in pos:
        # print(particle)
        file.write("N {0} {1} {2}\n".format(particle[0], particle[1], particle[2]))


def output_csv():
    global forces
    # implement output for kinetic and potential energy
    # calculate kinetic energy:
    E_kinetic = 0.5 * mass * vel**2
    kinetic_particle = np.sum(E_kinetic, axis=1)
    kinetic_sum = np.sum(kinetic_particle)
    E_potential = forces * vel
    potential_particle = np.sum(E_potential, axis=1)
    potential_sum = np.sum(potential_particle)

    temp = (2.0/3.0)*(kinetic_sum/n)
    file.write("{0};{1};{2}\n".format(kinetic_sum,potential_sum,temp))  # separated by ; (can be read by excel)
    #print E_potential, potential_particle, potential_sum


def distancevector(loc1, loc2, boundaries=True):
    difference = loc2 - loc1
    if not boundaries:  # when no PBC needed
        return difference
    else:  # when periodic boundary conditions are in effect:
        if abs(difference[0]) > xwidth / 2:
            difference[0] = - np.sign(difference[0]) * (xwidth - difference[0])

        if abs(difference[1]) > ywidth / 2:
            difference[1] = - np.sign(difference[1]) * (ywidth - difference[1])

        if abs(difference[2]) > zwidth / 2:
            difference[2] = - np.sign(difference[2]) * (zwidth - difference[1])
        return difference


def apply_periodic_boundary_conditions(position):
    if position[0] < 0:
        position[0] = xwidth - position[0]
    elif position[0] > xwidth:
        position[0] = position[0] - xwidth
    if position[1] < 0:
        position[1] = ywidth - position[1]
    elif position[1] > ywidth:
        position[1] = position[1] - ywidth
    if position[2] < 0:
        position[2] = zwidth - position[2]
    elif position[2] > zwidth:
        position[2] = position[2] - zwidth


def calculate_LJ_between_all():
    # add LJ forces for every particle
    for i in range(n - 1):  # from first to second last
        for i2 in range((i + 1), n):  # from current particle to last
            LJ_between = lennard_jones_force(pos[i], pos[i2])
            forces[i] += LJ_between
            forces[i2] -= LJ_between  # "add" the opposite vector


def mimimise():
    global pos
    global forces
    # start with zero forces:
    forces = np.zeros((n, 3), dtype=np.float64)
    # calculate forces
    calculate_LJ_between_all()
    # steepest
    maxforce_scalar=np.max(np.linalg.norm(forces, axis=1))
    pos += forces / (mass * maxforce_scalar) * dt


def progress_list():
    global pos
    global forces
    global vel
    global alpha

    calculate_LJ_between_all()
    vel = (vel + forces / mass * dt) * alpha
    pos += vel * dt

    for particle in pos:
        i = 0
        for i in range (len(pos)) :
            apply_periodic_boundary_conditions(pos[i])
            i += 1

def resamplekin_sumnoises(nn):
    if (nn==0.0):
        return 0.0
    elif (nn==1.0):
        rr=gasdev()
        return rr*rr
    elif (nn%2==0):
        return 2.0*gamdev(nn/2)
    else:
        rr=gasdev()
        return 2.0*gamdev((nn-1)/2) + rr*rr


def gasdev():
    iset=0
    rsq=0.0
    if (iset == 0):
        while (rsq >= 1.0 or rsq == 0.0):
            v1=2.0*random.random()-1.0
            v2=2.0*random.random()-1.0
            rsq=v1*v1+v2*v2
        fac=np.sqrt(-2.0*np.log(rsq)/rsq)
    	gset=v1*fac
    	iset=1
    	return v2*fac
    else:
        iset=0
        return gset

def gamdev(ia):
    if (ia < 6):
        x=1.0
        for j in range(0,int(ia)):
            x *= random.random()
        x = -np.log(x)
    else:
        while (random.random() > e):
            e=(1.0+y*y)*np.exp(am*np.log(x/am)-s*y)
            while (x<=0.0):
                y=v2/v1
                am=ia-1
                s=np.sqrt(2.0*am+1.0)
                x=s*y+am
                while (v1*v1+v2*v2 > 1.0):
                    v1=random.random()
                    v2=2.0*random.random()-1.0
    return x

def resamplekin(kk,sigma,ndeg,taut):
    # kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
    # sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
    # ndeg:  number of degrees of freedom of the atoms to be thermalized
    # taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
    if (taut>0.1):
        factor=np.exp(-1.0/taut)
    else:
        factor = 0.0
    rr = gasdev()
    ke_adj = kk + (1.0-factor)*(sigma*(resamplekin_sumnoises(ndeg-1)+rr**2)/ndeg-kk) + 2.0*rr*np.sqrt(kk*sigma/ndeg*(1.0-factor)*factor)
    #print ke_adj
    return ke_adj

def progress_linked_cell():
    global pos
    global forces
    global vel
    global alpha

    # Create a linked array for particles in boxes
    linkedArray = []
    linkedArrayLength = []
    indexBox = 0

    #Put particles into boxes
    for indexBox in range (len(boxPairs)) :
        boxArray = []

        #Go through particles
        indexParticle = 0
        for indexParticle in range (len(pos)) :
            if ((boxPosition[indexBox][0] <= pos[indexParticle][0] < (boxPosition[indexBox][0] + rcl)) and (boxPosition[indexBox][1] <= pos[indexParticle][1] < (boxPosition[indexBox][1] + rcl)) and (boxPosition[indexBox][2] <= pos[indexParticle][2] < (boxPosition[indexBox][2] + rcl))):
                boxArray.append(indexParticle)
                indexParticle += 1

        linkedArray.append(boxArray)
        length = len(boxArray)
        linkedArrayLength.append(length)
        indexBox += 1

    #Go through particles and calculate forces
    indexBox = 0
    while indexBox < numberBoxes:

        #calculate forces for particles
        indexParticle = 0
        for indexParticle in range(linkedArrayLength[indexBox] - 1):

            #particles in the same box
            indexParticle2 = 1
            for indexParticle2 in range(linkedArrayLength[indexBox]):
                LJ_between = lennard_jones_force(pos[linkedArray[indexBox][indexParticle]], pos[linkedArray[indexBox][indexParticle2]])
                indexParticle2 += 1
                forces[indexParticle] += LJ_between
                forces[indexParticle2] -= LJ_between


            #particles in the neighbouring boxes for the same particle
            indexNeightbour = 0
            neightboursAfterFirst = 12
            for indexNeightbour in range(neightboursAfterFirst):

                indexParticle2 = 0
                while indexParticle2 < linkedArrayLength[indexNeightbour]:

                    LJ_between = lennard_jones_force(pos[linkedArray[indexBox][indexParticle]], pos[linkedArray[indexNeightbour][indexParticle2]])
                    forces[indexParticle] += LJ_between
                    forces[linkedArray[indexNeightbour][indexParticle2]] -= LJ_between

                    indexParticle2 += 1

                indexNeightbour += 1

            indexParticle += 1

        indexBox += 1
    #Replace old position
    vel = (vel + forces / mass * dt)*alpha
    pos += vel * dt
    for particle in pos:
        i = 0
        for i in range (len(pos)) :
            apply_periodic_boundary_conditions(pos[i])
            i += 1
    #print(linkedArray)

#_________INI_________

#divide space i < (L/rcl)**3

#Splitting the cube into boxes
numberBoxes = int((L/rcl)**3)

boxesSplit = []
zBoxID = 0
while (zBoxID < L/rcl):
    yBoxID = 0
    while (yBoxID < L/rcl):
        xBoxID = 0
        while (xBoxID < L/rcl):
            boxesSplit.append((xBoxID,yBoxID,zBoxID))
            xBoxID += 1
        yBoxID += 1
    zBoxID += 1

#Create a list of BoxPairs
# [(L/rcl)**3, 13]
i = 0

boxPairs = []
while i < numberBoxes:
    xBoxID = boxesSplit[i][0]
    yBoxID = boxesSplit[i][1]
    zBoxID = boxesSplit[i][2]
    boxID1 =  (xBoxID    ,yBoxID + 1,zBoxID    )
    boxID2 =  (xBoxID + 1,yBoxID + 1,zBoxID    )
    boxID3 =  (xBoxID + 1,yBoxID    ,zBoxID    )
    boxID4 =  (xBoxID + 1,yBoxID - 1,zBoxID    )
    boxID5 =  (xBoxID - 1,yBoxID + 1,zBoxID + 1)
    boxID6 =  (xBoxID - 1,yBoxID    ,zBoxID + 1)
    boxID7 =  (xBoxID - 1,yBoxID - 1,zBoxID + 1)
    boxID8 =  (xBoxID    ,yBoxID + 1,zBoxID + 1)
    boxID9 =  (xBoxID    ,yBoxID    ,zBoxID + 1)
    boxID10 = (xBoxID    ,yBoxID - 1,zBoxID + 1)
    boxID11 = (xBoxID + 1,yBoxID + 1,zBoxID + 1)
    boxID12 = (xBoxID + 1,yBoxID    ,zBoxID + 1)
    boxID13 = (xBoxID + 1,yBoxID - 1,zBoxID + 1)
    boxIDs = ([boxID1,boxID2,boxID3,boxID4,boxID5,boxID6,boxID7,boxID8,boxID9,boxID10,boxID11,boxID12,boxID13])
    #print(boxIDs)
    boxPairs.append(boxIDs)
    i += 1

#Store displament vectors for each of the 13 boxes in the cell; vector: Box to Box
displacementVectors = []
displacements1 = [0,rcl,0]
displacements2 = [rcl,rcl,0]
displacements3 = [rcl,0,0]
displacements4 = [rcl,-rcl,0]
displacements5 = [-rcl,rcl,rcl]
displacements6 = [-rcl,0,rcl]
displacements7 = [-rcl,-rcl,rcl]
displacements8 = [0,rcl,rcl]
displacements9 = [0,0,rcl]
displacements10 = [0,-rcl,rcl]
displacements11 = [rcl,rcl,rcl]
displacements12 = [rcl,0,rcl]
displacements13 = [rcl,-rcl,rcl]
displacements = (displacements1, displacements2, displacements3, displacements4, displacements5, displacements6, displacements7, displacements8, displacements9, displacements10, displacements11,displacements12,displacements13)
displacementVectors.append(displacements)


#create a list of corner vectors [(L/rcl)**3]
# [(L/rcl)**3] vectors
boxCornerVectors = []
z = 0
while (z < L):
    y = 0
    while (y < L):
        x = 0
        while (x < L):
            boxCornerVectors.append([x,y,z])
            x += rcl
        y += rcl
    z += rcl

# cornerVectors.append(x,y,z)

#create an Array for all (L/rcl)**3
#[(x,y,z)] cartesian
boxPosition = []
z = 0
while (z < L):
    y = 0
    while (y < L):
        x = 0
        while (x < L):
            boxPosition.append((x,y,z))
            x += rcl
        y += rcl
    z += rcl

#_________MAIN_________

print("Generating {0} initial particle positions".format(n))
# generate_particles(n, pos, vel) # generate 10 random particles
generate_random_with_enough_spacing(n, pos, vel)

step = 0
# simulation loop
print("Going to simulate {0} steps!".format(steps))
while step < steps:

    # output state
    if ((step % saveInterval) == 0):
        print("Step {0}".format(step))
        if vmd:
            output_vmd()
        else:
            output_csv()

    # zero the forces before calculating new ones!
    forces = np.zeros((n, 3), dtype=np.float64)
    # calculate harmonic forces for every "molecule"
    #for i in range(int(n/2)):
    #    force_between_pair = harmonic_force(pos[i], pos[i+1])
    #    forces[i*2] = force_between_pair
    #    forces[i*2+1] = -force_between_pair
    # add LJ forces for every particle
    #for i in range(n - 1):  # from first to second last
    #    for i2 in range((i + 1), n):  # from current particle to last
    #        LJ_between = lennard_jones_force(pos[i], pos[i2])
    #        forces[i] += LJ_between
    #        forces[i2] -= LJ_between  # "add" the opposite vector
    # update velocity
    #vel = vel + forces / mass * dt
    # update positions
    #pos += vel * dt
    # apply periodic boundary conditions
    #for i in range(n):
    #    if periodic_boundaries:
    #        apply_periodic_boundary_conditions(pos[i])
    #print("Step {0}".format(step))

    # DEFINE SCRIPT

    #mimimise()

    tau = 1
    ndeg = 6
    E_kinetic = 0.5 * mass * vel**2
    kinetic_particle = np.sum(E_kinetic, axis=1)
    kinetic_sum = np.sum(kinetic_particle)
    ke_adj = resamplekin(kinetic_sum,ke,3,tau)
    alpha = ke_adj/kinetic_sum
    #alpha = 1
    #progress_list()
    progress_linked_cell()
    step += 1
print("\nFinished! Time: {0} s".format(time.time() - starttime))
