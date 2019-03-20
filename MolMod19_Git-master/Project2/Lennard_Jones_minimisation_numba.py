import numpy as np
import time
import numba
from numba import njit, jit, prange

# if numba causes problems, remove line above and every @njit line
# possibly also change prange to range!

starttime = time.time()

# some testings in this file
# give initial values here for forces and particles

mass = 1  # masses of particles
r_0 = 1  # equilibrium bond length for harmonic force
U = 100  # some parameter for harmonic force
sigma = 1  # Lennard jones parameter / "bond" length
epsilon = 1  # Lennard jones parameter / well depth
LJ_cutoff = 2.5 * sigma
periodic_boundaries = True  # using periodic boundaries or not?
xwidth = 25  # boundary sizes
ywidth = 25
limits = np.array([xwidth, ywidth, 2 ** 63])  # 2**63 or any large number representing infinite dimension!
# simulation parameters below

dt = 0.1 * r_0  # timestep / minimisation movement step
steps = 1000  # number of steps, 1000 and 400 particles takes 70s
n = 400  # number of particles should be even
# even number for building "diatomic molecules"

# used files
vmd = False
csv = True
if vmd:
    file = open("vmd.xyz", "w")
elif csv:
    file = open("csv.csv", "w")

# initialise arrays for values with n times 3dimensional coordinates, type is a 64bit float
pos = np.zeros((n, 3), dtype=np.float64)
vel = np.zeros((n, 3), dtype=np.float64)
forces = np.zeros((n, 3), dtype=np.float64)


def generate_random_with_enough_spacing(n_generate, positions, velocities, mindist=0.99 * sigma, maxdist=2.3 * sigma):
    made = 0
    while made < n_generate:
        random_position = np.random.rand(3)
        random_position[0] *= xwidth
        random_position[1] *= ywidth
        random_position[2] = 0  # place particles inside boundaries
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
                velocities[made] = np.zeros(3)
                made += 1


# function for harmonic force
def harmonic_force(pos1, pos2):
    vector_distance = distancevector(pos1, pos2)
    scalar_distance = np.linalg.norm(vector_distance)
    scalar_force = 2 * U * (scalar_distance - r_0)
    vector_force = scalar_force * vector_distance / scalar_distance
    return vector_force


@njit
def lennard_jones_force(pos1, pos2):
    vector_distance = distancevector(pos1, pos2)
    scalar_distance = np.linalg.norm(vector_distance)
    if scalar_distance > LJ_cutoff:
        return np.zeros(3)  # return a zerovector
    r = scalar_distance
    scalar_force = 4 * epsilon * (7 * (sigma ** 6) / (r ** 7) - (12 * (sigma ** 12) / (r ** 13)))
    vector_force = vector_distance / scalar_distance * scalar_force
    return vector_force


@njit
def lennard_jones_potential(pos1, pos2):
    vector_distance = distancevector(pos1, pos2)
    scalar_distance = np.linalg.norm(vector_distance)
    if scalar_distance > LJ_cutoff:
        return 0
    r = scalar_distance
    scalar_potential = 4 * epsilon * ((sigma / r) ** 12 - ((sigma / r) ** 6))
    return scalar_potential


def output_vmd():
    file.write("{0} \n\n".format(n))
    # print(pos)
    for particle in pos:
        # print(particle)
        file.write("N {0} {1} {2}\n".format(particle[0], particle[1], particle[2]))


def output_csv():
    # implement output for kinetic and potential energy
    # calculate kinetic energy:
    # E_kinetic = np.cumsum(0.5 * mass * np.power(vel, 2))
    # E_potential = 0
    energy = calculate_LJ_potential_all(pos)
    file.write("{0};\n".format(energy))  # separated by ; (can be read by excel)


@njit
def distancevector(loc1, loc2, boundaries=False):
    difference = loc2 - loc1
    if not boundaries:  # when no PBC needed
        return difference
    else:  # when periodic boundary conditions are in effect:
        if abs(difference[0]) > xwidth / 2:
            difference[0] = - np.sign(difference[0]) * (xwidth - difference[0])

        if abs(difference[1]) > ywidth / 2:
            difference[1] = - np.sign(difference[1]) * (xwidth - difference[1])
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


@njit(parallel=True)
def calculate_LJ_between_all(positions):
    # add LJ forces for every particle
    forces_LJ = np.zeros((n, 3), dtype=np.float64)
    for i in prange(n - 1):  # from first to second last
        for i2 in prange((i + 1), n):  # from current particle to last
            LJ_between = lennard_jones_force(positions[i], positions[i2])
            forces_LJ[i] += LJ_between
            forces_LJ[i2] -= LJ_between  # "add" the opposite vector
    return forces_LJ


@njit(parallel=True, fastmath=True)
def calculate_LJ_potential_all(positions):
    potential_sum = 0
    for i in prange(n - 1):
        for i2 in prange((i + 1), n):
            potential_sum += lennard_jones_potential(positions[i], positions[i2])
    return potential_sum


def mimimise():
    global pos
    global forces
    # start with zero forces:
    # calculate forces
    forces = calculate_LJ_between_all(pos)
    # steepest
    maxforce_scalar = np.max(np.linalg.norm(forces, axis=1))
    pos += forces / (mass * maxforce_scalar) * dt


print("Generating {0} initial particle positions".format(n))
# generate_particles(n, pos, vel) # generate 10 random particles
generate_random_with_enough_spacing(n, pos, vel)

step = 0
# simulation loop
print("Going to simulate {0} steps!".format(steps))
while step < steps:

    # output state
    if vmd:
        output_vmd()
        print(step / steps)
    else:
        output_csv()
        print(step / steps)
    mimimise()
    step += 1
print("\nFinished! Time: {0} s".format(time.time() - starttime))

# zero the forces before calculating new ones!
# forces = np.zeros((n, 3), dtype=np.float64)
# calculate harmonic forces for every "molecule"
# for i in range(int(n/2)):
#     force_between_pair = harmonic_force(pos[i], pos[i+1])
#     forces[i*2] = force_between_pair
#     forces[i*2+1] = -force_between_pair
# add LJ forces for every particle
# for i in range(n - 1):  # from first to second last
#     for i2 in range((i + 1), n):  # from current particle to last
#         LJ_between = lennard_jones_force(pos[i], pos[i2])
#         forces[i] += LJ_between
#         forces[i2] -= LJ_between  # "add" the opposite vector
# update velocity
# vel = vel + forces / mass * dt
# update positions
# pos += vel * dt
# apply periodic boundary conditions
# for i in range(n):
#    if periodic_boundaries:
#        apply_periodic_boundary_conditions(pos[i])
# print("Step {0}".format(step))
