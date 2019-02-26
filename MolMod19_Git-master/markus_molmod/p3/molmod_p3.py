import numpy as np
import matplotlib.pyplot as  plt
#import sys
# For timing
from timeit import default_timer as timer

timedict = {}
timedict["Energies"] = 0
timedict["Forces"] = 0
timedict["List generation"] = 0
timedict["I/O"] = 0

def plot_energies():
    outfile = "N{}_dt{}_gamma{}_T{}".format(N,dt,gamma,T)
    plt.figure(3)
    t = np.arange(S)
    plt.plot(t,kinetic-kinetic[0],label="Kinetic energy")
    plt.plot(t,potential-potential[0], label = "Potential energy")
    plt.plot(t,potential - (potential)[0], label="Total energy")
    plt.xlabel("Step")
    plt.ylabel("Energy")
    plt.title(outfile)
    #plt.legend(loc="upper left", bbox_to_anchor=(1,1))
    plt.legend()
    plt.savefig(outfile + ".pdf")

def print_times():
    time = 0
    for component in timedict:
        time += timedict[component]
    print("Total time spent: ", time)
    print("Setup: ", timedict["Setup"])
    print("I/O: ", timedict["I/O"])
    print("Iterations in total: ", timedict["Iterations"])
    print("Forces: ", timedict["Forces"])
    print("Energies: ", timedict["Energies"])

def closestImage(dq):
    """
    Goes through all vector components of dq, removes/adds L if particles are over the box
    """
    for i in range(3):
        if dq[i] > L/2:
            dq[i] -= L
        elif dq[i] < -L/2:
            dq[i] += L
    return dq

def calcForce(dq, vi, langevin=True):
    """
    Computes forces between two particles interacting via LJ potential
    """
    dq = closestImage(dq)
    r = np.linalg.norm(dq)
    f = -dq* 4*epsilon*(6*sigma**6*r**(-6)- 12*sigma**(12)*r**(-12))
    if langevin:
        f = f - gamma*mass*vi + np.random.normal(0, 2*mass*gamma*kb*T/dt)
    return f


def to_unit_cell(q):
    """
    Takes care of periodic boundary conditions
    """
    for i in range(N):
        for j in range(3):
            if abs(q[i,j] < -L/2):
                print("Ball is out of unit box, break")
            if q[i,j] < - L/2:
                q[i,j] += L
            elif q[i,j] > L/2:
                q[i,j] -= L
    return q

# Functions for the linked cell algorithm

def icell(ix,iy,iz):
    """
    Calculates a serial index for a cell.
    Used to index cells so that one can use a 2d array instead of a 4d array.
    """
    icell = (ix+M)%M + ((iy+M)%M)*M + ((iz+M)%M)*M*M
    return int(icell + M/2 - 2)

def get_lists(q):
    """
    Sorts atoms to cells by constructing two lists:
    1) head: contains one element per cell
    2) lists: list-of-lists, contains a list of atoms for each cell


    As the implementation uses lists as implemented in python, the first list is redundant.
    This is because the list data structure of python is already extendable.
    """
    start = timer()
    #head = [None]*M**3
    lists = [None]*M**3
    for iat, atom_coordinates in enumerate(q):
        # ix,iy,iz for cell
        cellvector = [int(i/cutoff) for i in atom_coordinates]
        cell_index = icell(*cellvector)
        # head list is not used
        #if head[cell_index] == None:
        #    head[cell_index] = iat
        lists[cell_index] = []
        lists[cell_index].append(iat)
    end = timer()
    timedict["List generation"] += end - start
    return lists

def neighbours_list(ix,iy,iz):
    """
    Calculates the neighbours for a cell and appends them to a list
    """
    neigh_cells = []
    neigh_cells.append(icell(ix+1,iy,iz)    )
    neigh_cells.append(icell(ix+1,iy+1,iz)  )
    neigh_cells.append(icell(ix,iy+1,iz)    )
    neigh_cells.append(icell(ix-1,iy+1,iz)  )
    neigh_cells.append(icell(ix+1,iy,iz-1)  )
    neigh_cells.append(icell(ix+1,iy+1,iz-1))
    neigh_cells.append(icell(ix,iy+1,iz-1)  )
    neigh_cells.append(icell(ix-1,iy+1,iz-1))
    neigh_cells.append(icell(ix+1,iy,iz+1)  )
    neigh_cells.append(icell(ix+1,iy+1,iz+1))
    neigh_cells.append(icell(ix,iy+1,iz+1)  )
    neigh_cells.append(icell(ix-1,iy+1,iz+1))
    neigh_cells.append(icell(ix,iy,iz+1)    )
    return neigh_cells

def calc_forces(q,v,lists,F):
    """
    Calculates the forces between atoms within a cell and within the list of neighbouring cells.
    Returns forces
    """
    start = timer()
    for cell_i in lists:
        if cell_i == None:
            continue
        # Sort the list
        cell_i.sort()
        ats = cell_i
        # Calculate forces within the cell
        for ind_i, at_i in enumerate(ats):
            for ind_j in range(ind_i+1, len(ats)): # No self-interaction or double counting
                at_j = ats[ind_j]
                dq = q[at_j,:] - q[at_i,:]
                F[at_i] = calcForce(dq, v[at_i,:])
        # Calculate forces within atoms in a cell and its neighbours

        # Get the neighbours
        cell_vec = [int(i/cutoff) for i in q[ats[0],:]]
        cell_index = icell(*cell_vec)
        neighs = neighbours_list(*cell_vec)

        for j in neighs:
            try:
                cell_j = lists[j]
            except:
                print(j, len(lists))
            if cell_j == None:
                continue
            for ind_i, at_i in enumerate(cell_i):
                for ind_j in range(ind_i+1, len(cell_j)):
                    at_j = cell_j[ind_j]

                    dq = q[at_j,:] - q[at_i,:]
                    F[at_i] = calcForce(dq, v[at_i,:])
    end = timer()
    timedict["Forces"] += end - start
    return F


def compute_energies(q,v):
    """
    Computes the kinetic and potential energies from q, v.

    Input: q, v: N*3 position and velocity numpy arrays
    Output: Tuple of kinetic, potential energies
    """
    start = timer()
    kinetic = 0.5*mass*np.linalg.norm(v)**2
    potential = 0
    for i in range(N):
        for j in range(N):
            if j>i:
                r = np.linalg.norm(closestImage(q[i,:]-q[j,:]))
                potential += 4*epsilon*((sigma/r)**12-(sigma/r)**6)
    end = timer()
    timedict["Energies"] += end - start
    return kinetic, potential

def print_q(k,q):
    """
    Append coordinates to the output file
    """
    start = timer()
    with open(outfile, "a") as out:
        out.write("{} \nStep {}/{} \n".format(N, k, S))
        for at in q:
            out.write("Ar {:.6f} {:.6f} {:.6f}\n".format(at[0],at[1],at[2]))
    end = timer()
    timedict["I/O"] += end - start

def print_parameters():
    print("""
    Molecular Modeling, project 3
    Parameters:
       N = {}
       Steps = {}
       dt = {}
       L = {}
       cutoff = {}
    LJ parameters:
       sigma = {}
       epsilon = {}

    Thermostat parameters:
       T = {}
       Gamma = {}

    Linked cell parameters:
       Number of cells = {}
       Cell length = {}
    """.format(N, S, dt, L, cutoff,sigma, epsilon, T, gamma, M**3, L_cell)
    )
###############################################################################
# Initialization of parameters
###############################################################################


N = 5 # number of particles
S = 100 # steps
dt = 0.001 # timestep
sigma = 0.001 # LJ sigma
epsilon = 0.01 # LJ epsilon
mass = 1
cutoff = 1
L = 5
verbose = True
# Parameters for the Langevin thermostat
T = 300
gamma = 0.1

# Parameters for printing the results
print_every_nth = 1
outfile = "N{}_S{}.xyz".format(N,S)

# Constants
kb = 1.38064852e-23
start = timer()
# Initialization of arrays
q = np.random.uniform(-L/2, L/2, size=(N,3)) # coordinates of particles, [-L/2, L/2]
v = np.random.rand(N,3)
kinetic = np.zeros(S)
potential = np.zeros(S)
dq = np.zeros(3)
F = np.zeros((N,3))

# Linked cell parameters

# Number of cells per edge
M = int(L/cutoff)
# Length of the cell's edge
L_cell = L/M




###############################################################################
# The MD simualtion routine
###############################################################################

if verbose:
    print_parameters()

# Center and remove the total momentum
q = q - np.mean(q)*L
v = v - np.mean(v)
q = to_unit_cell(q)

end = timer()

timedict["Setup"] = end - start

startit = timer()
for k in range(S):
    #print("Round ", k)
    # Results are written every n'th step
    if S%print_every_nth == 0:
        print_q(k,q)
    # Compute energies
    kinetic[k], potential[k] = compute_energies(q,v)
    # Sort atoms into cells, return lists
    lists = get_lists(q)
    # Compute forces of the atoms
    F = calc_forces(q, v, lists, F)
    # Evolve in time using the Leapfrog algorithm
    v += F/mass*dt
    q += v*dt
    # Take care of periodic boundary conditions
    to_unit_cell(q)
endit = timer()
timedict["Iterations"] = endit - startit
print_times()
plot_energies()
