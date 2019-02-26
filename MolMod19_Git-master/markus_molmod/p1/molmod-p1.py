import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
from numpy import linalg as LA
#%matplotlib inline

class Ndim_harmonic_oscillator:
    """
    Molecular Modeling, 2/2017
    Markus Rauhalahti, Juho Karhu, Alex Odiyo

    Simulates the dynamics of a m-dimensional harmonic oscillator with n-particles
    using the leapfrog algorithm.

    dim = number of dimensions
    n   = number of particles

    Input parameters:

    r0   : 2*dim array of particles
    v0   : 2*dim array of initial velocities

    r_eq : equilibrium distance, scalar
    k    : spring constant, scalar

    dt   : timestep
    max_t: simulation time

    """

    def __init__(self, r0, r_eq, v0, k, dt, max_t):
        self.r0 = r0
        self.r_eq = r_eq
        self.v0 = v0
        self.k = k
        self.dt = dt
        self.max_t = max_t


        self.n_steps = int(self.max_t/self.dt)
        self.n   = r0.shape[0]
        self.dim = r0.shape[1]

        # This is used for outputting information
        self.output_filename = "{}p_{}dim_{}dt_{}steps".format(self.n, self.dim, self.dt,self.n_steps)
        # plotting_interval: frequency of writing coordinates
        self.plot_interval = 10


    def variables_ok(self):
        # Check that input arrays have a proper form
        try:
            self.n == 2

        except:
            print("Simulation of only two particles allowed")
            return False
        try:
            self.r0.shape == self.v0.shape
        except:
            print("Initial coordinate array and velocity have different dimensions!")
            return False
        for var in [self.k, self.r_eq]:
            try:
                k > 0.
            except:
                print("Spring constant and eq. distance must be > 0")
                return False
        return True

    def write_coordinates(self,r,step):
        outfile = self.output_filename + ".xyz"
        with open(outfile, "a") as out:
            out.write("{}\n".format(self.n))
            out.write("Step {}/{}\n".format(step,self.n_steps))
            out.write("H {}\n".format(np.array2string(r[0,:],precision=5)[1:-1]))
            out.write("H {}\n".format(np.array2string(r[1,:],precision=5)[1:-1]))

    def force(self, r):
        """
        Calculates the force
        F_{ij} = -\nabla.U
               = k(r_[ij}-r_{ij,eq})
        Input:
        - r : n*dim coordinate array
        Output:
        - F : n*dim array with forces
        """
        force = np.zeros((self.n, self.dim))

        eq_vec = (r[0,:] - r[1,:])/LA.norm((r[0,:] - r[1,:])) * self.r_eq
        force[0,:] = -k*(r[0,:] - r[1,:] -  eq_vec)
        force[1,:] = -force[0,:]
        # Matches analytical with force/2
        return force

    def potential(self, r):
        """
        Calculates the potential energy
        U =  0.5*k(r_{ij}-r_{ij,eq})^2
        Input:
        - r : n*dim coordinate array
        Output:
        - potential energy, float
        """
        potential = 0.5*k*(LA.norm(r[0,:] - r[1,:]) - self.r_eq)**2
        return potential

    def kinetic(self, v):
        """
        Calculates the kinetic energy
        T_{ij} = 0.5*m*v**2, m=1
        Input:
        - v : n*dim speed array
        Output:
        - kinetic energy, float
        """
        velocity_vector = LA.norm(v, axis=1)
        kinetic = 0.5*np.dot(velocity_vector.T, velocity_vector)
        return kinetic
    """
    Runge-Kutta, not used
    def dv(self, r, dt):
        k1=self.force(r)
        k2=self.force(r+k1*dt/2)
        k3=self.force(r+k2*dt/2)
        k4=self.force(r+k3*dt)
        return dt/6*(k1+2*k2+2*k3+k4)
    """
    def leapfrog(self):

        # Initialize n*dim*t_step arrays for position, velocity, and force
        r = np.zeros((self.n, self.dim, self.n_steps))

        v = np.zeros((self.n, self.dim, self.n_steps))
        F = np.zeros((self.n, self.dim, self.n_steps))
        # Initialize n*n_steps arrays for kinetic potential energy
        T = np.zeros(self.n_steps)
        U = np.zeros(self.n_steps)

        # Set the initial conditions
        r[:,:,0] = self.r0
        v[:,:,0] = self.v0
        F[:,:,0] = self.force(r0)
        T[0] = self.kinetic(v0)
        U[0] = self.potential(r0)

        # Run leapfrog for wanted number of steps
        for i in range(self.n_steps-1):
            r[:,:,i+1] = r[:,:,i] + v[:,:,i]*dt
            v[:,:,i+1] = v[:,:,i] + self.force(r[:,:,i+1])*dt
            #v[:,:,i+1] = v[:,:,i] + self.dv(r[:,:,i+1],dt) Runge-Kutta call, not used
            T[i+1] = self.kinetic(v[:,:,i+1])
            U[i+1] = self.potential(r[:,:,i+1])
            if i%self.plot_interval: # plot every n'th
                self.write_coordinates(r[:,:,i],i)
        return r,T,U


    def plot_distance(self, r,show_analytical=True):
        t = np.arange(self.n_steps)*self.dt
        dr=(r[1,:,:]-r[0,:,:])
        dr_sc=LA.norm(dr, axis=0) - self.r_eq
        plt.figure(1)
        plt.plot(t,dr_sc,label="Numerical")
        if show_analytical:
            plt.plot(t,self.analytical(t),label="Analytical")
            std_numeric_analytical = np.std(dr_sc-self.analytical(t))
            print("Standard deviation between numerical and analytical solution: {0:.2e}".format(std_numeric_analytical))

        plt.xlabel("Time ")
        plt.ylabel("Distance of particles")
        plt.legend(loc="upper left", bbox_to_anchor=(1,1))
        #plt.savefig(self.output_filename + "-distance.png")

    def analytical(self, t):
        omega = np.sqrt(self.k)*np.sqrt(2)
        r0_ = LA.norm(self.r0[0,:]-self.r0[1,:]) - self.r_eq
        v0_ = LA.norm(self.v0[0,:]-self.v0[1,:])
        return r0_*np.cos(omega*t) + v0_/omega * np.sin(omega*t)

    def plot_energy(self,T,U):
        plt.figure(2)
        t = np.arange(self.n_steps)*dt
        plt.plot(t,T-T[0],label="Kinetic energy")
        plt.plot(t,U-U[0], label = "Potential energy")
        plt.plot(t,T+U - (T+U)[0], label="Total energy")
        plt.xlabel("Time")
        plt.ylabel("Energy")
        plt.legend(loc="upper left", bbox_to_anchor=(1,1))
        #plt.savefig(self.output_filename + "-energies.png")
    def print_parameters(self):
        print("Parameters:")
        print("- Equilibrium distance: {}".format(self.r_eq))
        print("- Spring constant: {}".format(self.k))
        print("- Timestep: {}".format(self.dt))
        print("- Simulation time: {}".format(self.max_t))

    def metrics(self, r,T,U):
        std_energy = np.std(T+U)
        print("Standard deviation for total energy: {0:.2e}".format(std_energy))

    def run(self):
        print("Harmonic oscillator with the leapfrog algorithm.")
        self.print_parameters()
        if self.variables_ok():
            r, t, u = self.leapfrog()
            self.plot_distance(r)
            self.plot_energy(t,u)
            self.metrics(r,t,u)
        else:
            print("There is something wrong with the parameters")

r0 = np.array([[0., 0., -.3],
               [0., 0., 1.,]])
v0 = np.array([[0., 0.1,0.],
               [0., 0., 0.]])
k = 1
r_eq = 2
dt = 0.005
max_t = 100

calc1 = Ndim_harmonic_oscillator(r0=r0, r_eq=r_eq,v0=v0, k=k,dt=dt,max_t = max_t)
calc1.run()
