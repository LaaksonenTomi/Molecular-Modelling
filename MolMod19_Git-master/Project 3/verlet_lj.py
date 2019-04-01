import random, numpy as np

outFile = open('lj_3d.xyz', 'w')

Ntotal = 50

eps = 1
sigma = 1

timestep = 0.0001
steps = 1000

systemSize = 5

def pick_coordinates(Ntotal, systemSize):
    global x, y, z
    x = np.zeros(Ntotal)
    y = np.zeros(Ntotal)
    z = np.zeros(Ntotal)

    for i in range(0, Ntotal):
        x[i] = random.uniform((-(systemSize-1)), (systemSize-1))
        y[i] = random.uniform((-(systemSize-1)), (systemSize-1))
        z[i] = random.uniform((-(systemSize-1)), (systemSize-1))

    for i in range(0, Ntotal):
        for j in range(0,Ntotal):
            if (i >= j):        #excludes paris that have already been calculated and the same particles
                continue
            else:
                g = ((x[i] - x[j])**2 + (y[i]-y[j])**2 + (z[i]-z[j])**2)
                r = np.sqrt(g)
		if (r<=1):
		        x[i] = random.uniform((-(systemSize-1)), (systemSize-1))
		        y[i] = random.uniform((-(systemSize-1)), (systemSize-1))
     			z[i] = random.uniform((-(systemSize-1)), (systemSize-1))

def write_xyz(Ntotal, x, y, z):
	print >> outFile, Ntotal, "\n"
	for i in range(0, Ntotal):
		print >> outFile, "H ",x[i], y[i], z[i]

def lennard_jones(Ntotal, x, y, z,systemSize):
    global Ax, Ay, Az, Fx, Fy, Fz
    Fx = np.zeros(Ntotal)
    Fy = np.zeros(Ntotal)
    Fz = np.zeros(Ntotal)

    Ax = np.zeros(Ntotal)
    Ay = np.zeros(Ntotal)
    Az = np.zeros(Ntotal)

    for i in range(0, Ntotal):
        for j in range(0,Ntotal):
            if (i >= j):        #excludes paris that have already been calculated and the same particles
                continue
            else:
		xvec_prelim = x[j] - x[i]
		if xvec_prelim > systemSize:
			xvec =   xvec_prelim - (systemSize*2)
		elif xvec_prelim <= -systemSize:
                                        xvec = xvec_prelim + (systemSize*2)
		else:
			xvec = xvec_prelim

                yvec_prelim = y[j] - y[i]
                if yvec_prelim > systemSize:
                        yvec =   yvec_prelim - (systemSize*2)
                elif yvec_prelim <= -systemSize:
			yvec = yvec_prelim + (systemSize*2)
		else:
			yvec = yvec_prelim
                zvec_prelim = z[j] - z[i]
                if zvec_prelim > systemSize:
                        zvec =   zvec_prelim - (systemSize*2)
                elif zvec_prelim <= -systemSize:
			zvec = zvec_prelim + (systemSize*2)
		else:
			zvec = zvec_prelim
		r = (xvec**2 + yvec**2 + zvec**2)**0.5

		if (r <= 2.5):
			g = r**2
			U = 4 * eps*(((sigma)**12 / (g)**6) - ((sigma)**6 / (g)**3))
                	Fx[i] = 24 * eps * (( (2 * sigma**12) / g**7 ) - (sigma**6 / g**4) ) * (x[i] - x[j])
                	Fy[i] = 24 * eps * (( (2 * sigma**12) / g**7 ) - (sigma**6 / g**4) ) * (y[i] - y[j])
                	Fz[i] = 24 * eps * (( (2 * sigma**12) / g**7 ) - (sigma**6 / g**4) ) * (z[i] - z[j])

                Ax[i] = Ax[i] + Fx[i]
                Ay[i] = Ay[i] + Fy[i]
                Az[i] = Az[i] + Fz[i]

                Ax[j] = Ax[j] - Fx[i]
                Ay[j] = Ay[j] - Fy[i]
                Az[j] = Az[j] - Fz[i]

def pbc(n, x, y,z,systemSize):
        boxsize = systemSize * 2
        for n in range(0,n):
                if x[n] < (-boxsize*0.5):
                        x[n] = x[n] + boxsize
                if x[n] >= (boxsize*0.5):
                        x[n] = x[n] - boxsize
                if y[n] < (-boxsize*0.5):
                        y[n] = y[n] + boxsize
                if y[n] >= (boxsize*0.5):
                        y[n] = y[n] - boxsize
                if z[n] < (-boxsize*0.5):
                        z[n] = z[n] + boxsize
                if y[n] >= (boxsize*0.5):
                        z[n] = z[n] - boxsize

def steepest_descent(Ax, Ay, Az, Ntotal, x, y, z):
    global Dmax
    Dfunc = (Ax**2 + Ay**2 + Az**2)**0.5
    Dmax = np.amax(Dfunc)
    for i in range(0, Ntotal):
        x[i] = x[i] + (Ax[i]/[Dmax]) * t
        y[i] = y[i] + (Ay[i]/[Dmax]) * t
        z[i] = z[i] + (Az[i]/[Dmax]) * t

pick_coordinates(Ntotal, systemSize)
write_xyz(Ntotal, x, y, z)
lennard_jones(Ntotal, x, y, z,systemSize)

for i in range(0,steps):
    # performs velocity verlet algorithm, step 1
	Ax = Ax + (Fx * 0.5 * timestep)
	Ay = Ay + (Fx * 0.5 * timestep)
	Az = Az + (Fz * 0.5 * timestep)
	pbc(Ntotal, x, y,z,systemSize)
	# performs velocity verlet algorithm, step 2
	x = x + (Ax * timestep)
	y = y + (Ay * timestep)
	z = z + (Az * timestep)
	pbc(Ntotal, x, y,z,systemSize)
	if i % 10 == 0:
		write_xyz(Ntotal, x, y, z)

	lennard_jones(Ntotal, x, y, z,systemSize)

	# performs velocity verlet algorithm, step 3
	Ax = Ax + (Fx * 0.5 * timestep)
	Ay = Ay + (Fx * 0.5 * timestep)
	Az = Az + (Fz * 0.5 * timestep)

#for n in range(0, Ntotal)


outFile.close()
