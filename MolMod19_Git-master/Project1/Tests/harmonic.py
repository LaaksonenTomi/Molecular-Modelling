import math

#constants
timestep = 0.001
steps = 1000
mass1 = 1
mass2 = 1
massCM = (mass1 * mass2)/(mass1 + mass2)
spring_constant = 1
eq_distance = (r2 - r1)/linis.norm((r2 - r1)) * r_eq

#starting values
velx1 = 0
vely1 = 0
velx2 = 0
vely2 = 0
x1 = 0.75
y1 = 0
x2 = -0.75
y2 = 0
r = r2 - r1
theta = math.atan((y2 - y1)/(x2 - x1))
if (velx1**2 + vely1**2) == 0:
    alpha1 = 0
    alpha2 = 0
else:
    alpha1 = math.acos(((x2 -x1)*velx1 + (y2 - y1)*vely1)/(math.sqrt((x2 - x1)**2 + (y2 - y1)**2)*math.sqrt(velx1**2 + vely1**2)))
    alpha2 = math.acos(((x2 -x1)*velx2 + (y2 - y1)*vely2)/(math.sqrt((x2 - x1)**2 + (y2 - y1)**2)*math.sqrt(velx2**2 + vely2**2)))
vel1_par = math.cos(alpha1) * math.sqrt(velx1**2 + vely1**2)
vel1_per = math.sin(alpha1) * math.sqrt(velx1**2 + vely1**2)
vel2_par = math.cos(alpha2) * math.sqrt(velx1**2 + vely1**2)
vel2_per = math.sin(alpha2) * math.sqrt(velx1**2 + vely1**2)
#r = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
ang_velocity1 = math.sqrt(velx1**2 + vely1**2) / (r/2)
ang_velocity2 = math.sqrt(velx2**2 + vely2**2) / (r/2)
K1 = 1/2 * mass1 * math.sqrt(velx1 ** 2 + vely1 **2)
K2 = 1/2 * mass2 * math.sqrt(velx2 ** 2 + vely2 **2)
U = (r - eq_distance) * spring_constant
E_rot = math.sqrt(ang_velocity2**2 + ang_velocity1**2) / (4 * r)

#opening a file
#3 ja hiili kuvaa CM:C$C$, helppo poistaa
f = open("harmonic.txt", "a")
#defining writing step
def writestep():
    f.write("2\nstep: " + str(i) + " K1: " + str(K1) + " K2: " + str(K2) + " U: " + str(U) + " E_J: " + str(E_rot) + " E_tot: " + str(K1 + K2 + U + E_rot) + "\nH "+ str(x1) + " " + str(y1) + " 0\nH "+ str(x2) + " " + str(y2) + " 0\n")
  #  f.write(str(velx1) + " vely1: " + str(vely1) + " velx2: " + str(velx2) + " r: " + str(r)+ " theta: " + str(theta) + " alpha1: " + str(alpha1) + " alpha2" + str(alpha2) + " vel1_per: " + str(vel1_per) + " vel1_par: " + str(vel1_par) + " vel2_per: " + str(vel2_per) + " vel2_par: " + str(vel2_par) + " ang1: " + str(ang_velocity1) + " ang2: " + str(ang_velocity2) + "\nH")
#leapfrog (?)
i = 0
writestep()

for i in range(steps):
    
    x1 += velx1 * timestep
    y1 += vely1 * timestep
    x2 += velx2 * timestep
    y2 += vely2 * timestep
    r = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    #ang_velocity1 = (ang_velocity1 + vel1_per)/ (r/2)
    #ang_velocity2 = (ang_velocity2 + vel2_per)/ (r/2)
    theta += vel1_per * timestep
    # r 1 -> 2
    vel1_par -= (spring_constant/mass1) * (r - eq_distance)
    velx1 = math.cos(alpha1)*math.cos(theta)*vel1_par + math.sin(alpha1)*math.cos(theta)*vel1_per
    vely1 = math.cos(alpha1)*math.sin(theta)*vel1_par + math.sin(alpha1)*math.sin(theta)*vel1_per
    vel2_par += (spring_constant/mass2) * (r - eq_distance)
    velx2 = math.cos(alpha2)*math.cos(theta)*vel2_par + math.sin(alpha2)*math.cos(theta)*vel2_per
    vely2 = math.cos(alpha2)*math.sin(theta)*vel2_par + math.sin(alpha2)*math.sin(theta)*vel2_per
    p1 = mass1 * math.sqrt(vel1_par**2+ vel1_per**2)
    p2 = mass2 * math.sqrt(vel2_par**2+ vel2_per**2)
 
    #energy 
    U = (r - eq_distance) * spring_constant
    K1 = p1**2 /(2 *mass1)
    K2 = p2**2 /(2 *mass2)
    linear_momentum = p1 + p2
    angular_momentum = massCM * math.sqrt(ang_velocity1**2 + ang_velocity2**2) * r/2
    E_rot = math.sqrt(ang_velocity2**2 + ang_velocity1**2) / (4 * r)
    E_tot = U + K1 + K2 + E_rot
    writestep()

f.close
