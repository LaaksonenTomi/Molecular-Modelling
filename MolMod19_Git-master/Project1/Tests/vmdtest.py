import numpy
import random
import math

atoms = 100
steps = 50
matrix = numpy.ones((atoms, 2), dtype = numpy.float32)
       
for i in range(atoms):
    matrix[i][0] = 0
    matrix[i][1] = 0

def timestep():
    for i in range(atoms):
        a = random.random()
        matrix[i][0] += math.cos(a * 2 * math.pi)
        matrix[i][1] += math.sin(a * 2 * math.pi)

f = open("xyz.txt", "a")
f.write(str(atoms) + '\n')
for g in range(atoms):
    f.write('\n' + str(g + 1) + '\t0\t0\t0')
f.write('\n')
for i in range(steps):
    timestep()
    f.write(str(atoms) + '\n')
    for h in range(atoms):
        d = round((matrix[h][0]), 3)
        e = round((matrix[h][1]), 3)
        f.write('\n' + str(h + 1) + '\t' + str(d) + '\t' + str(e) + '\t 0')
    f.write('\n')    
