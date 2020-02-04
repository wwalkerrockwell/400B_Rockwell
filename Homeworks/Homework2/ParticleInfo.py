import numpy as np
import math
import astropy as u
from astropy import units as g
from Readfile import Read
def ParticleInfo(filename, type, number): 

    data, time, total = Read(filename) #Read in data from file
    Dist_Mag = math.sqrt((data['x'][number])**2 + (data['y'][number])**2 + (data['z'][number])**2) #Calculate distance magnitude
    Vel_Mag = math.sqrt((data['vx'][number])**2 + (data['vy'][number])**2 + (data['vz'][number])**2) #Calculate Velocity magnitude
    mass = data['m'][number]

    Dist_Mag *= g.kpc #Convert Distance to kpc
    Vel_Mag *= g.km/g.s #Convert Velocity to km/s
    Dist_Mag = np.around(Dist_Mag, 3) #Round
    Vel_Mag = np.around(Vel_Mag, 3)

    index = np.where(data['type'] == type) #catalogue values of desired type
    mnew = data['m'][index] * 1e10 * g.Msun
    xnew = data['x'][index] * g.kpc
    ynew = data['y'][index] * g.kpc
    znew = data['z'][index] * g.kpc
    vxnew = data['vx'][index] * g.km/g.s
    vynew = data['vy'][index] * g.km/g.s
    vznew = data['vz'][index] * g.km/g.s

    return Dist_Mag, Vel_Mag, mass;
filename = "MW_000.txt"
type = int(input("Enter particle type: "))
index = int(input("Enter particle number: "))
distance, velocity, mass = ParticleInfo(filename, type, index)

distance = distance.to(g.lyr) #Convert to light years
print("Distance Magnitude = ", distance)
print("Velocity Magnitude  = ", velocity)

