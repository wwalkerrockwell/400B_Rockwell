import numpy as np
import math
import astropy as u
from astropy import units as g
from Readfile import Read
def ParticleInfo(filename, type, number): 

    data = Read(filename) #Read in data from file
    Dist_Mag = math.sqrt((data['x'][number])**2 + (data['y'][number])**2 + (data['z'][number])**2) #Calculate distance magnitude
    Vel_Mag = math.sqrt((data['vx'][number])**2 + (data['vy'][number])**2 + (data['vz'][number])**2) #Calculate Velocity magnitude

    Dist_Mag = np.around(Dist_Mag, 3) #Round distance and velocity
    Vel_Mag = np.around(Vel_Mag, 3)

    Dist_Mag *= g.kpc #Convert Distance to kpc

    print("Distance Magnitude = ", Dist_Mag) #Print out desired values
    print("Velocity Magnitude  = ", Vel_Mag)
    print("Mass = ", data['m'][number])

    Dist_Mag = Dist_Mag.to(g.lyr) #Convert to light years
    Dist_Mag = np.around(Dist_Mag, 3) #Round
    print("Distance Magnitude = ", Dist_Mag)


    return;
filename = "MW_000.txt"
ParticleInfo(filename, 1, 99)

