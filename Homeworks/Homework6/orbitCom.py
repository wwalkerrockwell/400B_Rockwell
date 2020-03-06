import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from Readfile import Read
from CenterOfMass import CenterOfMass


def OrbitCom(galaxy, start, end, n):
    
    fileout = "Orbit_M31.txt" #initialize file to contain all values
    
    delta = 0.1 #initialize delta and VolDec
    VolDec = 2
    
    snapID = np.arange(start, end, n) #Array for deciding files to be read in
    orbit = np.zeros([snapID.size, 7]) #Array for containing all COM's
    
    
    for i, snapID in enumerate(snapID):
        ilbl = '000' + str(snapID) #Initialize new file to be read in from desired galaxy
        ilbl = ilbl[-3:]
        filename = "%s_"%(galaxy) + ilbl + ".txt"

        print(snapID)
        print(i)
        print(filename)
        
        COM = CenterOfMass(filename, 2) #Take in COM values
        orbitCOMP = COM.COM_P(VolDec, delta) #COM Position
        orbitCOMV = COM.COM_V(orbitCOMP[0], orbitCOMP[1], orbitCOMP[2]) #COM Velocity

        orbit[i, 0] = COM.time.value #assign time and COM values in order.
        orbit[i, 1] = orbitCOMP[0].value
        orbit[i, 2] = orbitCOMP[1].value
        orbit[i, 3] = orbitCOMP[2].value
        orbit[i, 4] = orbitCOMV[0].value
        orbit[i, 5] = orbitCOMV[1].value
        orbit[i, 6] = orbitCOMV[2].value

    #Output orbit array to file
    np.savetxt(fileout, orbit, header="time  x  y  z  vx  vy  vz", comments='#', fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])

    return orbit

def MagDiff(Vector1, Vector2, index1, index2, index3): #Calculate Magnitude difference of desired vectors
    magnitude = np.arange(np.size(Vector1)) #Initialize array for magnitdue calculation

    for i in range(np.size(magnitude)): #For loop through all indexes to calculate Magnitude for each index
        magnitude[i] = np.sqrt((Vector1[i][index1] - Vector2[i][index1])**2 + (Vector1[i][index2] - Vector2[i][index2])**2 +(Vector1[i][index3] - Vector2[i][index3])**2) 
    return magnitude



M33 = np.genfromtxt("Orbit_M33.txt", dtype = None, names = True, skip_header = 0)
MW = np.genfromtxt("Orbit_MW.txt", dtype = None, names = True, skip_header = 0)
M31 = np.genfromtxt("Orbit_M31.txt", dtype = None, names = True, skip_header = 0)

Diff = MagDiff(M33, M31, 1, 2, 3) #Enter 1-3 for Position values, 4-6 for Velocity values

plt.plot(M31['time'], Diff)
plt.legend(['M31-M33'])

plt.xlabel("Time (Myr)")
plt.ylabel("Distance (kpc)") 
plt.show()


#1. 3 Close encounters will be experiences, with the third one resulting in the final fusion
#2. Spikes in velocity represent close encounters due to gravitational pull experienced by both galaxies
#3. In about 6.1 Gyrs they will merge. M33 slows down exponentially upon MW and M31's merger, 
#and it will begin to orbit the merged galaxies more closely before stabilizing.
