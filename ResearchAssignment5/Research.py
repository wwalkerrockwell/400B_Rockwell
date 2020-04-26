#What is the time evolution of the inner dark matter density profile of M33?  
#Does it  become  more  or  less  concentrated  with  time?   
#Is  it  well  fit  by  a  Hernquistprofile?  
#What might it mean if there is or isn’t evolution?
import numpy as np
import astropy.units as u
from astropy.constants import G
import matplotlib.pyplot as plt
import matplotlib
from Readfile import Read
from CenterOfMass2 import CenterOfMass2
from MassProfile import MassProfile


#Do what we did with orbitCom, but with MassProfile instead of COM
#Initialize an array to hold all values calculated in each MassProfile snapshot iteration


#Let galaxy be the name 
#Let start-end be the snapshots, and n be the step interval between them
#Let radius by an array of radii to be the profile to be passed in for each snapshot

#Use Local Density rather than Mean Density-----------------------------------------------------------------------------------------------
def DensityProfileTime(galaxy, start, end, n, radius):
    fileout1 = "M33_Density.txt"
    fileout2 = "M33_Mass.txt"
    snapID = np.arange(start, end, n) #Array for deciding files to be read in
    density = np.zeros([snapID.size, radius.size + 1]) #Array for containing all Profiles
    Mhalo = np.zeros([snapID.size, radius.size + 1])

    for i, snapID in enumerate(snapID):
        Profile = MassProfile(galaxy, snapID)
        Mhalo[i, 1:radius.size+1], time = Profile.MassEnclosed(1, radius) #Initialize Mass array to be plugged into density
        print("i = ", i) 
        Mhalo[i, 0] = time
        density[i, 0] = time

        for j in range(np.size(radius)):       
            if (j != 0):
                density[i, j + 1] = Mhalo[i, j + 1] / ((4./3.)*3.141592*((radius[j]**3) - (radius[j - 1]**3))) * 10e2 #Calculate local density
            else:
                density[i, j + 1] = Mhalo[i, j + 1] / ((4./3.)*3.141592*(radius[j]**3)) * 10e2 #Calculate local density for first index

        print("Mhalo = ", Mhalo[i])
        print("Density = ", density[i])

    np.savetxt(fileout1, density, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
    np.savetxt(fileout2, Mhalo, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
    return
#Each row of the file represents a snapshot; what would normally be a 1d output array of masses at varying radii
#This means that each column represents a different radius, save for the first which is time
#The plan is to plot out each radius as a line, evolving over time. 
#So radius 1 will be a line of Time vs Density values over that column, and so on for each other column.
#We could also just make some simple radius vs density lines

#Make a function that will run HernquistMass for all the snapshots
#It should compare each row 

def HernquistMass(scale, Mhalo, R):
    # Determine the mass enclosed using Hernquist 1990 Mass profile 
    # Input:   R   Radius  
    #         scale   Scale Length  
    #         Mhalo  Total Halo Mass (Msun)
    # Returns: Mass in units Msun. 
    
        # Hernquist 1990 Mass profile
    M = Mhalo*R**2/(R+scale)**2
    Density = (M * scale) / ((2 + 3.1415 * R) * (r + scale)**3)
    return Density

def HernquistDensityProfile(a, Mhalo, radius, start, end, n):
    snapID = np.arange(start, end, n)
    HernquistDensity = np.zeros([snapID.size, radius.size + 1])
    print(HernquistDensity)
    
    for i in range(np.size(Mhalo)):
        print(Mhalo[i][0])
        print(Mhalo[i][radius.size])
        HernquistDensity[i, 0] = Mhalo[i][0]
        HernquistDensity[i, 1:radius.size+1] = HernquistMass(a, Mhalo[i][1:radius.size + 1], radius)

    np.savetxt("M33HernquistDensity", density, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
    return

radius = np.arange(1, 32, 5)
a = 1 #Scale Factor
#DensityProfileTime("M33", 0, 801, 1, radius)

Mhalo = np.genfromtxt("M33_Mass_Local.txt", dtype = None, names = True, skip_header = 0) #Create double array for all masses
HernquistDensityProfile(a, Mhalo, radius, 0, 10, 1)
print(radius)


"""
plt.plot(Mhalo["Time"], Mhalo["Second"], label="R = 6kpc")
plt.plot(Mhalo["Time"], Mhalo["Third"], label="R = 11kpc")
plt.plot(Mhalo["Time"], Mhalo["Fourth"], label="R = 16kpc")
plt.plot(Mhalo["Time"], Mhalo["Fifth"], label="R = 21kpc")
plt.plot(Mhalo["Time"], Mhalo["Sixth"], label="R = 26kpc")
plt.plot(Mhalo["Time"], Mhalo["Seventh"], label="R = 31kpc")
plt.legend()

plt.xlabel("Time (Myr)")
plt.ylabel("Density (10^2 Msun / kpc^3)") 

plt.show()
"""



