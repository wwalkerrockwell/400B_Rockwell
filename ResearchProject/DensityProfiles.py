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

#Let galaxy be the name 
#Let start-end be the snapshots, and n be the step interval between them
#Let radius by an array of radii to be the profile to be passed in for each snapshot
def DensityProfile(galaxy, start, end, n, radius):
    fileout1 = "M33_Density.txt"
    fileout2 = "M33_Mass.txt"
    snapID = np.arange(start, end, n) #Array for deciding files to be read in
    density = np.zeros([snapID.size, radius.size + 1]) #Array for containing all Profiles
    Mhalo = np.zeros([snapID.size, radius.size + 1])

    for i, snapID in enumerate(snapID): #Loop through snapshots and calculate their density profiles
        Profile = MassProfile(galaxy, snapID)
        Mhalo[i, 1:radius.size+1], time = Profile.MassEnclosed(1, radius) #Initialize Mass array to be plugged into density
        print("i = ", i) 
        Mhalo[i, 0] = time #initialize time
        density[i, 0] = time

        for j in range(np.size(radius)): #Loop through radii and calculate their local density values
            if (j != 0):
                density[i, j + 1] = Mhalo[i, j + 1] / ((4./3.)*3.141592*((radius[j]**3) - (radius[j - 1]**3))) * 10e2 #Calculate local density
            else:
                density[i, j + 1] = Mhalo[i, j + 1] / ((4./3.)*3.141592*(radius[j]**3)) * 10e2 #Calculate local density for first index

    np.savetxt(fileout1, density, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f']) #Save density and mass arrays to files
    np.savetxt(fileout2, Mhalo, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
    return
#Each row of the file represents a snapshot; what would normally be a 1d output array of masses at varying radii
#This means that each column represents a different radius, save for the first which is time
#Plot out each radius as a line, evolving over time. 

#a is scale length, Mhalo is the halo mass, R is the radius.
def HernquistMass(a, Mhalo, R):
    M = Mhalo*R**2/(R+a)**2
    Density = (M * a) / ((2 + 3.1415 * R) * (R + a)**3) #Density profile
    return Density

def HernquistDensityProfile(a, Mhalo, radius, start, end, n): #Calculate Hernquist density for each snapshot
    snapID = np.arange(start, end, n)
    HernquistDensity = np.zeros([snapID.size, radius.size + 1])

    for i in range(np.size(snapID)): #Loop through indices of the density profiles to calculate their hernquist counterparts
        HernquistDensity[i, 0] = Mhalo[i][0] #Initialize time
        HernquistDensity[i, 1:radius.size+1] = HernquistMass(a, Mhalo[i][1:radius.size+1], radius) * 10e4 

    np.savetxt("M33HernquistDensity.txt", HernquistDensity, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f']) #Save to file
    return

#This function exists so I could look at the mass profiles as well to ensure consistency
#It uses the same inputs as HernquistDensityProfile
#Simply change the HernquistMass function to return mass instead of density to use this
def HernquistMassProfile(a, Mhalo, radius, start, end, n): 
    snapID = np.arange(start, end, n)
    Mass = np.zeros([snapID.size, radius.size + 1])

    for i in range(np.size(snapID)):
        Mass[i, 0] = Mhalo[i][0]
        Mass[i, 1:radius.size+1] = HernquistMass(a, Mhalo[i][1:radius.size+1], radius)

    np.savetxt("M33HernquistMass.txt", Mass, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f']) #Save to file
    return


radius = np.arange(1, 32, 5)
a = 1 #Scale Factor
#DensityProfile("M33", 0, 801, 1, radius)

#Mhalo = np.loadtxt("M33_Mass.txt", skiprows=1)
#HernquistDensityProfile(a, Mhalo, radius, 0, 801, 1)
print(radius)


Density = np.genfromtxt("M33_Density_Local.txt", dtype = None, names = True, skip_header = 0) #PLOT DENSITY PROFILES  OVER TIME
plt.plot(Density["Time"], Density["Second"] * 10e-2, label="R = 6kpc") 
plt.plot(Density["Time"], Density["Third"] * 10e-2, label="R = 11kpc")
plt.plot(Density["Time"], Density["Fourth"] * 10e-2, label="R = 16kpc")
plt.plot(Density["Time"], Density["Fifth"] * 10e-2, label="R = 21kpc")
plt.plot(Density["Time"], Density["Sixth"] * 10e-2, label="R = 26kpc")
plt.plot(Density["Time"], Density["Seventh"] * 10e-2, label="R = 31kpc")




#HernquistDensity = np.genfromtxt("M33HernquistDensity.txt", dtype = None, names = True, skip_header = 0) #PLOT HERNQUIST PROFILES OVER TIME
#plt.plot(HernquistDensity["Time"], HernquistDensity["Second"] * 10e-4 * 7, label="R = Hernquist6kpc") 
#plt.plot(HernquistDensity["Time"], HernquistDensity["Third"] * 10e-4 * 9, label="R = Hernquist11kpc")
#plt.plot(HernquistDensity["Time"], HernquistDensity["Fourth"], label="R = Hernquist16kpc")
#plt.plot(HernquistDensity["Time"], HernquistDensity["Fifth"], label="R = Hernquist21kpc")

plt.legend()
plt.xlabel("Time (Myr)")
plt.ylabel("Density (Msun / kpc^3)") 
plt.show()


"""
Mhalo = np.loadtxt("M33_Mass.txt", skiprows=1) #Plot Hernquist Mass
HernquistMass = np.loadtxt("M33HernquistMass.txt", skiprows=1)
#HernquistMassProfile(a, Mhalo, radius, 0, 801, 1)
plt.plot(radius, Mhalo[800][1:radius.size+1], label="Mass")
plt.plot(radius, HernquistMass[800][1:radius.size+1], label="HernquistMass")

plt.ylabel("Mass (Msun)")
plt.xlabel("Radius (kpc)")
plt.semilogy()
plt.legend()
plt.show()
"""


snapshot = 000#Plot critical Density Snapshots with Hernquist Profiles. Change snapshot value to see different snapshots
Density = np.loadtxt("M33_Density_Local.txt", skiprows=1) 
HernquistDensity = np.loadtxt("M33HernquistDensity.txt", skiprows=1)
plt.plot(radius, Density[snapshot][1:radius.size+1] * 10e-2, label="Density")
plt.plot(radius, HernquistDensity[snapshot][1:radius.size+1] * 10e-4, label="HernquistDensity")

plt.ylabel("Density ( Msun / kpc^3)")
plt.xlabel("Radius (kpc)")
plt.legend()
plt.semilogy()
plt.show()

