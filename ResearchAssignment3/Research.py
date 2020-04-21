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
from CenterOfMass2 import CenterOfMass2 #Use COM2
from MassProfile import MassProfile


#Do what we did with orbitCom, but with MassProfile instead of COM
#Initialize an array to hold all values calculated in each MassProfile snapshot iteration


#Let galaxy be the name 
#Let start-end be the snapshots, and n be the step interval between them
#Let radius by an array of radii to be the profile to be passed in for each snapshot
def DensityProfileTime(galaxy, start, end, n, radius):
    fileout1 = "M33_Density.txt"
    fileout2 = "M33_Mass.txt"
    snapID = np.arange(start, end, n) #Array for deciding files to be read in
    density = np.zeros([snapID.size, radius.size]) #Array for containing all Profiles
    Mhalo = np.zeros([snapID.size, radius.size])

    for i, snapID in enumerate(snapID):
        Profile = MassProfile(galaxy, snapID)
        Mhalo[i] = Profile.MassEnclosed(1, radius) #Initialize Mass array to be plugged into density
        print("i = ", i)
        for j in range(np.size(radius)):
            density[i, j] = Mhalo[i, j] / (4/3)*(radius[j]**3) #Calculate density

        print("Mhalo = ", Mhalo[i])
        print("Density = ", density[i])
    print(density)
    np.savetxt(fileout1, density, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
    np.savetxt(fileout2, Mhalo, fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
    return
#Each row of the file represents a snapshot; what would normally be a 1d output array of masses at varying radii
#This means that each column represents a different radius
#The plan is to plot out each radius as a line, evolving over time. 
#So radius 1 will be a line of Time vs Density values over that column, and so on for each other column.
#We could also just make some simple radius vs density lines

#Make a function that will run HernquistMass for all the snapshots
#It should compare each row 

def HernquistMass(a, Mhalo, radius):
    M = np.zeros(radius.size) #initialize arrays
    Density = np.zeros(radius.size)

    for i in range(np.size(radius)):
        M[i] = (Mhalo[i] * (radius[i]**2)) / (a + radius[i])**2 #Mass calculation
        Density[i] = ((M[i]*a) / ((2 * 3.1415 * radius[i]) * (radius[i] + a)**3)) #Density Equation

    print("Hernquist Mass = ", M)
    return Density

def HernquistProfile(a, Mhalo, radius):

    for i in range(np.size(Mhalo[i])):
        HernquistMass(a, Mhalo[i], radius)


    return

radius = np.arange(1, 32, 5)
a = 1 #Scale Factor
#DensityProfileTime("M33", 0, 801, 1, radius)

Mhalo = np.genfromtxt("M33_Density3.75.txt", dtype = None, names = True, skip_header = 0) #Create double array for all masses
time = np.arange(0, 801)

plt.plot(time, Mhalo["Third"], label="R = 11")
plt.plot(time, Mhalo["Fourth"], label="R = 16")
plt.plot(time, Mhalo["Fifth"], label="R = 21")
plt.plot(time, Mhalo["Sixth"], label="R = 26")
plt.legend()

plt.xlabel("Snapshot")
plt.ylabel("Density") 

plt.show()



