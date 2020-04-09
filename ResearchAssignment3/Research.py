#What is the time evolution of the inner dark matter density profile of M33?  
#Does it  become  more  or  less  concentrated  with  time?   
#Is  it  well  fit  by  a  Hernquistprofile?  
#What might it mean if there is or isn’t evolution?
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from Readfile import Read
from CenterOfMass import CenterOfMass


#Do what we did with orbitCom, but with MassProfile instead of COM
#Initialize an array to hold all values calculated in each MassProfile snapshot iteration


#Let galaxy be the name 
#Let start-end be the snapshots, and n be the step interval between them
#Let radius by an array of radii to be the profile to be passed in for each snapshot
def DensityProfileTime(galaxy, start, end, n, radius):


    snapID = np.arange(start, end, n) #Array for deciding files to be read in
    density = np.zeros([snapID.size, radius.size + 1]) #Array for containing all Profiles

    for i, snapID in enumerate(snapID):
        Profile = MassProfile(galaxy, snapID)
        Mbulge = Profile.MassEnclosed(3, radius) #Initialize Mass array to be plugged into density
        density[i, 0] = i #Initialize column of ID's

        for j in range(np.size(radius)):
            density[i, j + 1] = Mbulge[j] / (radius[j]**3) #Calculate density

        #Save to a file
    return

#Make a function that will run HernquistMass for all the snapshots





