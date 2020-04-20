import numpy as np
import math
import astropy.units as u
from astropy.constants import G
from Readfile import Read
from array import *
import matplotlib.pyplot as plt
from CenterOfMass2 import CenterOfMass2

class MassProfile:

    def __init__(self, galaxy, snap):
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + ".txt"
        self.delta = 0.1 #initialize delta and VolDec
        self.VolDec = 4

        
        self.data, self.time, self.total  = Read(self.filename) #Read in the file and assign desired values
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc
        self.m = self.data['m']
        self.gname = galaxy
        #print(self.gname)
        self.Gravity = G 
        self.Gravity = self.Gravity.to(u.kpc*u.km**2/u.s**2/u.Msun) #convert Gravity to desired units


  
    def MassEnclosed(self, type, radius):
        print(self.filename)
        print("type = ", type)
        COM = CenterOfMass2(self.filename, type) #calculate center of Mass
        galCOMP = COM.COM_P(self.VolDec, self.delta)
        
        
        index = np.where(self.data['type'] == type) #Select desired indexes and adjust values
        xIndex = self.x[index] - galCOMP[0]
        yIndex = self.y[index] - galCOMP[1]
        zIndex = self.z[index] - galCOMP[2]

        RIndex = np.sqrt(xIndex**2 + yIndex**2 + zIndex**2) #Calculate the distance magnitudes
        mIndex = np.zeros(radius.size) #initialize an array for storing all mass profiles
        for i in range(np.size(radius)): #loop through each radius and sum up each total mass profile value
            index2 = np.where(radius[i] > RIndex.value)
            mIndex[i] = np.sum(self.m[index][index2])

        mIndex *= u.Msun
        #print(mIndex)

        return mIndex 

    def MassEnclosedTotal(self, radius):
        rows, cols = (3, np.size(radius)) #initialize arrays for total Mass calculation
        Mass = [[0] * cols] * rows
        totalMass = np.zeros(radius.size) * u.Msun
        i = 0

        while i < 3: #Calculate the 3 masses
            if self.gname == 'M33' and i == 2: #Contingency for galaxy with no bulge
                break
            Mass[i] = self.MassEnclosed(i + 1, radius) 
            totalMass += Mass[i]
            i += 1

        #print("TotalMass = ", Mass)
        #print("TotalMassValue = ", totalMass)
        return totalMass

    """
    def HernquistMass(self, scale, Mhalo, R):
    # Determine the mass enclosed using Hernquist 1990 Mass profile 
    # Input:   R   Radius  
    #         scale   Scale Length  
    #         Mhalo  Total Halo Mass (Msun)
    # Returns: Mass in units Msun. 
    
        # Hernquist 1990 Mass profile
        print(Mhalo*R**2/(R+scale)**2*)
        return Mhalo*R**2/(R+scale)**2

    """

    def HernquistMass(self, a, Mhalo, radius):
        M = np.zeros(radius.size) * u.Msun #initialize arrays
        Density = np.zeros(radius.size) * u.Msun / u.kpc**3

        for i in range(np.size(radius)):
            M[i] = (Mhalo[i] * (radius[i]**2)) / (a + radius[i])**2 #Mass calculation
            Density[i] = ((M[i]*a) / ((2 * 3.1415 * radius[i]) * (radius[i] + a)**3)) #Density Equation

        #print("Hernquist Mass = ", M)
        return Density
  

    def CircularVelocity(self, type, radius):
        Mass = self.MassEnclosed(type, radius) #take in mass values
        Velocity = np.zeros(radius.size) * u.km / u.s #initialize velocity array
        for i in range(np.size(radius)): #Calculate velocity profiles
            Velocity[i] = np.sqrt(self.Gravity*Mass[i] / radius[i])


        #print(Velocity)
        return Velocity

    def CircularVelocityTotal(self, radius):
        rows, cols = (3, np.size(radius)) 
        Mass = self.MassEnclosedTotal(radius) #Calculate total mass for each type
        Velocity = np.zeros(radius.size) * u.km / u.s
        for i in range(np.size(radius)): #Repeat Circular Velocity but with total mass
            Velocity[i] = np.sqrt(self.Gravity*Mass[i] / radius[i])

        #print(Velocity)
        return Velocity

    def HernquistVCirc(self, a, Mhalo, radius): #just do the same thing as in Circular velocity but with Hernquist Mass function
        Mass = self.HernquistMass(a, Mhalo, radius)
        Velocity = np.zeros(radius.size) * u.km / u.s #initialize velocity array
        for i in range(np.size(radius)): #Calculate velocity profiles
            Velocity[i] = np.sqrt(self.Gravity*Mass[i] / radius[i])


        #print(Velocity)
        return Velocity

#Profile = MassProfile("M31", 0)

#radius = np.arange(1, 30) * u.kpc
#a = 1 * u.kpc #Scale Factor
#MHalo = Profile.MassEnclosed(1, radius)
#Mdisk = Profile.MassEnclosed(2, radius)
#Mbulge = Profile.MassEnclosed(3, radius)
#TotalMass = Profile.MassEnclosedTotal(radius)
#Hernquist = Profile.HernquistMass(a, MHalo, radius)