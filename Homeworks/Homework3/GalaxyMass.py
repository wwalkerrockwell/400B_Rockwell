import numpy as np
import astropy as u
from astropy import units as g
from Readfile import Read

def ComponentMass(filename, particleType): #Function for summing up mass of a particular type in a galaxy
    data, time, particles = Read(filename) #Read in the data from the given file

    i = 0 #initialize totalMass and loop index values
    totalMass = 0

    #data.size represents the total number of data indexes in the file
    while i < data.size: #While the index is of the desired type, sum up the masses of all values in dataset
        if (data['type'][i] == particleType):
            totalMass += data['m'][i]
        i += 1

    totalMass *= 1e2 * g.Msun #Mass is in 1e10, so multiply by 1e2 and g.Msun to make units proper
    totalMass = np.around(totalMass, 3) #Round

    return totalMass; #Return the desired mass

file1 = "M31_000.txt" #initialize each of the the desired datasets
file2 = "M33_000.txt"
file3 = "MW_000.txt"

Mass = ComponentMass(file1, 1)
print(Mass)
Mass = ComponentMass(file1, 2) 
print(Mass)
Mass = ComponentMass(file1, 3)
print(Mass)
Mass = ComponentMass(file2, 1)
print(Mass)
Mass = ComponentMass(file2, 2) 
print(Mass)
Mass = ComponentMass(file2, 3)
print(Mass)
Mass = ComponentMass(file3, 1)
print(Mass)
Mass = ComponentMass(file3, 2) 
print(Mass)
Mass = ComponentMass(file3, 3)
print(Mass)
