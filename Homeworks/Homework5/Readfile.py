import numpy as np
import astropy.units as u
def Read(filename):
    file = open(filename, "r") #Open file

    line1 = file.readline() #Read first line
    label, value = line1.split() #Split up the string and float into respective 
    time = float(value) * u.Myr #Return time 

    line2 = file.readline() #Repeat for second line
    label, value = line2.split()
    particles = float(value) #Total Number of Particles

    file.close() #close file
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3) #Create array for all other data values
    #print(data['x'][2])
    #print(data)
    
    return data, time, particles;
filename = "MW_000.txt" #Declare filename for repeated use

Read(filename) #run function
