# Homework 4
# Center of Mass Position and Velocity
# William Rockwell

###############################################################################
# Keep in mind this is just a template; you don't need to follow every step and feel free to change anything.
# We also strongly encourage you to try to develop your own method to solve the homework.
###############################################################################

# import modules
import numpy as np
import math
import astropy.units as u
from Readfile import Read


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.data, self.time, self.total  = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

        


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)
        
        # write your own code to compute the generic COM using Eq. 1 in the homework instructions

        i = 0
        totalMass = 0
        for i in m: #Sum up total mass
            totalMass += i

        xSum = 0 #Initialize variables for summation
        ySum = 0
        zSum = 0
        i = 0
        j = 0
        k = 0

        for i, val in enumerate(a): #Sum up numerators of center of Mass equations for x, y, and z
            xSum += a[i] * m[i]

        for j, val2 in enumerate(b):
            ySum += b[j] * m[j]

        for k, val3 in enumerate(c):
            zSum += c[k] * m[k]

        # xcomponent Center of mass
        Acom = xSum / totalMass
        # ycomponent Center of mass
        Bcom = ySum / totalMass
        # zcomponent Center of mass
        Ccom = zSum / totalMass
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    
        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below

        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)


        # iterative process to determine the center of mass                                                            
        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)
        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW < RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            #print ("maxR", maxR)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            xNew = x2 - XCOM2
            yNew = y2 - YCOM2
            zNew = z2 - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        COMP *= u.kpc
        np.around(COMP, 3)
        return COMP
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(self.vx, self.vy, self.vz, self.m)
        xV = self.vx - VXCOM
        yV = self.vy - VYCOM
        zV = self.vz - VZCOM
        RV = np.sqrt(xV**2 + yV**2 + zV**2) * u.kpc
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV < RVMAX)

        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV] 
        vznew = self.vz[indexV]
        mnew = self.m[indexV]
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # write your own code below
        COMV = [VXCOM, VYCOM, VZCOM]
        COMV *= u.km/u.s
        np.around(COMV, 3)

        # return the COM vector                                                                                        
        return COMV
    


# ANSWERING QUESTIONS
#######################


# Create  a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("M33_000.txt", 2)
MWCOM1 = CenterOfMass("M31_000.txt", 2)

# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM 
MW_COMP = MWCOM.COM_P(0.1)
MW_COMV = MWCOM.COM_V(MW_COMP[0],MW_COMP[1],MW_COMP[2])

MW_COMP1 = MWCOM1.COM_P(0.1)
MW_COMV1 = MWCOM1.COM_V(MW_COMP1[0],MW_COMP1[1],MW_COMP1[2])

COMP_Mag = np.sqrt(MW_COMP[0]**2 + MW_COMP[1]**2 + MW_COMP[2]**2)
COMV_Mag = np.sqrt(MW_COMV[0]**2 + MW_COMV[1]**2 + MW_COMV[2]**2)

sepMag = np.sqrt((MW_COMP[0] - MW_COMP1[0])**2 + (MW_COMP[1] - MW_COMP1[1])**2 + (MW_COMP[2] - MW_COMP1[2])**2) #Calculate separation magnitude
VelSepMag = np.sqrt((MW_COMV[0] - MW_COMV1[0])**2 + (MW_COMV[1] - MW_COMV1[1])**2 + (MW_COMV[2] - MW_COMV1[2])**2) #Calculate separation magnitude

print("Position Magnitude = ", COMP_Mag)
print("Velocity Magnitude = ", COMV_Mag)
print("Separation Magnitude = ", sepMag)
print("Velocity Separation Magnitude = ", VelSepMag)

# now write your own code to answer questions
#1
#MW
#Com_P = [-0.874  2.394 -1.422 ] kpc
#Position Magnitude = 2.9193 kpc
#Com_V = [-2.524  5.067  1.512] km / s
#Velocity Magnitude = 5.859 km / s

#M31
#Com_P = [-476.220  491.440  -412.402] kpc
#Position Magnitude = 798.9833 kpc
#Com_V = [ 45.518 101.956 141.343] km / s
#Velocity Magnitude = 180.1248 km / s

#M33
#Com_P = [-377.658  611.429 -284.637] kpc
#Position Magnitude = 772.974 kpc
#Com_V = [ 75.489 -72.208  47.322  ] km / s
#Velocity Magnitude =  114.682 km / s

#2
#Distance Separation between MW and M31 = 770.1293 kpc
#Velocity Separation between MW and M31 = 118.9794 km / s

#3
#Distance Separation between M33 and M31 = 201.0865 kpc
#Velocity Separation between M33 and M31 = 200.1788 km / s

#4
#The critical part of the collision is the center of mass and velocity, as those affect the dynamic
#of the collision between galaxies the most