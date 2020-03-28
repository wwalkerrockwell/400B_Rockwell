import numpy as np
import astropy.units as u
from astropy.constants import G
import matplotlib.pyplot as plt
import matplotlib
from Readfile import Read
from CenterOfMass import CenterOfMass


class M33AnalyticOrbit:

    def __init__(self, filename):
        self.filename = filename
        M33COM = CenterOfMass("M33_000.txt", 2) #Initialize COMS for M33 and M31
        M31COM = CenterOfMass("M31_000.txt", 2)

        M33_COMP = M33COM.COM_P(0.1) #Initialize COMP and COMV for both M33 and M31
        M33_COMV = M33COM.COM_V(M33_COMP[0],M33_COMP[1],M33_COMP[2])
        M31_COMP = M31COM.COM_P(0.1)
        M31_COMV = M31COM.COM_V(M31_COMP[0],M31_COMP[1],M31_COMP[2])
        
        self.x = M33_COMP[0].value - M31_COMP[0].value #Initialize position and velocity vectors
        self.y = M33_COMP[1].value - M31_COMP[1].value
        self.z = M33_COMP[2].value - M31_COMP[2].value
        self.vx = M33_COMV[0].value - M31_COMV[0].value
        self.vy = M33_COMV[1].value - M31_COMV[1].value
        self.vz = M33_COMV[2].value - M31_COMV[2].value

        self.rdisk = 5  #Initialize mass and radius values for particle types
        self.Mdisk = 0.12 
        self.rbulge = 1 
        self.Mbulge = 0.019 
        self.rhalo = 25 
        self.Mhalo = 1.921

    def HernquistAccel(self, M, r_a, x, y, z):
        accel = np.zeros(3)
        accel[0] = -(G.value*M / (x*(r_a + x)**2))
        accel[1] = -(G.value*M / (y*(r_a + y)**2))
        accel[2] = -(G.value*M / (z*(r_a + z)**2))
        return accel

    def MiyamotoNagaiAccel(self, M, rd, x, y, z):
        R = np.sqrt(x**2 + y**2)
        B = rd + np.sqrt(z**2 + (rd/5.)**2)
        accel = np.zeros(3)
        accel[0] = -((G.value*M)/(R**2 + B**2)**1.5) * x
        accel[1] = -((G.value*M)/(R**2 + B**2)**1.5) * y
        accel[2] = -((G.value*M*B)/((R**2 + B**2)**1.5 * (B - rd))) * z
        return accel

    def M31Accel(self, x, y, z):
        bulgeAccel = self.HernquistAccel(self.Mbulge, self.rbulge, self.x, self.y, self.z)
        haloAccel = self.HernquistAccel(self.Mhalo, self.rhalo, self.x, self.y, self.z)
        diskAccel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, self.x, self.y, self.z)

        totalAccel = np.zeros(3)
        for i in range(np.size(totalAccel)):
            totalAccel[i] = bulgeAccel[i] + haloAccel[i] + diskAccel[i]
        return totalAccel

    def LeapFrog(self, dt, x, y, z, vx, vy, vz):
        r = [x, y, z] #initialize arrays of distance and velocity for convenience
        v = [vx, vy, vz]
        print("LeapFrog Check")
        for i in range(np.size(r)): #Calculate radius half steps
            r[i] += v[i]*(dt/2.)

        a = self.M31Accel(r[0], r[1], r[2]) #Calculate half step acceleration values
        for i in range(np.size(v)): #Calculate next velocity step
            v[i] += a[i]*dt

        for i in range(np.size(r)): #Calculate full radius step
            r[i] += v[i]*(dt/2.)

        return r, v

    def OrbitIntegrator(self, t0, dt, tf): #t0 is starting time, tf is final time, dt is time interval

        ID = np.arange(t0, tf, dt)
        orbit = np.zeros([ID.size, 7])
        orbit[0] = [t0, self.x, self.y, self.z, self.vx, self.vy, self.vz]
        r, v = self.LeapFrog(dt, self.x, self.y, self.z, self.vx, self.vy, self.vz)   
        
        i = 1
        t = t0 + dt
        while (t < tf):
            orbit[i] = [t, r[0], r[1], r[2], v[0], v[1], v[2]]
            r, v = self.LeapFrog(dt, r[0], r[1], r[2], v[0], v[1], v[2])
            i += 1
            t += dt

        np.savetxt(self.filename, orbit, header="time  x  y  z  vx  vy  vz", comments='#', fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
        return

Orbit = M33AnalyticOrbit("Orbit.txt")

Orbit.OrbitIntegrator(0, 1, 11)