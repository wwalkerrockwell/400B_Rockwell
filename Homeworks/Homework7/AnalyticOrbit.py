import numpy as np
import astropy.units as u
from astropy.constants import G
import matplotlib.pyplot as plt
import matplotlib
from Readfile import Read
from CenterOfMass import CenterOfMass


class M33AnalyticOrbit:

    def __init__(self, filename):
        self.G = 4.498768e-6
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
        self.Mdisk = 0.12e12
        self.rbulge = 1 
        self.Mbulge = 0.019e12
        self.rhalo = 25 
        self.Mhalo = 1.921e12

    def HernquistAccel(self, M, r_a, x, y, z):
        accel = np.zeros(3)
        r = np.sqrt(x**2 + y**2 + z**2) #Magnitude of position vector
        accel[0] = -((self.G*M) / (r*(r_a + r)**2)) * x
        accel[1] = -((self.G*M) / (r*(r_a + r)**2)) * y
        accel[2] = -((self.G*M) / (r*(r_a + r)**2)) * z


        return accel

    def MiyamotoNagaiAccel(self, M, rd, x, y, z):
        R = np.sqrt(x**2 + y**2)
        B = rd + np.sqrt(z**2 + (rd/5.)**2)
        accel = np.zeros(3)
        accel[0] = -((self.G*M)/(R**2 + B**2)**1.5) * x
        accel[1] = -((self.G*M)/(R**2 + B**2)**1.5) * y
        accel[2] = -((self.G*M*B)/((R**2 + B**2)**1.5 * (B - rd))) * z
        return accel

    def M31Accel(self, x, y, z):
        bulgeAccel = self.HernquistAccel(self.Mbulge, self.rbulge, x, y, z)
        haloAccel = self.HernquistAccel(self.Mhalo, self.rhalo, x, y, z)
        diskAccel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, x, y, z)


        totalAccel = np.zeros(3)
        for i in range(np.size(totalAccel)):
            totalAccel[i] = bulgeAccel[i] + haloAccel[i] + diskAccel[i]
        return totalAccel

    def LeapFrog(self, dt, x, y, z, vx, vy, vz):
        r = [x, y, z] #initialize arrays of distance and velocity for convenience
        v = [vx, vy, vz]
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
        print(ID.size)
        orbit[0] = [t0, self.x, self.y, self.z, self.vx, self.vy, self.vz]
        r, v = self.LeapFrog(dt, self.x, self.y, self.z, self.vx, self.vy, self.vz)   
        
        i = 1
        t = t0 + dt
        while (i < ID.size):
            print(i)
            orbit[i] = [t, r[0], r[1], r[2], v[0], v[1], v[2]]
            r, v = self.LeapFrog(dt, r[0], r[1], r[2], v[0], v[1], v[2])
            i += 1
            t += dt


        np.savetxt(self.filename, orbit, header="time  x  y  z  vx  vy  vz", comments='#', fmt=['%.2f','%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
        return

Orbit = M33AnalyticOrbit("Orbit.txt")

#Orbit.OrbitIntegrator(0, 0.1, 10)

Diff = np.genfromtxt("Orbit.txt", dtype = None, names = True, skip_header = 0)


plt.plot(Diff['time'], Diff['vx'], label="X Position")
plt.plot(Diff['time'], Diff['vy'], label="Y Position")
plt.plot(Diff['time'], Diff['vz'], label="Z Position")


plt.xlabel("Time (Gyr)")
plt.ylabel("Velocity (km/s)")

plt.show()

#2. The plots show the x, y, and z positions converging later in its lifespan. 
#We also see spikes in velocity at two points, likely representing collisions with the other galaxies

#3. Gravity, all these graphs account for is relative positions and velocities. It's merely an estimate using some known values.

#4. I'd have to account for M33's relative position and velocities with MW, and combine the effects with respect to M31 as well.