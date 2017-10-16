import numpy as np
import analysis as an

"""This module is used for the double pendulum. Only RK4 is implemented here.
Plotting, energy calculations and stability checks are imported from the analysis.py module"""


def RKDouble(h, R, G, tMax):
    """Calculates the double pendulum solution using the RK4 method, and plots against time"""
    
    t = np.arange(0, tMax, h)
    theta = np.zeros(t.size)
    phi = np.zeros(t.size)
    w = np.zeros(t.size)
    v = np.zeros(t.size)
    theta[0] = 0.1 #Using small initial angle
    
    for i in range(1, t.size): #This method requires four different ODEs to be solved simultaneously using RK4
        
        kv1 = (R+1) * theta[i-1] - (R+1) * phi[i-1] + G * (1-(R**-1)) * w[i-1] - (G*R**-1) * v[i-1]
        kw1 = -(R+1) * theta[i-1] + R * phi[i-1] - G * w[i-1]
        kphi1 = v[i-1]
        ktheta1 = w[i-1]
        
        v1 = v[i-1] + kv1 * h/2.
        w1 = w[i-1] + kw1 * h/2.
        phi1 = phi[i-1] + kphi1 * h/2.
        theta1 = theta[i-1] + ktheta1 * h/2.
        
        kv2 = (R+1) * theta1 - (R+1) * phi1 + G * (1-(R**-1)) * w1 - (G*R**-1) * v1
        kw2 = -(R+1) * theta1 + R * phi1 - G * w1
        kphi2 = v1
        ktheta2 = w1
        
        v2 = v[i-1] + kv2 * h/2.
        w2 = w[i-1] + kw2 * h/2.
        phi2 = phi[i-1] + kphi2 * h/2.
        theta2 = theta[i-1] + ktheta2 * h/2.
        
        kv3 = (R+1) * theta2 - (R+1) * phi2 + G * (1-(R**-1)) * w2 - (G*R**-1) * v2
        kw3 = -(R+1) * theta2 + R * phi2 - G * w2
        kphi3 = v2
        ktheta3 = w2
        
        v3 = v[i-1] + kv3 * h
        w3 = w[i-1] + kw3 * h
        phi3 = phi[i-1] + kphi3 * h
        theta3 = theta[i-1] + ktheta3 * h
        
        kv4 = (R+1) * theta3 - (R+1) * phi3 + G * (1-(R**-1)) * w3 - (G*R**-1) * v3
        kw4 = -(R+1) * theta3 + R * phi3 - G * w3
        kphi4 = v3
        ktheta4 = w3
        
        v[i] = v[i-1] + h/6. * (kv1 + 2*kv2 + 2*kv3 + kv4)
        w[i] = w[i-1] + h/6. * (kw1 + 2*kw2 + 2*kw3 + kw4)
        phi[i] = phi[i-1] + h/6. * (kphi1 + 2*kphi2 + 2*kphi3 + kphi4) #Final phi value calculated here per step
        theta[i] = theta[i-1] + h/6. * (ktheta1 + 2*ktheta2 + 2*ktheta3 + ktheta4) #Final theta value calculated here per step
        
        #(Could have significantly shortened the approach by using matrices)
    
    fig = an.doublePlotting('RK4', theta, phi, t, R, G) #Plotting the numerical solution
    energies = an.doubleEnergy(theta, phi, w, v, R, G, t, fig) #Calculating energy values for numerical solution obtained
    an.stabilityCheck(energies) #Performing stability check on this FDM