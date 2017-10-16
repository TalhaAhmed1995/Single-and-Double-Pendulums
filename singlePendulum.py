# -*- coding: utf-8 -*-
import numpy as np
import analysis as an

"""This module is strictly used for the simple pendulum. The four FDMs used are:

1) Explicit Euler
2) Leapfrog
3) RK4
4) Implicit Euler

Plotting, energy calculations and stability checks are imported from the analysis.py module"""
   
def euler(h, D, tMax):
    """Calculates the solution using the Explicit Euler method, and plots against time"""
    
    t = np.arange(0, tMax, h) 
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1 #Using small initial angle
    
    for i in range(1, t.size): #Iteratively obtaining theta values using Explicit Euler
        v[i] = v[i-1] + (-theta[i-1] - D*v[i-1]) * h
        theta[i] = theta[i-1] + (v[i-1]) * h
        
    trueValues = an.realSolution(theta[0], D, t) #Obtaining the analytic solution
    
    fig = an.singlePlotting('Explicit Euler', theta, t, h, D, trueValues) #Plotting the numerical solution
    energies = an.singleEnergy(theta, v, t, fig) #Calculating energy values for numerical solution obtained
    an.stabilityCheck(energies) #Performing stability check on this FDM
    
def leapfrog(h, D, tMax):
    """Calculates the solution using the Leapfrog method, and plots against time"""
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1 #Using small initial angle
    
    v[1] = v[0] + (-theta[0] - D*v[0]) * h #First value in leapfrog method must actually be calculated using Euler Forward Scheme
    theta[1] = theta[0] + (v[0]) * h
    
    for i in range(2, t.size): #Iteratively obtaining theta values using Leapfrog
        v[i] = v[i-2] + 2 * (-theta[i-1] - D*v[i-1]) * h
        theta[i] = theta[i-2] + 2 * (v[i-1]) * h
        
    trueValues = an.realSolution(theta[0], D, t) #Obtaining the analytic solution
    
    fig = an.singlePlotting('Leapfrog', theta, t, h, D, trueValues) #Plotting the numerical solution
    energies = an.singleEnergy(theta, v, t, fig) #Calculating energy values for numerical solution obtained
    an.stabilityCheck(energies) #Performing stability check on this FDM
    
   
def RK4(h, D, tMax):
    """Calculates the solution using the RK4 method, and plots against time"""
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1 #Using small initial angle
    
    for i in range(1, t.size): #Iteratively obtaining theta values using RK4
        kv1 = -theta[i-1] - D * v[i-1]
        ktheta1 = v[i-1]
        v1 = v[i-1] + kv1 * h/2.
        theta1 = theta[i-1] + ktheta1 * h/2.
        kv2 = -theta1 - D * v1
        ktheta2 = v1
        v2 = v[i-1] + kv2 * h/2.
        theta2 = theta[i-1] + ktheta2 * h/2.
        kv3 = -theta2 - D * v2
        ktheta3 = v2
        v3 = v[i-1] + kv3 * h
        theta3 = theta[i-1] + ktheta3 * h
        kv4 = -theta3 - D * v3
        ktheta4 = v3
        
        v[i] = v[i-1] + h/6. * (kv1 + 2*kv2 + 2*kv3 + kv4)
        theta[i] = theta[i-1] + h/6. * (ktheta1 + 2*ktheta2 + 2*ktheta3 + ktheta4) #Final theta value calculated here per step
        
    trueValues = an.realSolution(theta[0], D, t) #Obtaining the analytic solution
    
    fig = an.singlePlotting('RK4', theta, t, h, D, trueValues) #Plotting the numerical solution
    energies = an.singleEnergy(theta, v, t, fig) #Calculating energy values for numerical solution obtained
    an.stabilityCheck(energies) #Performing stability check on this FDM

def implicit(h, D, tMax):
    """Calculates the solution using the Implicit Euler method, and plots against time"""
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.1 #Using small initial angle
    
    denominator = 1 / (h**2 + D*h +1) #Coefficient of inverse matrix used in implicit euler method
    for i in range(1, t.size): #Iteratively obtaining theta values using Implicit Euler
        v[i] = denominator * (-h*theta[i-1] + v[i-1])
        theta[i] = theta[i-1] + v[i]*h
        
    trueValues = an.realSolution(theta[0], D, t) #Obtaining the analytic solution
    
    fig = an.singlePlotting('Implicit Euler', theta, t, h, D, trueValues) #Plotting the numerical solution
    energies = an.singleEnergy(theta, v, t, fig) #Calculating energy values for numerical solution obtained
    an.stabilityCheck(energies) #Performing stability check on this FDM

def largeEuler(h, D, tMax):
    """Calculates the solution using the Explicit Euler method (FOR LARGE ANGLES), and plots against time"""
    
    t = np.arange(0, tMax, h)
    v = np.zeros(t.size)
    theta = np.zeros(t.size)
    theta[0] = 0.75 * np.pi #Using large initial angle
    
    for i in range(1, t.size): #Iteratively obtaining theta values using Explicit Euler
        v[i] = v[i-1] + (-np.sin(theta[i-1]) - D*v[i-1]) * h #This time the small angle. approx is not made, as seen by the sine function
        theta[i] = theta[i-1] + (v[i-1]) * h
    
    trueValues = an.realSolution(theta[0], D, t) #Obtaining the analytic solution
    
    fig = an.singlePlotting('Explicit Euler (Large Angle)', theta, t, h, D, trueValues) #Plotting the numerical solution
    energies = an.singleEnergy(theta, v, t, fig) #Calculating energy values for numerical solution obtained
    an.stabilityCheck(energies) #Performing stability check on this FDM