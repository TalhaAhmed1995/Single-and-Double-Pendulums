import numpy as np
from matplotlib import pyplot as plt

"""A module where all the functions required for analysis have been created.
   This includes analytics solutions, plotting, calculating energies and performing stability checks for both pendulum types"""
   
def realSolution(theta0, D, t):
    """Calculates the analytic solution for the simple pendulum problem"""
    
    alpha = np.sqrt(1 - (D**2)/4)
    trueValue = np.exp(-D*t/2) * (theta0*np.cos(alpha*t) + ((D*theta0)/(2*alpha))*np.sin(alpha*t)) 
    
    return trueValue #Outputting the analytically calculated theta values
   

def singlePlotting(method, theta, t, h, D, trueValues):
    """Plots numerical solutions of the SIMPLE pendulum"""
    
    #Plotting the FDM solution for single pendulum
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, label='Theta (Numerial Solution)')
    ax1.plot(t, trueValues, '--', label='Theta (Analytic Solution)') #Analytic solution for comparison
    plt.title(method + ' for Simple Pendulum with h=' +str(h) + ' and D=' +str(D))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    
    return fig #Returned to add a further subplot for the energy
    
def doublePlotting(method, theta, phi, t, R, G):
    """Plots numerical solutions of the DOUBLE pendulum"""
    
    #Plotting the FDM solution for double pendulum
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(t, theta, 'b', label='Theta')
    ax1.plot(t, phi, 'r', label='Phi')
    plt.title(method + ' for Double Pendulum with R = ' +str(R) + ' and G = ' +str(G))
    plt.xlabel('Time')
    plt.ylabel('Angle')
    plt.legend()
    plt.grid()
    plt.show()    
    
    return fig #Returned to add a further subplot for the energy

def singleEnergy(theta, v, t, fig):
    """Calculates energy values for numerical solutions of the SIMPLE pendulum"""
    
    #Creating arrays that are the same size as the array of time values
    KEs = np.zeros(t.size)
    PEs = np.zeros(t.size)
    energies = np.zeros(t.size)
    
    if theta[0] < 1: #Calculating the energy under the small angle approx. if the initial angle is less than unity. 
        for i in range(t.size):
            KE = (v[i]**2)/2.
            PE = (theta[i]**2)/2. #Small angle approx. made for the PE
            
            energies[i] = KE + PE #Total energy calculated using sum of KE and PE
            KEs[i] = KE
            PEs[i] = PE
    else: 
        for i in range(t.size): #This version of energy is used for large angles
            KE = (v[i]**2)/2.
            PE = (1 - np.cos(theta[i])) #Small angle approx. not made here
            
            energies[i] = KE + PE
            KEs[i] = KE
            PEs[i] = PE
            
    #Plotting the energies calculated over the specified range
    ax2 = fig.add_subplot(212)
    ax2.plot(t, energies, 'g', label='Total Energy')
    ax2.plot(t, KEs, 'b', label='Kinetic Energy')
    ax2.plot(t, PEs, 'r', label='Potential Energy')
    plt.ylim(0, max(energies)+0.0005)
    plt.title('Energy over time')
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.legend()
    plt.grid()
    
    return energies #Returning energy array to be using for stability checks
    
def doubleEnergy(theta, phi, w, v, R, G, t, fig):
    """Calculates energy values for numerical solutions of the DOUBLE pendulum"""
    
    #Creating arrays that are the same size as the array of time values
    energies = np.zeros(t.size)
    KEs = np.zeros(t.size)
    PEs = np.zeros(t.size)
    
    for i in range(t.size):
        
        KE = 1/2. * w[i]**2 + 1/2. * R * (w[i]**2 + v[i]**2 + 2*w[i]*v[i]) #Small angle approx. made for both KE and PE
        PE = (theta[i]**2)/2. + R*(theta[i]**2)/2. + R*(phi[i]**2)/2.
        
        energies[i] = KE + PE #Total energy calculated using sum of KE and PE
        KEs[i] = KE
        PEs[i] = PE
        
    #Plotting the energies calculated over specified range
    ax2 = fig.add_subplot(212)
    ax2.plot(t, energies, 'g', label='Total Energy')
    ax2.plot(t, KEs, 'b', label='Kinetic Energy')
    ax2.plot(t, PEs, 'r', label='Potential Energy')
    plt.ylim(0, max(energies)+0.0005)
    plt.title('Energy of Double Pendulum with R = ' +str(R) + ' and G = ' +str(G))
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.legend()
    plt.grid()
    plt.show()
    
    return energies #Returning energy array to be using for stability checks


def stabilityCheck(energies):
    """Performs a stability test on predetermined energy values of a particular numerical solution"""

    first20Peaks = [energies[0]] #A list of the first 20 peaks found if the energy oscillates. The initial energy is also included
    
    i = 0
    while len(first20Peaks) < 19: #Finding the first 20 peaks in any potential energy oscillations
        i += 1
        if i > energies.size - 2: #If no peaks are found in the energy, just append the first 20 energy values of the solution
            for i in range(1, 20):
                first20Peaks.append(energies[i])
        elif energies[i] == energies[i-1]: #If the energy value is the same as the previous value, it's the same peak
            pass
        elif energies[i] > energies[i-1] and energies[i] >= energies[i+1]: #If neighbouring energy values are lower, the current value is a peak
            first20Peaks.append(energies[i]) #Appending energy value of peak to list
        
    highestInitialEnergy = max(first20Peaks) #Finding the max energy from the first 20 energy oscillations
    
    maximumEnergy = max(energies) #The maximum energy across the whole time-span
    
    percentage = ((float(maximumEnergy) - float(highestInitialEnergy))/float(highestInitialEnergy)) *100. #Calculating the percentage difference
    
    #A tolerance of 1.0%. If higher than this, solution is unstable, otherwise it is stable
    stability = "unstable" if percentage > 1.0 else "stable"
    
    #Outputting key variables of stability analysis
    print 'This solution is',  stability
    print 'Maximum Energy:', maximumEnergy
    print 'Highest Energy in first 20 Energy Oscillations:', highestInitialEnergy
    print 'Percentage Difference:', percentage  