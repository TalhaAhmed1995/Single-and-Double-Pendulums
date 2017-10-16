#Testing Numerical Methods on Pendulum Models

The simple pendulum FDM solutions can be produced by running the singlePendulum.py module. This contains all four methods plus the large angle form of Explicit Euler.

All of these methods require 3 parameters, which are a step-size, h, damping coefficient, D, and time-span, tMax. Once inputted, the function will produce plots of the angle and energy against time. A stability analysis report will also be printed.

## Example Usage

An example of the parameters used to replicate RK4 for 0.2 damping as shown in the report is:

```
Input: RK4(h=0.02, D=0.2, tMax=80)
Output: This solution is stable
		Maximum Energy: 0.005
		Highest Energy in first 20 Energy Oscillations: 0.005
		Percentage Difference: 0.0
```

![Plots of angle and energy against time](/images/RK4021.png?raw=true)

Similarly, the doublePendulum.py module outputs results for the double pendulum, using just the RK4 method. The function requires step-size, h, mass ratio, R, damping, G, and time-span, tMax. A particular combination used in the report is `RKDouble(h=0.02, R=1.0, G=0.0, tMax = 80)`.

The analysis.py module doesn't require any manual inputting. It simply contains functions for plotting, energy calculating and stability testing.
