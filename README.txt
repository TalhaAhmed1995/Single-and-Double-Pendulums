{\rtf1\ansi\ansicpg1252\cocoartf1504\cocoasubrtf830
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 The simple pendulum FDM solutions can be produced by running the singlePendulum.py module. This contains all four methods plus the large angle form of Explicit Euler. \
All of these methods require 3 parameters, which are a step-size, h, damping coefficient, D, and time-span, tMax. Once inputted, the function will produce plots of the angle and energy against time. A stability analysis report will also be printed.\
An example of the parameters used to replicate RK4 for 0.2 damping as shown in the report is:\
\
\'93Input: RK4(h=0.02, D=0.2, tMax=80)\'94\
\'93Output: *Plots of angle and energy against time*\
		This solution is stable\
		Maximum Energy: 0.005\
		Highest Energy in first 20 Energy Oscillations: 0.005\
		Percentage Difference: 0.0\'94\
\
Similarly, the doublePendulum.py module outputs results for the double pendulum, using just the RK4 method. The function requires step-size, h, mass ratio, R, damping, G, and time-span, tMax. A particular combination used in the report is RKDouble(h=0.02, R=1.0, G=0.0, tMax = 80).\
\
The analysis.py module doesn\'92t require any manual inputting. It simply contains functions for plotting, energy calculating and stability testing.\
\
}