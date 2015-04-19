#!/usr/bin/env python
#
# Script calculates Cs and Cp for RFID antennas
# in Series-Parallel topology. 
#
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

w = 2.0*math.pi*13.56e6

# antenna topology is
# Zs(=Rs AND Cs) AND ( Zp(=Cp) OR Zl(=Rl AND L) ) 
#
# --Rs--Cs-------
# |         |   |
# |         |  Rl
# U0        Cp  |
# |         |   L
# |         |   |
# --------------
Rs = 50.0   # Generator U0 output impedance
Rin = 50.0  # This is the input inpedance of Cs/Cp/Rl/L circuit
Rl = 2    
L = 1.5e-6

#!!! do not change anython below here !!!#

# calculation for Cp (solve quadratic equation)
a = pow(w,2)*(pow(w,2)*pow(L,2)+pow(Rl,2))
b = -2*pow(w,2)*L
c = 1-Rl/Rin

Cp1_1 = (- b - math.sqrt(pow(b,2)-4*a*c))/(2*a)
Cp1_2 = (- b + math.sqrt(pow(b,2)-4*a*c))/(2*a)
#print("Cp1_1:",Cp1_1)
#print("Cp1_2:",Cp1_2)

# calculation for Cs
Cp = Cp1_1
Cs = 1-2*pow(w,2)*L*Cp+pow(w,2)*pow(Cp,2)*(pow(w,2)*pow(L,2)+Rl)
Cs = Cs/(pow(w,2)*(L-Cp*(pow(w,2)*pow(L,2)+pow(Rl,2))))

if Cs < 0:
    Cp=Cp_2
    Cs = 1-2*pow(w,2)*L*Cp+pow(w,2)*pow(Cp,2)*(pow(w,2)*pow(L,2)+Rl)
    Cs = Cs/(pow(w,2)*(L-Cp*(pow(w,2)*pow(L,2)+pow(Rl,2))))

# print results
print("Cp: %.2e F" % Cp)
print("Cs: %.2e F" % Cs)
Zin= Rl/(1-2*pow(w,2)*L*Cp+pow(w,2)*pow(Cp,2)*(pow(w,2)*pow(L,2)+pow(Rl,2)))
print("Verifikation Zin: %.3f Ohm"% Zin)

# calculate antenna current
s = complex(0,w)
Zp = 1/(s*Cp)
Zs = Rs+1/(s*Cs)
Zl = Rl+s*L
U0=1
Zabs = abs(Zs+Zl*Zs/Zp+Zl)
Il = U0/Zabs
print("Spulenstrom Il: %.3fmA" % Il)

