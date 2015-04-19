#!/usr/bin/env python
#
# Script calculates Cs and Cp for RFID antennas
# in Series-Parallel topology.
#
# In reference to equations (16) and (17) from paper:
# "Dependence of RFID Reader Antenna Design on Read Out Distance"
# IEEE Transactions on Ant. and Prop., Vol.56, Dec. 2008

import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

w = 2.0*math.pi*13.56e6

# antenna topology is
# Zs(=Rs AND Cs) AND ( Zp(=Cp) OR Zl(=Rl AND L) ) 
#
print("rfid-antenna-calculator.py")
print("--------------------------")
print("")
print("circuit:")
print(" --Rs--Cs-------")
print(" |         |   |")
print(" |         |  Rl")
print(" U0(=1V)   Cp  |")
print(" |         |   L")
print(" |         |   |")
print(" ---------------")

Rs = 50.0   # Generator U0 output impedance
Rin = 50.0  # This is the input inpedance of Cs/Cp/Rl/L circuit
Rl = 2    
L = 1.5e-6

#!!! do not change anython below here !!!#
print("")
print("using Rs=%.02f [Ohm], Rl=%.02f [Ohm], L=%.02e [Henry]" % (Rs,Rl,L))
print("")

# calculation for Cp (solve quadratic equation)
a = pow(w,2)*(pow(w,2)*pow(L,2)+pow(Rl,2))
b = -2*pow(w,2)*L
c = 1-Rl/Rin

Cp1 = (- b - math.sqrt(pow(b,2)-4*a*c))/(2*a)
Cp2 = (- b + math.sqrt(pow(b,2)-4*a*c))/(2*a)

# calculation for Cs
Cp = Cp1
Cs1 = 1-2*pow(w,2)*L*Cp+pow(w,2)*pow(Cp,2)*(pow(w,2)*pow(L,2)+Rl)
Cs1 = Cs1/(pow(w,2)*(L-Cp*(pow(w,2)*pow(L,2)+pow(Rl,2))))

Cp=Cp2
Cs2 = 1-2*pow(w,2)*L*Cp+pow(w,2)*pow(Cp,2)*(pow(w,2)*pow(L,2)+Rl)
Cs2 = Cs2/(pow(w,2)*(L-Cp*(pow(w,2)*pow(L,2)+pow(Rl,2))))

# print results
print("Solution 1: Cp=%.2e Cs=%.2e [F]" % (Cp1,Cs1))
print("Solution 2: Cp=%.2e Cs=%.2e [F]" % (Cp2,Cs2))

if Cp1<0 or Cs1<0:
    print("-> using solution 2")
    Cs = Cs2
    Cp = Cp2
else:
    print("-> using solution 1")
    Cs = Cs1
    Cp = Cp1

Zin= Rl/(1-2*pow(w,2)*L*Cp+pow(w,2)*pow(Cp,2)*(pow(w,2)*pow(L,2)+pow(Rl,2)))
print("Verification: Zin=%.3f [Ohm]"% Zin)

# calculate antenna current
s = complex(0,w)
Zp = 1/(s*Cp)
Zs = Rs+1/(s*Cs)
Zl = Rl+s*L
U0=1
Zabs = abs(Zs+Zl*Zs/Zp+Zl)
Il = U0/Zabs
print("Max. RFID coil (L) current: I=%.3f [A]" % Il)

