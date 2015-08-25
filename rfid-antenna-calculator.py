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

def Zres(R,w): return complex(R,0)
def Zind(L,w): return complex(0,w * L)
def Zcap(C,w): return 1.0 / (complex(0,w) * C)
def par(Z1,Z2): return Z1*Z2/(Z1+Z2)

# antenna topology is
# Zs(=Rs AND Cs) AND ( Zp(=Cp) OR Zl(=Rl AND L) ) 
#
print("rfid-antenna-calculator.py")
print("--------------------------")
print("")
print("circuit:")
print(" --Rq--Cs----------")
print(" |         |  |   |")
print(" |         |  |  Rl")
print(" U0(=1V)   Cp Rp  |")
print(" |         |  |   L")
print(" |         |  |   |")
print(" ------------------")

U0=1
Rq=50.0   # generator U0 output impedance
Rin=50.0  # this is the desired input inpedance of Cs/Cp/Rl/L circuit
Rl=1.5    # resistance of the coil
Rp=5e3    # parallel resistor
L=1.3e-6

w0=2.0*math.pi*13.56e6
s = complex(0,w0)
Rlp=(w0*L)**2/Rl
Q=w0*L/Rl

#!!! do not change anython below here !!!#
print("")
print("using Rq=%.02f [Ohm], Rl=%.02f  [Ohm], L=%.02e [Henry]" % (Rq,Rl,L))
print("with Rlp=%.02e [Ohm], Q=%.02f" % (Rlp,Q))
print("")
    
# calculation for Cp (solve quadratic equation)
a=(w0**2)*(w0**2)*(L**2)+(Rl**2)
b=-2*(w0**2)*L
c=1-Rl/Rin

Cp1=(- b - math.sqrt(pow(b,2)-4*a*c))/(2*a)
Cp2=(- b + math.sqrt(pow(b,2)-4*a*c))/(2*a)

# calculation for Cs
Cp=Cp1
Cs1=1-2*pow(w0,2)*L*Cp+pow(w0,2)*pow(Cp,2)*(pow(w0,2)*pow(L,2)+Rl)
Cs1=Cs1/(pow(w0,2)*(L-Cp*(pow(w0,2)*pow(L,2)+pow(Rl,2))))

Cp=Cp2
Cs2=1-2*pow(w0,2)*L*Cp+pow(w0,2)*pow(Cp,2)*(pow(w0,2)*pow(L,2)+Rl)
Cs2=Cs2/(pow(w0,2)*(L-Cp*(pow(w0,2)*pow(L,2)+pow(Rl,2))))

# print results
print("Solution 1: Cp=%.2e Cs=%.2e [F]" % (Cp1,Cs1))
print("Solution 2: Cp=%.2e Cs=%.2e [F]" % (Cp2,Cs2))

if Cp1<0 or Cs1<0:
    print("-> using solution 2")
    Cs=Cs2
    Cp=Cp2
else:
    print("-> using solution 1")
    Cs=Cs1
    Cp=Cp1

# define the antenna impedances
Zp=1/(s*Cp)
Zs=Rq+1/(s*Cs)
Zl=Rl+s*L

# calculate antenna input impedance
Zin=Rl/(1-2*pow(w0,2)*L*Cp+pow(w0,2)*pow(Cp,2)*(pow(w0,2)*pow(L,2)+pow(Rl,2)))
print("Input impedance (without Rs): Zin=%.2f [Ohm]"% Zin)

# calculate antenna current
Zabs=abs(Zs+Zl*Zs/Zp+Zl)
Il=U0/Zabs
print("Max. RFID coil (L) current: I=%.3f [A]" % Il)

# evaluation of frequency range
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step

A_f=[]
A_Zin_real=[]
A_Zin_imag=[]
A_Il=[]

def evaluate(w):
    Zl = Rl+Zind(L,w)
    Zp = Zcap(Cp,w)
    Zs = Rq+Zcap(Cs,w)

    Zin = Zs+Zp*Zl/(Zp+Zl)    
    Il = U0 / (Zl+Zs+Zs*Zl/Zp)

    A_Zin_real.append(Zin.real)
    A_Zin_imag.append(Zin.imag)
    A_Il.append(abs(Il)*1000)

    # print("f=%.03f:Il=%.03f,Zin=%.03f" % (w/(2*math.pi),abs(Il),Zin.real))


for w in frange(2*math.pi*13e6,2*math.pi*14e6,2*math.pi*0.01e6):
    evaluate(w)
    A_f.append(w/(2*math.pi))


plt.plot(A_f,A_Zin_real,A_f,A_Zin_imag,A_f,A_Il)
plt.show()
