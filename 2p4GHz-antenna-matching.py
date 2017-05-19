#!/usr/bin/env python
#
# Script calculates optimal S11 (reflection) values
# for a 2.4 GHz antenna matching

import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
import pylab as pl

Zw = complex(50,0)

def Zres(R,w): return complex(R,0)
def Zind(L,w): return complex(0,w * L)
def Zcap(C,w): return 1.0 / (complex(0,w) * C)
def par(Z1,Z2): return Z1*Z2/(Z1+Z2)
def ser(Z1,Z2): return Z1+Z2
def S11_abs(Z): return abs((Z-Zw)/(Z+Zw))

# antenna topology is
print("2p4GHz-antenna-matching.py")
print("--------------------------")
print("")
print("Matching circuit is PI-filter:")
print(" ----Ls1--------------")
print("  |          |       |")
print("  |          |       |")
print(" Cp1        Cp2    Zant")
print("  |          |       |")
print("  |          |       |")
print(" ---------------------")

Cp1_all = pl.frange(0.1e-12,4.0e-12,0.1e-12)
Ls1_all = pl.frange(0.1e-9,6.0e-9,0.2e-9)
Cp2_all = pl.frange(0.1e-12,4.0e-12,0.1e-12)

w_low =  2*math.pi*2.4e9
w_mid =  2*math.pi*2.44e9
w_high = 2*math.pi*2.48e9
Z0 = {w_low:complex(79,67),w_mid:complex(116,-70),w_high:complex(31,-55)}
S11 = {w_low:complex(0,0),w_mid:complex(0,0),w_high:complex(0,0)}
Z = {}
# goes to 31-j53 Cp1=1.0pF Ls1=4.8nH Cp2=0.7pF
S11_min = 1000

for Cp2 in Cp2_all:
    # print('-Cp2={:.2e}'.format(Cp2))
    for Ls1 in Ls1_all:
        # print('--Ls1={:.2e}'.format(Ls1))
        for Cp1 in Cp1_all:
            # print('---Cp1={:.2e}'.format(Cp1))
            S11_new = 0
            for w in Z0:
                Z[w] = par(Z0[w],Zcap(Cp2,w))
                Z[w] = ser(Z[w],Zind(Ls1,w))
                Z[w] = par(Z[w],Zcap(Cp1,w))
                S11[w] = S11_abs(Z[w])
                S11_new += S11[w]**2
                # print('    w={:.3e} Z[w].real={:.3f} S11[w]={:.3f}'.format(w,Z[w].real,S11[w]))
                # print('    S11_new=={:.3f}'.format(S11_new))

            if S11_new < S11_min:
                S11_min = S11_new
                print('New min: Cp1={:.2e} Ls1={:.2e} Cp2={:.2e} S11_min={:.5f}'.format(Cp1,Ls1,Cp2,S11_min))
        
    


