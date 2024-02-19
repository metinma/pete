#PETE 218 Homework 2
import numpy as np #Dependency. pip install numpy
import matplotlib.pyplot as plt #Dependency. pip install matplotlib
import PVTProps as pvt #Dependency. Provided

#Constants for Question 1
SN  = [2, 4, 4, 2, 0, 4, 4] #Student number
GG  = 0.65 + (SN[6] + 1) / 50 #Total gas gravity
CO2 = (SN[4] + 1) / 200 #CO2 mole fraction
H2S = (SN[5] + 1) / 200 #H2S mole fraction
N2  = (SN[3] + 1) / 200 #N2 mole fraction
T   = np.array([60, 200, 400]) #F Temperature
PR  = np.array([15, 10000]) #psi Pressure range
P   = np.linspace(PR[0], PR[1], 100) #psi Pressure vector
TR  = np.array([(T + 460),] * np.size(P)).T #R Temperature

#Constants for Question 2
API  = 20 + SN[6]*3 #API Gravity
RSB  = (API*6)**1.2 + (SN[5]-5)*20 #Rsb
PGG  = 0.65 + (SN[6]+1)/50 #Produced gas gravity
TO   = 120 + SN[3]*3 #Oil temperature
TSEP = 120 #Seperator temperature
PSEP = 125 #Seperator pressure

#Constants for Question 3
Y   = np.array([(SN[5]+1)*7/100, (SN[6]+1)*2/100, 1.0 - (((SN[6]+1)*2/100)+((SN[5]+1)*7/100))]) #Mole Fraction for CH4, C2H6, C3H8
TC  = np.array([-116.67+460, 89.92+460, 205.05 + 460]) #Critical Temperature for CH4, C2H6, C3H8
PC  = np.array([666.4, 705.5, 616.0]) #Critical Pressure for CH4, C2H6, C3H8
AF  = np.array([0.0104, 0.0979, 0.1522]) #Acentric Factor for CH4, C2H6, C3H8
R   = 10.731 #Universal gas constant

#Correction methods
sCorrCriticProp = "Standing"
sCorrNonHC      = "WichertAziz"
sCorrViscosity  = "LeeGonzalesEakin"
sCorrelation    = "Standing"
sCorrelationAlt = "VasquezBeggs"

#Calculated pseudo critical properties
tpc = pvt.Tpc(GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
ppc = pvt.Ppc(GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)

#Calculated pseudo reduced critical properties 
trpc = TR / tpc
prpc = P / ppc

#Z-Factor calculation
zFactor = np.zeros((np.size(T), np.size(P)))
for i in range(np.size(T)):
    for j in range(np.size(P)):
        zFactor[i,j] = pvt.ZFactor(trpc[i,j], prpc[j])

#Gas compressibility calculation
cg = np.zeros((np.size(T), np.size(P)))
for i in range(np.size(T)):
    for j in range(np.size(P)):
        cg[i,j] = pvt.Cg(T[i], P[j], GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)

#Gas viscosity calculation
ug = np.zeros((np.size(T), np.size(P)))
for i in range(np.size(T)):
    for j in range(np.size(P)):
        ug[i,j] = pvt.Ug(T[i], P[j], GG, N2, CO2, H2S, None, None, sCorrViscosity, sCorrCriticProp, sCorrNonHC)

#Solution Gas Oil Ratio calculation
rss = np.zeros((np.size(P)))
rsv = np.zeros((np.size(P)))
for i in range(np.size(P)):
    rss[i] = pvt.Rs(API, RSB, PGG, TO, P[i], sCorrelation, TSEP, PSEP)
    rsv[i] = pvt.Rs(API, RSB, PGG, TO, P[i], sCorrelationAlt, TSEP, PSEP)

#Oil Formation Volume Factor calculation
bos = np.zeros((np.size(P)))
bov = np.zeros((np.size(P)))
for i in range(np.size(P)):
    bos[i] = pvt.Bo(API, RSB, PGG, TO, P[i], sCorrelation, TSEP, PSEP)
    bov[i] = pvt.Bo(API, RSB, PGG, TO, P[i], sCorrelationAlt, TSEP, PSEP)

#Oil Viscosity calculation
uo = np.zeros((np.size(P)))
for i in range(np.size(P)):
    uo[i] = pvt.Uo(API, RSB, PGG, TO, P[i])

#Oil Compressibility calculation
cos = np.zeros((np.size(P)))
cov = np.zeros((np.size(P)))
for i in range(np.size(P)):
    cos[i] = pvt.Co(API, RSB, PGG, TO, P[i], sCorrelation, TSEP, PSEP)
    cov[i] = pvt.Co(API, RSB, PGG, TO, P[i], sCorrelationAlt, TSEP, PSEP)

#Z-Factor calculation using cubic form of SRK EOS
def zfactor_srk_eos():
    t = 200 + 460   #Temperature, R
    p = np.array([15, 1000, 2000, 3000, 4000, 5000, 6000]) #Pressure, psia
    tpc     = Y[0]*TC[0] + Y[1]*TC[1] + Y[2]*TC[2]
    ppc     = Y[0]*PC[0] + Y[1]*PC[1] + Y[2]*PC[2]
    tr = t / tpc    #Reduced temperature
    app     = Y[0]*AF[0] + Y[1]*AF[1] + Y[2]*AF[2]
    a       = 0.42727 * R**2 * tpc**2 / ppc
    b       = 0.08664 * R * tpc / ppc
    alpha   = (1 + (0.480 + 1.574 * app - 0.176 * app**2) * (1 - tr**0.5))**2
    aa      = (alpha * a * p) / (R**2 * t**2)
    bb      = (b * p) / (R * t)
    
    c = np.zeros((np.size(p), 4))
    for i in range(np.size(p)):
        c[i,0] = 1
        c[i,1] = -1
        c[i,2] = aa[i] - bb[i] - bb[i]**2
        c[i,3] = - aa[i] * bb[i]
    
    roots = np.zeros(np.size(p))
    for i in range(np.size(p)):
        roots[i] = np.sort(np.real(np.roots(c[i,:])))[-1]   #HOW ???
    
    return p, roots

#Figure generation, and cleanup of unnecessary subplots
FIGSIZE = 5 #Size of the figure    
fig, ax = plt.subplots(3, 4, figsize=(4*FIGSIZE, 3*FIGSIZE))
fig.delaxes(ax[0][3])
fig.delaxes(ax[2][1])
fig.delaxes(ax[2][2])
fig.delaxes(ax[2][3])

#Plots for Question 1
for i in range(np.size(T)):
    ax[0,0].semilogx(P, zFactor[i,:], label=f"{T[i]} °F")
    ax[0,0].set_xlabel("Pressure, psi")
    ax[0,0].set_ylabel("Z-Factor")    
    ax[0,0].legend()
    
for i in range(np.size(T)):
    ax[0,1].semilogx(P, cg[i,:], label=f"{T[i]} °F")
    ax[0,1].set_xlabel("Pressure, psi")
    ax[0,1].set_ylabel("Gas Compressibility")
    ax[0,1].legend()

for i in range(np.size(T)):
    ax[0,2].semilogx(P, ug[i,:], label=f"{T[i]} °F")
    ax[0,2].set_xlabel("Pressure, psi")
    ax[0,2].set_ylabel("Gas Viscosity")
    ax[0,2].legend()

#Plots for Question 2
ax[1,0].semilogx(P, rss, label="Standing")
ax[1,0].semilogx(P, rsv, label="Vasquez-Beggs")
ax[1,0].set_xlabel("Pressure, psi")
ax[1,0].set_ylabel("Solution Gas Oil Ratio")
ax[1,0].legend()

ax[1,1].semilogx(P, bos, label="Standing")
ax[1,1].semilogx(P, bov, label="Vasquez-Beggs")
ax[1,1].set_xlabel("Pressure, psi")
ax[1,1].set_ylabel("Oil Formation Volume Factor")
ax[1,1].legend()

ax[1,2].semilogx(P, uo)
ax[1,2].set_xlabel("Pressure, psi")
ax[1,2].set_ylabel("Oil Viscosity")

ax[1,3].semilogx(P, cos, label="Standing")
ax[1,3].semilogx(P, cov, label="Vasquez-Beggs")
ax[1,3].set_xlabel("Pressure, psi")
ax[1,3].set_ylabel("Oil Compressibility")
ax[1,3].legend()

#Plots for Question 3
ax[2,0].semilogx(zfactor_srk_eos()[0], zfactor_srk_eos()[1], marker="o")
ax[2,0].set_xlabel("Pressure, psi")
ax[2,0].set_ylabel("Z-Factor")

plt.show() #Render the figure