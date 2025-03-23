import numpy as np


def main():
    # Surface equipment information
    # Stand Pipe L, Hose L, Swivel L, Kelly L
    SEL = [40, 55, 5, 40]
    # Stand Pipe ID, Hose ID, Swivel ID, Kelly ID
    SEID = [3.5, 2.5, 2.25, 3.25]
    # Mud gallon per minute
    SEGPM = 480
    # Pressure loss
    SEPD = 66

    return "Lubricated!"


# Laziness Machine
class FluidFlowType:
    def __init__(self, n, re, rel, ret, na):
        self.n = n
        self.re = re
        self.rel = rel
        self.ret = ret
        self.na = na

    def laminar_flow(self):
        return self.n / self.re

    def transient_flow(self):
        return ((self.re - self.rel) / 800) * \
            ((((np.log10(self.na) + 3.94) / 50) / self.ret ** ((1.75 - np.log10(self.na)) / 7)) -
             (self.n / self.rel)) + (self.n / self.rel)

    def turbulent_flow(self):
        return ((np.log10(self.na) + 3.93) / 50) / \
            self.re ** ((1.75 - np.log10(self.na)) / 7)


def apd(Q, THETA_300, THETA_3, RHO, TVD, OD, ID, L):
    for i in range(L):
        # Va = annular fluid velocity for the interval (ft/sec)
        va = 0.408 * Q / (OD[i] ** 2 - ID[i] ** 2)
        # na = annular flow behavior index (dimensionless)
        na = 0.5 * np.log10(THETA_300 / THETA_3)
        # Ka = annular consistency factor (poise)
        ka = 5.11 * THETA_300 / 511 ** na
        # μea = annular effective viscosity (cp)
        muea = 100 * ka * (144 * va / (OD[i] - ID[i])) ** (na - 1)
        # Rea = Annular Reynolds number (dimensionless)
        re = 928 * va * (OD[i] - ID[i]) * RHO / \
            (muea * ((2 * na + 1) / (3 * na)) ** na)
        # ReL = the laminar to transitional flow Reynolds number (dimensionless)
        rel = 3470 - 1370 * na
        # ReT = the transitional to turbulent flow Reynolds number (dimensionless)
        ret = 4270 - 1370 * na

        # fa = the annular fanning friction factor (dimensionless)
        if re < rel:
            fa = FluidFlowType(24, re, rel, ret, na).laminar_flow()
        elif ret < re:
            fa = FluidFlowType(24, re, rel, ret, na).turbulent_flow()
        elif rel < re < ret:
            fa = FluidFlowType(24, re, rel, ret, na).transient_flow()

        # Pa = the interval pressure drop (psi)
        apd += fa * va ** 2 * RHO * L[i] / 25.81 / (OD[i] - ID[i])

    # ECD = equivalent circulating density (lbs/gal)
    ecd = (apd / 0.052 / TVD) + RHO

    return apd, ecd


def dspd(Q, RHO, ID, THETA_600, THETA_300, L):
    for i in range(L):
        # Vp = pipe fluid velocity (ft/sec)
        vp = 0.408 * Q / ID[i] ** 2
        # np = pipe flow behavior index (dimensionless)
        np = 3.32 * np.log(THETA_600 / THETA_300)
        # Kp = pipe consistency factor (poise)
        kp = 5.11 * THETA_600 / 1022 ** np
        # μep = pipe effective viscosity (cp)
        muep = 100 * kp * (96 * vp / ID[i]) ** (np - 1)
        # Rep = pipe Reynolds number (dimensionless)
        rep = 928 * vp * ID[i] * RHO / (muep * ((3 * np + 1) / (4 * np)) ** np)
        # ReL = the laminar to transitional flow Reynolds number (dimensionless)
        rel = 3470 - 1370 * np
        # ReT = the laminar to transitional flow Reynolds number (dimensionless)
        ret = 4270 - 1370 * np
        # fp = the pipe fanning friction factor (dimensionless)
        if rep < rel:
            fp = FluidFlowType(16, rep, rel, ret, np).laminar_flow()
        elif ret < rep:
            fp = FluidFlowType(16, rep, rel, ret, np).transient_flow()
        elif rel < rep < ret:
            fp = FluidFlowType(16, rep, rel, ret, np).turbulent_flow()
        # Pp = the pipe pressure drop (psi)
        dspd += fp * vp ** 2 * RHO * L[i] / (25.81 * ID[i])

    return dspd


""" Needs confirmation!
def dcpd(Q, RHO, ID, THETA_600, THETA_300, L):
    for i in range(L):
        # Vc = collar fluid velocity (ft/sec)
        vc = 0.408 * Q / ID[i] ** 2
        # nc = collar flow behavior index (dimensionless)
        nc = 3.32 * np.log(THETA_600 / THETA_300)
        # Kc = collar consistency factor (poise)
        kc = 5.11 * THETA_600 / 1022 ** nc
        # μec = collar effective viscosity (cp)
        muec = 100 * kc * (96 * vc / ID[i]) ** (nc - 1)
        # Rec = collar Reynolds number (dimensionless)
        rec = 928 * vc * ID[i] * RHO / (muec * ((3 * nc + 1) / (4 * nc)) ** nc)
        # ReL = the laminar to transitional flow Reynolds number (dimensionless)
        rel = 3470 - 1370 * nc
        # ReT = the laminar to transitional flow Reynolds number (dimensionless)
        ret = 4270 - 1370 * nc
        # fp = the collar fanning friction factor (dimensionless)
        if rec < rel:
            fc = FluidFlowType(16, rec, rel, ret, nc).laminar_flow()
        elif ret < rec:
            fc = FluidFlowType(16, rec, rel, ret, nc).transient_flow()
        elif rel < rec < ret:
            fc = FluidFlowType(16, rec, rel, ret, nc).turbulent_flow()
        # Pc = the collar pressure drop (psi)
        dcpd += fc * vc ** 2 * RHO * L[i] / (25.81 * ID[i])

    return dcpd
"""


def hsi_dpb_jif(Q, MW, J1, J2, J3, B):
    # ΔPb (Pressure Loss at Bit)
    dpb = 156.5 * Q ** 2 * MW / (J1 ** 2 + J2 ** 2 + J3 ** 2) ** 2
    # HHPb (Hydraulic Horsepower at the Bit)
    hhp = Q * dpb / 1714

    # HHPb/in2 (Hydraulic Horsepower per Square Inch of the Bit Area)
    # Optimized Drilling Range: 2.5 - 5.0
    hsi = hhp * 1.27 / B ** 2

    # Vn (Bit Nozzle Velocity)
    vn = 417 * Q / (J1 ** 2 + J2 ** 2 + J3 ** 2)

    # I.F. (Impact Force)
    jif = vn * Q * MW / 1930

    return hsi, dpb, jif


# Not Standpipe Pressure
def sppd(Q, THETA_600, THETA_300, THETA_3, RHO, TVD, OD, ID, L, MW, J1, J2, J3, B):
    # PaT = total annular pressure loss (psi)
    apd = apd(Q, THETA_300, THETA_3, RHO, TVD, OD, ID, L)[0]
    # PpT = total drill string pressure loss (psi)
    dspd = dspd(Q, RHO, ID, THETA_600, THETA_300, L)
    # PcT = total drill collar pressure loss (psi)
    # dcpd = dcpd(Q, RHO, ID, THETA_600, THETA_300, L)
    # PB = the bit pressure loss (psi)
    dpb = hsi_dpb_jif(Q, MW, J1, J2, J3, B)[1]

    # PMAX = maximum standpipe pressure (psi)
    # return apd + dspd + dcpd + dpb
    return apd + dspd + dpb


main()
