import numpy as np


def main():
    # Flow rate, gal / min
    Q = 160
    # Viscometer @ 600 RPM
    THETA_600 = 38
    # Viscometer @ 300 RPM
    THETA_300 = 28
    # Viscometer @ 3 RPM
    THETA_3 = 18
    # Density, lbm / gal
    RHO = 12
    # True vertical depth, ft
    TVD = 7277
    # Outer diameter of each section, in
    OD = np.array([6.5, 2.875])
    # Inner diameter of each section, in
    ID = np.array([3, 2.151])
    # Length of each section, ft
    L = np.array([4000, 3426])
    # Mud weight
    MW = 13
    # Bit size, in
    B = 8.375
    # Jet nozzle optimization constant
    # HHP = 0.65, IF = 0.48, balanced = 0.59
    C = 0.59
    # Jet nozzle amount
    N = 3
    # Surface equipment with: ([40, 55, 5, 40], [3.5, 2.5, 2.25, 3.25]) @ 160 GPM
    # Surface equipment pressure loss (psi)
    SEPD = 8

    print(
        f"Standpipe Pressure (psi): {ons(Q, C, RHO, THETA_600, THETA_300, THETA_3, TVD, OD, ID, L, SEPD, N)[1]:.2f}")
    print(
        f"Equivalent Circulating Density (lbs/gal): {apd_ecd(Q, THETA_300, THETA_3, RHO, TVD, OD, ID, L)[1]:.2f}")
    print(
        f"Bit Pressure Drop (psi): {hsi_jif_pbd(Q, MW, B, C, RHO, THETA_600, THETA_300, THETA_3, TVD, OD, ID, L, SEPD, N)[2]:.2f}")
    print(
        f"Horsepower per square inch (HHPb/in2): {hsi_jif_pbd(Q, MW, B, C, RHO, THETA_600, THETA_300, THETA_3, TVD, OD, ID, L, SEPD, N)[0]:.2f}")
    print(
        f"Jet Impact Force (lbf): {hsi_jif_pbd(Q, MW, B, C, RHO, THETA_600, THETA_300, THETA_3, TVD, OD, ID, L, SEPD, N)[1]:.2f}")
    print(
        f"Optimum Nozzle Size (1/32 in): {ons(Q, C, RHO, THETA_600, THETA_300, THETA_3, TVD, OD, ID, L, SEPD, N)[0][-1]:.2f}")

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
        return ((np.log10(self.na) + 3.93) / 50) / self.re ** ((1.75 - np.log10(self.na)) / 7)


def apd_ecd(Q, THETA_300, THETA_3, RHO, TVD, OD, ID, L):
    apd = 0

    for i in range(len(L)):
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
    dspd = 0

    for i in range(len(L)):
        # Vp = pipe fluid velocity (ft/sec)
        vp = 0.408 * Q / ID[i] ** 2
        # np = pipe flow behavior index (dimensionless)
        np_val = 3.32 * np.log(THETA_600 / THETA_300)
        # Kp = pipe consistency factor (poise)
        kp = 5.11 * THETA_600 / 1022 ** np_val
        # μep = pipe effective viscosity (cp)
        muep = 100 * kp * (96 * vp / ID[i]) ** (np_val - 1)
        # Rep = pipe Reynolds number (dimensionless)
        rep = 928 * vp * ID[i] * RHO / \
            (muep * ((3 * np_val + 1) / (4 * np_val)) ** np_val)
        # ReL = the laminar to transitional flow Reynolds number (dimensionless)
        rel = 3470 - 1370 * np_val
        # ReT = the laminar to transitional flow Reynolds number (dimensionless)
        ret = 4270 - 1370 * np_val
        # fp = the pipe fanning friction factor (dimensionless)
        if rep < rel:
            fp = FluidFlowType(16, rep, rel, ret, np_val).laminar_flow()
        elif ret < rep:
            fp = FluidFlowType(16, rep, rel, ret, np_val).turbulent_flow()
        elif rel < rep < ret:
            fp = FluidFlowType(16, rep, rel, ret, np_val).transient_flow()
        # Pp = the pipe pressure drop (psi)
        dspd += fp * vp ** 2 * RHO * L[i] / (25.81 * ID[i])

    return dspd


# Jet optimization
def ons(Q, C, RHO, THETA_600, THETA_300, THETA_3, TVD, OD, ID, L, SEPD, N):
    apd_val = apd_ecd(Q, THETA_300, THETA_3, RHO, TVD, OD, ID, L)[0]
    dspd_val = dspd(Q, RHO, ID, THETA_600, THETA_300, L)

    # Bit pressure drop
    # Good hydraulics normally result when
    # between 50% to 65% of surface pressure is utilized at the bit
    # We used 60%
    dpb_val = (apd_val + dspd_val + SEPD) * (6 / 4)

    at = Q / (2.96 * (1238.5 * C * dpb_val / RHO) ** 0.5)
    # Nozzle sizes
    j = np.zeros(N)

    for i in range(N):
        j[i] = ((1303.797 / (N - i)) * (at - (np.sum(j ** 2) / 1303.797))) ** 0.5

    return j, dpb_val / 0.6


# HSI and JIF
def hsi_jif_pbd(Q, MW, B, C, RHO, THETA_600, THETA_300, THETA_3, TVD, OD, ID, L, SEPD, N):

    # Jet nozzle sizes
    j_val = ons(Q, C, RHO, THETA_600, THETA_300,
                THETA_3, TVD, OD, ID, L, SEPD, N)[0]

    # ΔPb (Pressure Loss at Bit)
    pbd_val = 156.5 * Q ** 2 * MW / (np.sum(j_val ** 2)) ** 2

    # HHPb (Hydraulic Horsepower at the Bit)
    hhp = Q * pbd_val / 1714

    # HHPb/in2 (Hydraulic Horsepower per Square Inch of the Bit Area)
    # Optimized Drilling Range: 2.5 - 5.0
    hsi = hhp * 1.27 / B ** 2

    # Vn (Bit Nozzle Velocity)
    vn = 417 * Q / (np.sum(j_val ** 2))

    # I.F. (Impact Force)
    jif = vn * Q * MW / 1930

    return hsi, jif, pbd_val


main()
