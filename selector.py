#Algorithms to hold the information that is required to run the main program
import numpy as np #Dependency => pip install numpy
import sys #Dependency => Installed by default

#Component = ["Pc" - psi, "Tc" - R, "Vc" - ft3/lbmole, "OMEGA" - 1]
CONSTANTS = np.array([["methane", 666.4, -116.67+460.0, 0.0958*16.043, 0.010],
                      ["ethane", 706.5, 89.92+460.0, 0.0783*30.070, 0.098],
                      ["propane", 616.0, 206.06+460.0, 0.0727*44.097, 0.153],
                      ["n-butane", 550.6, 305.62+460, 0.0703*58.123, 0.199],
                      ["n-pentane", 488.5, 385.8+460, 0.0875*72.150, 0.251]],
                      dtype=object)

#Equation of State = ["Number", "Name of the model"]
EOS = np.array([["1", "srk_eos"],
                ["2", "vdw_eos"],
                ["3", "pr_eos"]],
                dtype=object)

#Temperature = ["Name of the component", "T" - F]
TEMP = np.array([["methane", -220.0, -180.0, -150.0, -135.0, -117.0, -105.0],
                 ["ethane", -10.0, 20.0, 40.0, 60.0, 89.9, 120.0],
                 ["propane", 60.0, 100.0, 140.0, 180.0, 206.6, 220.0],
                 ["n-butane", 140.0, 180.0, 220.0, 260.0, 205.6, 330.0],
                 ["n-pentane", 160.0, 230.0, 280.0, 340.0, 386.6, 430.0]],
                 dtype=object)

#Estimated Vapor Pressure = ["Name of the component", "P" - psi]
PSAT = np.array([["methane", 63.44, 189.37, 359.89, 478.50],
                 ["ethane", 187.24, 296.68, 387.26, 499.20],
                 ["propane", 108.83, 191.35, 311.86, 478.87],
                 ["n-butane", 93.84, 155.73, 243.26, 363.21],
                 ["n-pentane", 42.93, 108.54, 187.01, 330.54]], 
                 dtype=object)

def constants(component):
    while True:
        try:
            row = np.where(CONSTANTS[:,0] == component)[0]
        except ValueError:
            print("Invalid entry. Terminating program...")
            sys.exit()
        if len(row) == 0:
            print("Invalid entry. Terminating program...")
            sys.exit()
        else:
            break

    return CONSTANTS[row.tolist()[0], 1:]

def eos(model):
    while True:
        try:
            row = np.where(EOS[:,0] == model)[0]
        except IndexError:
            print("Invalid entry. Terminating program...")
            sys.exit()
        else:
            break

    return EOS[row.tolist()[0], 1]

def temp(component):
    row = np.where(TEMP[:,0] == component)[0]
    return TEMP[row.tolist()[0], 1:].astype(float)

def psat(component):
    row = np.where(PSAT[:,0] == component)[0]
    return PSAT[row.tolist()[0], 1:].astype(float)