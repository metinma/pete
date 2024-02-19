import matplotlib.pyplot as plt
import numpy as np

Pressure = np.array([29.6, 29.5, 29.4, 29.3, 29.2, 29.1])
time = np.array([0, 2, 4, 6, 8, 10])
temp = np.array([40.0, 35.0, 30.0, 25.0, 20.0, 15.0])
ln_Pressure = np.log(Pressure)
over_temp = 1 / temp

R = 8.3144598

slope = (ln_Pressure[-1] - ln_Pressure[0]) / (over_temp[-1] - over_temp[0])

hv = -(slope * R) 
print(hv)

fig, ax = plt.subplots(1,3)
FSIZE = 5
fig.set_size_inches(FSIZE*3, FSIZE)

ax[0].plot(time, Pressure, marker="o")
ax[0].set_xlabel("Time, min")
ax[0].set_ylabel("Pressure, kPa")

ax[1].plot(temp, Pressure, marker="o")
ax[1].invert_xaxis()
ax[1].set_xlabel("Temperature, °C")
ax[1].set_ylabel("Pressure, kPa")

ax[2].plot(over_temp, ln_Pressure, marker="o")
ax[2].set_xlabel("1 / T")
ax[2].set_ylabel("ln P")

plt.show()