import matplotlib.pyplot as plt
import numpy as np

''' Plot.py

Short script for Plotting the output from SpatialProfile if used in cylindrical mode.
	
Can be used to plot all SpatialProfile output files that are written in standard format.
Numer of skip_header lines might need to be adjusted if other / changed headers are written.
To use this plotting script:
	- Change the fileName to the proper name
	- set skip_header number for your specific header lines
	- set H to the size of the domain in Y-axis
	- set R to the smaller size of (X-axis,Z-axis)
	- set max_res to the number of bins you used in R/H direction (currently only same number supported, could be easily implemented)
	- set the view parameter y_min, y_max, R2_max to bin the bin numbers you want to plot (0 to max_res)
	- set levels to the number of contours desired
	- set datasets to the number of accumulated datasets, should be found in header or config.xml
'''



# FileName
fileName = 'profileName.NDpr'

# Read File
my_data = np.genfromtxt(fileName, skip_header=7, usecols=range(1, 800))

# Size of Domain in Y-Direction
H = 80
# Size of Domain in R-Direction from Center of XZ-Plane
R = 40
# Bin number
max_res = 800

# View in bin units
y_min = 30
y_max = 250
R2_max = 400

# Levels for Contour plot
levels = 100
# Number of accumulated timesteps (optional, set 1 to disregard)
datasets = 50

R_max = int(np.floor(np.sqrt(R2_max)))


# get R coordinates from R2, y is already linear
try:
	plt.contourf([np.sqrt(i) for i in range(R2_max-1)], [j for j in range(y_min, y_max)], my_data[y_min:y_max, 0:R2_max]/datasets, levels)
except:
	plt.contourf([np.sqrt(i) for i in range(R2_max)], [j for j in range(y_min, y_max)], my_data[y_min:y_max, 0:R2_max]/datasets, levels)

cb = plt.colorbar(ticks=[])

# x levels r^2
plt.xticks([i for i in range(0, R_max, int(R_max/5))], [np.round(i*(R/np.sqrt(max_res)),1) for i in range(0, R_max, int(R_max/5))])
# y levels are linear
plt.yticks([j for j in range(y_min, y_max, 40)], [j*(H/max_res) for j in range(y_min, y_max, 40)])
plt.axes().set_aspect("auto")

plt.show()