'''
Daniel Vogler
produce spatially correlated aperture fields
'''

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import csv
import time
import sys
import glob, os
import re
import math

## dd/mm/yyyy format
currentDate =  (time.strftime("%Y%m%d"))

print("System input variables")
print(str(sys.argv))

# read in filepath from command line
filename =  sys.argv[1]

# trim file path
filenameCore = re.sub('.txt', '', filename)

print("Saving file as %s: ", filenameCore)

apertureValue = []
xCoordinate = []
yCoordinate = []

### plot variables
fontSize = 15

saveFileFolder = './'
fieldName = 'normal'

# define measurement locations on sample
fileType = ['X', 'Y', 'Z']
# convert from meters to mm
conversionScale = 1
# maximum aperture
maxAperture = 0.001 #/conversionScale
# smallest contactOffset apertureValue
contactOffsetMinimum = 4e-6
# set number of contact points in contact
fractionContact = 0.1
# inlet radius
inletRadius = 3.0
# inletRadiusMaxAperture
inletRadiusMaxAperture = 0.33*3.0
# colormap
colormap = cm.gnuplot


with open(filename) as f:
	reader = csv.reader(f, delimiter=' ')
	for line in reader:
		apertureValue.append(float(line[0]))
		xCoordinate.append(float(line[1]))
		yCoordinate.append(float(line[2]))

apertureValue = [1.0*math.exp(x)+0 for x in apertureValue]
apertureValue = [x/np.max(apertureValue)*0.01 for x in apertureValue]

print("Mean aperture: %f \n" %np.mean(apertureValue) )
print("Variance aperture: %f \n" %np.var(apertureValue) )
print("Mean x-Coordinate: %f \n" %np.mean(xCoordinate) )
print("Mean y-Coordinate: %f \n" %np.mean(yCoordinate) )

### Contact creation
# number of points in contact for given percentContact
numberContact = round( fractionContact*np.size(apertureValue) )
print("\t Total number of points: %d" %np.size(apertureValue) )
print("\t Number of points in contact: %d" %numberContact )

# split array in contact and not in contact
splitValuesNoContact = np.partition(apertureValue, numberContact)[numberContact]

# find maximum apertureValue of points not in contact
maxValueNoContact = np.max(splitValuesNoContact)
print("\n \t Maximum global aperture before bringing in contact: %f" %np.max(apertureValue) )
print("\t Minimum global aperture before bringing in contact: %f" %np.min(apertureValue) )

print("\t Maximum local aperture: %f" %maxValueNoContact)

# move aperture values by maximum apertureValue for contact
#apertureValue = apertureValue - maxValueNoContact

# set values lower than contactOffsetMinimum to contactOffsetMinimum
apertureValue = [float(x) for x in apertureValue]
apertureValue = [contactOffsetMinimum if i < contactOffsetMinimum else i for i in apertureValue]
apertureValue = [maxAperture if i > maxAperture else i for i in apertureValue]

print("\n\t Maximum global aperture after bringing in contact: %f" %np.max(apertureValue) )
print("\t Minimum global aperture after bringing in contact: %f" %np.min(apertureValue) )
print("\n\t Mean global aperture after bringing in contact: %f" %np.mean(apertureValue) )
print("\t Variance global aperture after bringing in contact: %f" %np.var(apertureValue) )


# find nearest neighbour
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

coordLength = len(apertureValue)**0.5

maxReservoirGridZ = maxAperture

# number of grid points in x and y
nX = int(coordLength)
nY = int(coordLength)

minX = np.min(xCoordinate)
maxX = np.max(xCoordinate)
minY = np.min(yCoordinate)
maxY = np.max(yCoordinate)

# pick x and y reservoir locations
reservoirX = np.linspace(minX, maxX, num=nX)
reservoirY = np.linspace(minY, maxY, num=nY)

# initialize reservoir grid
reservoirGridX = [[0 for x in range(nX)] for x in range(nY)]
reservoirGridY = [[0 for x in range(nX)] for x in range(nY)]
reservoirGridZ = [[0 for x in range(nX)] for x in range(nY)]

# reshape vector to 2D array
reservoirGridX = np.reshape( xCoordinate, (-1, nX))
reservoirGridY = np.reshape( yCoordinate, (-1, nY))
reservoirGridZ = np.reshape( apertureValue, (-1, nY))

# total number of points
totalNumberOfPoints = len(reservoirX)*len(reservoirY)

# convert to float
reservoirGridX = np.array(reservoirGridX).astype(np.float)
reservoirGridY = np.array(reservoirGridY).astype(np.float)
reservoirGridZ = np.array(reservoirGridZ).astype(np.float)

# inlet location
inlet = [[-inletRadius, inletRadius],[-inletRadius, inletRadius]]

# find nearest neighbour
value = [[-1, -1],[-1, -1]]
index = [[-1, -1],[-1, -1]]

# findest nearest indeces to inlet
value[0][0], index[0][0] = find_nearest( reservoirX, inlet[0][0] )
value[0][1], index[0][1] = find_nearest( reservoirX, inlet[0][1] )
value[1][0], index[1][0] = find_nearest( reservoirY, inlet[1][0] )
value[1][1], index[1][1] = find_nearest( reservoirY, inlet[1][1] )

# replace values near inlet with high aperture
for i in range(index[0][0], index[0][1]):
	for j in range(index[1][0], index[1][1]):
		radius = (reservoirGridX[i][j]**2 + reservoirGridY[i][j]**2 )**(0.5)
		ratio = radius/inletRadius

		if ratio <= 0.4:

			reservoirGridZ[i][j] = maxAperture

		if (ratio < 0.4 and ratio <= 1.0):

			inletValue = (1.0 - ratio + inletRadiusMaxAperture)*maxAperture

			reservoirGridZ[i][j] = max(inletValue, reservoirGridZ[i][j])


### plot aperture distribution histogram
print("Plot aperture distribution \n")
pl.figure()
pl.hist(apertureValue, bins=100, range=[0.0, 0.0011])
pl.title("aperture distribution")
pl.xlabel("aperture [m]")
pl.ylabel("frequency [-]")

### plot aperture field
print("Plot reservoir aperture field\n")
pl.figure( num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k' )
# plot
pl.pcolormesh(reservoirGridX, reservoirGridY, reservoirGridZ, cmap=colormap, vmin=0, vmax=maxReservoirGridZ )
# limits and labels
pl.xlabel('x [m]', fontsize = fontSize)
pl.ylabel('y [m]', fontsize = fontSize)
pl.grid(b=True, which='major', color='lightgrey', linestyle='-')
pl.tick_params(axis='x', which='major', labelsize=fontSize)
pl.tick_params(axis='y', which='major', labelsize=fontSize)
# colorbar
cbar = pl.colorbar()
cbar.ax.tick_params(labelsize=fontSize)
cbar.ax.set_title('aperture [m]')
pl.savefig( str( saveFileFolder+filenameCore+fieldName+"_"+currentDate+".png"), bbox_inches='tight' )

pl.show()
