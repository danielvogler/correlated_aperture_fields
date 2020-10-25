### 
### Daniel Vogler
### generate aperture fields with given stochastic properties
###

library(gstat)
library(sp)
library(tidyr)

date <- Sys.Date()
print("generating aperture field with gstat package\n")

### conversion factor micron to m
conversionMicronToM = 1e-6
### aperture
apertureMeanMicron = 200
apertureMean = apertureMeanMicron*conversionMicronToM

### variance
apertureVarianceMicron = 50
apertureVariance = apertureVarianceMicron*conversionMicronToM

### correlation length
apertureCorrelationLength = 10

### nugget
apertureNugget = 0
maximumRealizations = 100

# field dimensions
minX = -50
maxX = 50

minY = -50
maxY = 50

# spacing
spacingX = 1.0
spacingY = 1.0

epsilon = 1e-6

### print aperture values
cat("aperture mean [m]:", apertureMean, "\n")
cat("aperture variance [m]:", apertureVariance, "\n")
cat("aperture correlation length [m]:", apertureCorrelationLength, "\n")
cat("dimensions in x-direction [m]: min - ", minX, ",max - ", maxX, "\n")
cat("dimensions in y-direction [m]: min - ", minY, ",max - ", maxY, "\n")
cat("spacing x-direction [m]:", spacingX, "\n")
cat("spacing y-direction [m]:", spacingY, "\n")

# create structure
xy <- expand.grid( seq(minX, maxX, spacingX), seq(minY, maxY, spacingY) )
names(xy) <- c("x","y")

# define the gstat object (spatial model)
g.dummy <- gstat(formula=log(z)~1, locations=~x+y, dummy=T, beta=0,
	model=vgm(1,'Exp',apertureCorrelationLength,nugget=apertureNugget), nmax=maximumRealizations)

# make four simulations based on the stat object
yy <- predict(g.dummy, newdata=xy, nsim=1)

# show one realization
gridded(yy) = ~x+y

names(yy)
names(xy)

print(yy)
print(mean(yy@data))

# strings to save plot and figure
filenameString = paste("./reservoir_aperture_data_-_ACL-", apertureCorrelationLength, "m_AM-", apertureMeanMicron, "um_AV-", apertureVarianceMicron, "um_", date, ".txt", sep = "")
figurenameString = paste("./reservoir_aperture_data_-_ACL-", apertureCorrelationLength, "m_AM-", apertureMeanMicron, "um_AV-", apertureVarianceMicron, "um_", date, ".png", sep = "")

# show all four simulations:
tmpplot <- spplot(yy, scales=list(draw=T), colorkey=list(labels=list(at=c(0.0,0.00025,0.0005,0.00075,0.001))))
png(filename=figurenameString)
print(tmpplot)

# save table
write.table(yy, file=filenameString, row.names=FALSE, col.names=FALSE)
