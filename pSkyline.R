# Spatial pSkyline

library(oce)
localPskylineFolder = "/home/welder/Documentos/projetos/pskyline"
source(paste(localPskylineFolder, "sphericalDistance.R", sep = "/"))




pSkyline <- function (tab, latCol, lonCol, meanError, sdError, p = 0.5, 
                      refLat = NULL, refLon = NULL, degreeNotation = TRUE, 
                      pdfSpatial = "normal") {
  # Default of reference point
  if (is.null(refLat)) {
    refLat = tab[1, latCol]
    refLon = tab[1, lonCol]
  }
    
  # Cast coordinates to numeric
  tab[, latCol] = as.numeric(as.character((tab[, latCol])))
  tab[, lonCol] = as.numeric(as.character((tab[, lonCol])))
  
  P = c()   # set of Pareto efficient points 
  countOfParetoEfficientPoints = 0
  probParetoEfficient = numeric(nrow(tab))
  for (i in 1:nrow(tab)) {
    probBeDominatedByThis = numeric(nrow(tab)) # ignore position i (It is itself)
    for (j in (1:nrow(tab))) {
      probBeDominatedByThis[j] = 
        probDominance(pointEvaluated = list(lat = tab[i, latCol], lon = tab[i, lonCol]), 
                      concurrent = list(lat = tab[j, latCol], lon = tab[j, lonCol]), 
                      refLat, refLon, meanError, sdError, pdfSpatial)
    }
    # compute the prob of being Pareto efficient
    probParetoEfficient[i] = prod(1 - probBeDominatedByThis[-i])
    # cat("skyline prob of object "); cat(i); cat(": "); cat(probParetoEfficient[i]); cat("\n")
    if (probParetoEfficient[i] >= p) {
      countOfParetoEfficientPoints = countOfParetoEfficientPoints + 1
      P[countOfParetoEfficientPoints] = i
    }
  }
  
  cat("Sorted skyline probabities: "); 
  cat(sort(probParetoEfficient[probParetoEfficient >= p], decreasing = T)); cat("\n")
  
  # Ploting pskylines points
  plot(tab[, latCol] ~ tab[, lonCol], pch = 20, main = "p-skyline (red points)",  
       xlab = "longitude", ylab = "latitude")  # all points
  points(refLat ~ refLon, pch = 8, cex = 2, col = 2)
  points(tab[P, latCol] ~ tab[P, lonCol], pch = 20, col = 2)
  
  return(P)
}


deriveDistFromRefPoints <- function (tab, latCol, lonCol, refLat = NULL, refLon = NULL) {
  if (is.null(refLat)) {
    derivedDistance = matrix(0, nrow(tab), 1)
    for (i in 1:nrow(tab)) {
      derivedDistance[i, 1] = 
        sphericalDist(tab[i, latCol], tab[i, lonCol], tab[1, latCol], tab[1, lonCol])  
    }
  } else {
    derivedDistance = matrix(0, nrow(tab), length(refLat))
    for (j in 1:length(refLat)) {
      for (i in 1:nrow(tab)) {
        derivedDistance[i, j] = 
          sphericalDist(tab[i, latCol], tab[i, lonCol], refLat[j], refLon[j])  
      }
    }
  }
  return (derivedDistance) 
}

probDominance <- function(pointEvaluated, concurrent, refLat, refLon, 
                          meanError, sdError, pdfSpatial) {
  # simulating real positions for 'pointEvaluated' and 'concurrent' points
  simulatedDistFromEvaluatedToRefPoint = simulateDistToRefPoint(pointEvaluated, n = 20,
                                                            refLat, refLon, meanError, 
                                                            sdError, pdfSpatial)
  simulatedDistFromConcurrentToRefPoint = simulateDistToRefPoint(concurrent, n = 20, 
                                                            refLat, refLon, meanError, 
                                                            sdError, pdfSpatial)
  
  # prob = mean(simulatedDistFromConcurrentToRefPoint <  # até o dia 20/12
  #               simulatedDistFromEvaluatedToRefPoint)
  prob = 1
  for (j in 1:length(refLat)) {
    prob = prob * mean(simulatedDistFromConcurrentToRefPoint[, j] <  
                         simulatedDistFromEvaluatedToRefPoint[, j])
  }
  return (prob)
}

probSkyline <- function(u, S) {
  return (prob)
}

simulateDistToRefPoint <- function(databasePoint, n = 20, refLat, refLon, meanError, 
                               sdError, distribution = "normal") {
  
  # simulating the size of the positional error
  if (distribution %in% c("normal", "norm", "n"))
    sizePositionalError = rnorm(n, mean = meanError, sd = sdError)
  if (distribution %in% c("exponencial", "expo", "exp", "e"))
    sizePositionalError = rexp(n, rate = 1 / meanError)
  if (distribution %in% c("chisquare", "chisq", "c"))
    sizePositionalError = rchisq(n, df = meanError)
  
  angle = runif(n, 0, 2*pi)
  errorInLon = sizePositionalError * cos(angle)
  errorInLat = sizePositionalError * sin(angle)
  new_lat = databasePoint$lat + metersToDegrees(errorInLat)
  new_lon = databasePoint$lon + metersToDegrees(errorInLon)
  
  # distances = numeric(n) # até o dia 20/12
  distances = matrix(0, n, length(refLat))
  for (i in 1:n) {
    for (j in 1:length(refLat)) {
      # distances[i] = sphericalDist(new_lat[i], new_lon[i], refLat, refLon) # até o dia 20/12
      distances[i, j] = sphericalDist(new_lat[i], new_lon[i], refLat[j], refLon[j])
    }
  }
  return (distances)
}

metersToDegrees <- function(x) {
  # consider that each degree possess 111.111 km
  return(x/111111)
}

# Example
pathSchool = paste(localPskylineFolder, "sampleData_36gynSchools", sep = "/")
tab = read.csv2(pathSchool)
latCol = 2
lonCol = 3
refLat = c(tab[1, 2], tab[2, 2])
refLon = c(tab[1, 3], tab[2, 3])
refLat = refLon = NULL
degreeNotation = TRUE
p = 0.5

meanError = 60
sdError = 10
pdfSpatial = "normal"


pSkyline(tab, latCol, lonCol, meanError = 100, sdError = 1000, p = 0.50, 
         refLat = c(-16.68, -16.71, -16.695), refLon = c(-49.21, -49.275, -49.35), degreeNotation = TRUE, 
         pdfSpatial = "normal")
