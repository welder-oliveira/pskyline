# Testes do projeto pSkyline
localPskylineFolder = "/home/welder/Documentos/projetos/pskyline"
source(paste(localPskylineFolder, "sphericalDistance.R", sep = "/"))
source(paste(localPskylineFolder, "pSkyline.R", sep = "/"))


probDominanceTest <- function() {
  consistencyInFindingTheProbableNearestPointToReference = FALSE
  distFromAPointWithErrorToItselfTest = FALSE
  distFromAPointWithErrorToAnotherFarAwayTest = FALSE
  sphericalDistTest = FALSE
  
  resultconsistencyInFindingTheProbableNearestPointToReference = 
                probDominance(pointEvaluated = list(lat = -24.00, lon = -49.00), 
                concurrent = list(lat = -24.00, lon = -50.00), 
                refLat = -24.01, refLon = -50.01, meanError = 100, sdError = 100, 
                pdfSpatial = "normal")
  
  distFromAPointWithErrorToItself = simulateDistToRefPoint(list(lat = -24.00, lon = -49.00), 
                          n = 20, refLat = -24.00, refLon = -49.00, meanError = 1, 
                          sdError = 0.01, pdfSpatial)

  distFromAPointWithErrorToAnotherFarAway = simulateDistToRefPoint(list(lat = -24.00, lon = -49.00), 
                          n = 20, refLat = -24.00, refLon = -50.00, meanError = 1, 
                          sdError = 0.01, pdfSpatial)
  
  geodesicDist = sphericalDist(lat1 = -24.00, long1 = -49.00, lat2 = -25.00, long2 = -49.00)
  
  if (resultconsistencyInFindingTheProbableNearestPointToReference == 1)
    consistencyInFindingTheProbableNearestPointToReference = TRUE

  if ( mean(distFromAPointWithErrorToItself) < 1) 
    distFromAPointWithErrorToItselfTest = TRUE
  
  if (mean(distFromAPointWithErrorToAnotherFarAway) > 1e5 && 
      mean(distFromAPointWithErrorToAnotherFarAway < 2e5))
    distFromAPointWithErrorToAnotherFarAwayTest = TRUE
  
  if (geodesicDist > 1e5 && geodesicDist < 1e5 + 2e4)
    sphericalDistTest = TRUE
  
  return(list(consistencyInFindingTheProbableNearestPointToReference = 
                consistencyInFindingTheProbableNearestPointToReference, 
              distFromAPointWithErrorToItselfTest = distFromAPointWithErrorToItselfTest, 
              distFromAPointWithErrorToAnotherFarAwayTest = distFromAPointWithErrorToAnotherFarAwayTest, 
              sphericalDistTest = sphericalDistTest
              ))
}


probDominanceTest()
