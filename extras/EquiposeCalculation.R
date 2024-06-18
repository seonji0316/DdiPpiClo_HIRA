.libPaths("  ") #Please set the path where the DdiPpiClo package library is located. ex) C:/DdiPpiClo/renv/library/R-4.1/x86_64-w64-mingw32

outputFolder <- "  " # Please set the path to the "extras/result" folder of the DdiPpiClo package
cmOutputFolder <- file.path(outputFolder, "cmOutput")
omr <- readRDS(file.path(cmOutputFolder, "outcomeModelReference.rds"))

getOutcomesOfInterest <- function() {
  pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "DdiPpiClo")
  tcosOfInterest <- read.csv(pathToCsv, stringsAsFactors = FALSE) 
  outcomeIds <- as.character(tcosOfInterest$outcomeIds)
  outcomeIds <- do.call("c", (strsplit(outcomeIds, split = ";")))
  outcomeIds <- unique(as.numeric(outcomeIds))
  return(outcomeIds)
}

computePreferenceScore <- function (data, unfilteredData = NULL) {
  if (is.null(unfilteredData)) {
    proportion <- sum(data$treatment)/nrow(data)
  }
  else {
    proportion <- sum(unfilteredData$treatment)/nrow(unfilteredData)
  }
  propensityScore <- data$propensityScore
  propensityScore[propensityScore > 0.9999999] <- 0.9999999
  x <- exp(log(propensityScore/(1 - propensityScore)) - log(proportion/(1 - proportion)))
  data$preferenceScore <- x/(x + 1)
  return(data)
}

outcomesOfInterest <- getOutcomesOfInterest()

subset <- omr[omr$outcomeId %in% outcomesOfInterest,]
subset <- subset[subset$strataFile != "", ]
subset <- split(subset, seq(nrow(subset)))

computeEquipoise <- function(row) {
  ps <- readRDS(file.path(outputFolder, "cmOutput", row$psFile))
  ps <- computePreferenceScore(ps)
  auc <- CohortMethod::computePsAuc(ps)
  equipoise <- mean(ps$preferenceScore >= 0.3 & ps$preferenceScore <= 0.7)
  results <- list(targetId = row$targetId, 
                  comparatorId = row$comparatorId, 
                  outcomeId = row$outcomeId,
                  analysisId = row$analysisId,
                  auc = auc, 
                  equipoise = equipoise
  )
  return(results)
}

results <- lapply(subset, computeEquipoise)
results <- do.call(rbind.data.frame, results)

write.csv(results, file.path(outputFolder, "export", "ps_result.csv"), row.names = F)
results$databaseId <- "HIRA"