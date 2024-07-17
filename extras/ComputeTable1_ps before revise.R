library(DdiPpiClo)
library(dplyr)



# Please upload "FunctionsForTable1.R" file in the "extras" folder of the DdiPpiClo package
source("./extras/functionsForTable1.R")  # Please set the path to the "extras" folder of the DdiPpiClo package

#### Need to define #### 
outputFolder <- " " # Please set the path to the "extras/result" folder of the DdiPpiClo package
databaseId <- " " # Please write your databaseId

################################

minCellCount <- 5
maxCores <- parallel::detectCores()
exportFolder <- file.path(outputFolder, "export")
balanceFolder <- file.path(outputFolder, "balanceRevise")

if(!file.exists(balanceFolder)) {
  dir.create(balanceFolder)
}

cmOutputFolder <- file.path(outputFolder, "cmOutput")
results <- readRDS(file.path(cmOutputFolder, "outcomeModelReference.rds"))

outcomesOfInterest <- getOutcomesOfInterest()
  

subset <- results[results$outcomeId %in% outcomesOfInterest,]
subset <- subset[subset$strataFile != "", ]

subset <- split(subset, seq(nrow(subset)))
lapply(subset, computeCovariateBalanceRevise, cmOutputFolder = cmOutputFolder, balanceFolder = balanceFolder)


ParallelLogger::logInfo("Exporting diagnostics")
ParallelLogger::logInfo("- covariate_balance table")
fileName <- file.path(exportFolder, "covariate_balance_revise.csv")
if (file.exists(fileName)) {
  unlink(fileName)
}
first <- TRUE
balanceFolder <- file.path(outputFolder, "balanceRevise")
files <- list.files(balanceFolder, pattern = "bal_.*.rds", full.names = TRUE)
pb <- txtProgressBar(style = 3)
for (i in 1:length(files)) {
  ids <- gsub("^.*bal_t", "", files[i])
  targetId <- as.numeric(gsub("_c.*", "", ids))
  ids <- gsub("^.*_c", "", ids)
  comparatorId <- as.numeric(gsub("_[aso].*$", "", ids))
  if (grepl("_s", ids)) {
    subgroupId <- as.numeric(gsub("^.*_s", "", gsub("_a[0-9]*.rds", "", ids)))
  } else {
    subgroupId <- NA
  }
  if (grepl("_o", ids)) {
    outcomeId <- as.numeric(gsub("^.*_o", "", gsub("_a[0-9]*.rds", "", ids)))
  } else {
    outcomeId <- NA
  }
  ids <- gsub("^.*_a", "", ids)
  analysisId <- as.numeric(gsub(".rds", "", ids))
  balance <- readRDS(files[i])
  inferredTargetBeforeSize <- mean(balance$beforeMatchingSumTarget/balance$beforeMatchingMeanTarget,
                                   na.rm = TRUE)
  inferredComparatorBeforeSize <- mean(balance$beforeMatchingSumComparator/balance$beforeMatchingMeanComparator,
                                       na.rm = TRUE)
  inferredTargetAfterSize <- mean(balance$afterMatchingSumTarget/balance$afterMatchingMeanTarget,
                                  na.rm = TRUE)
  inferredComparatorAfterSize <- mean(balance$afterMatchingSumComparator/balance$afterMatchingMeanComparator,
                                      na.rm = TRUE)
  
  balance$databaseId <- databaseId
  balance$targetId <- targetId
  balance$comparatorId <- comparatorId
  balance$outcomeId <- outcomeId
  balance$analysisId <- analysisId
  balance$interactionCovariateId <- subgroupId
  balance <- balance[, c("databaseId",
                         "targetId",
                         "comparatorId",
                         "outcomeId",
                         "analysisId",
                         "interactionCovariateId",
                         "covariateId",
                         "beforeMatchingSumTarget",
                         "beforeMatchingMeanTarget",
                         "beforeMatchingSumComparator",
                         "beforeMatchingMeanComparator",
                         "beforeMatchingStdDiff",
                         "afterMatchingSumTarget",
                         "afterMatchingMeanTarget",
                         "afterMatchingSumComparator",
                         "afterMatchingMeanComparator",
                         "afterMatchingStdDiff")]
  colnames(balance) <- c("databaseId",
                         "targetId",
                         "comparatorId",
                         "outcomeId",
                         "analysisId",
                         "interactionCovariateId",
                         "covariateId",
                         "targetSumBefore",
                         "targetMeanBefore",
                         "comparatorSumBefore",
                         "comparatorMeanBefore",
                         "stdDiffBefore",
                         "targetSumAfter",
                         "targetMeanAfter",
                         "comparatorSumAfter",
                         "comparatorMeanAfter",
                         "stdDiffAfter")
  
  balance$targetMeanBefore[is.na(balance$targetMeanBefore)] <- 0
  balance$comparatorMeanBefore[is.na(balance$comparatorMeanBefore)] <- 0
  balance$stdDiffBefore <- round(balance$stdDiffBefore, 3)
  balance$targetMeanAfter[is.na(balance$targetMeanAfter)] <- 0
  balance$comparatorMeanAfter[is.na(balance$comparatorMeanAfter)] <- 0
  balance$stdDiffAfter <- round(balance$stdDiffAfter, 3)
  balance <- enforceMinCellValue(balance,
                                 "targetMeanBefore",
                                 minCellCount/inferredTargetBeforeSize,
                                 TRUE)
  balance <- enforceMinCellValue(balance,
                                 "comparatorMeanBefore",
                                 minCellCount/inferredComparatorBeforeSize,
                                 TRUE)
  balance <- enforceMinCellValue(balance,
                                 "targetMeanAfter",
                                 minCellCount/inferredTargetAfterSize,
                                 TRUE)
  balance <- enforceMinCellValue(balance,
                                 "comparatorMeanAfter",
                                 minCellCount/inferredComparatorAfterSize,
                                 TRUE)
  balance$targetMeanBefore <- round(balance$targetMeanBefore, 3)
  balance$comparatorMeanBefore <- round(balance$comparatorMeanBefore, 3)
  balance$targetMeanAfter <- round(balance$targetMeanAfter, 3)
  balance$comparatorMeanAfter <- round(balance$comparatorMeanAfter, 3)
  balance <- balance[balance$targetMeanBefore != 0 & balance$comparatorMeanBefore != 0 & balance$targetMeanAfter !=
                       0 & balance$comparatorMeanAfter != 0 & balance$stdDiffBefore != 0 & balance$stdDiffAfter !=
                       0, ]
  balance <- balance[!is.na(balance$targetId), ]
  colnames(balance) <- SqlRender::camelCaseToSnakeCase(colnames(balance))
  write.table(x = balance,
              file = fileName,
              row.names = FALSE,
              col.names = first,
              sep = ",",
              dec = ".",
              qmethod = "double",
              append = !first)
  first <- FALSE
  setTxtProgressBar(pb, i/length(files))
}
close(pb)
