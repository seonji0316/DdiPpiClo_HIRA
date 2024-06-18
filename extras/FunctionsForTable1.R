getOutcomesOfInterest <- function() {
  pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "DdiPpiClo")
  tcosOfInterest <- read.csv(pathToCsv, stringsAsFactors = FALSE) 
  outcomeIds <- as.character(tcosOfInterest$outcomeIds)
  outcomeIds <- do.call("c", (strsplit(outcomeIds, split = ";")))
  outcomeIds <- unique(as.numeric(outcomeIds))
  return(outcomeIds)
}

computeMeansPerGroup <- function (cohorts, cohortMethodData) 
{
  hasStrata <- "stratumId" %in% colnames(cohorts)
  if (hasStrata) {
    stratumSize <- cohorts %>% group_by(.data$stratumId, 
                                        .data$treatment) %>% count() %>% ungroup()
  }
  if (hasStrata && any(stratumSize %>% pull(.data$n) > 1)) {
    w <- stratumSize %>% mutate(weight = 1/.data$n) %>% 
      inner_join(cohorts, by = c("stratumId", "treatment")) %>% 
      select(.data$rowId, .data$treatment, .data$weight)
    wSum <- w %>% group_by(.data$treatment) %>% summarize(wSum = sum(.data$weight, 
                                                                     na.rm = TRUE)) %>% ungroup()
    cohortMethodData$w <- w %>% inner_join(wSum, by = "treatment") %>% 
      mutate(weight = .data$weight/.data$wSum) %>% select(.data$rowId, 
                                                          .data$treatment, .data$weight)
    sumW <- 1
    result <- cohortMethodData$covariates %>% inner_join(cohortMethodData$w, 
                                                         by = c("rowId")) %>% group_by(.data$covariateId, 
                                                                                       .data$treatment) %>% summarise(sum = sum(.data$covariateValue, 
                                                                                                                                na.rm = TRUE), mean = sum(.data$weight * .data$covariateValue, 
                                                                                                                                                          na.rm = TRUE), sumSqr = sum(.data$weight * .data$covariateValue^2, 
                                                                                                                                                                                      na.rm = TRUE), sumWSqr = sum(.data$weight^2, na.rm = TRUE)) %>% 
      mutate(sd = sqrt(abs(.data$sumSqr - .data$mean^2) * 
                         sumW/(sumW^2 - .data$sumWSqr))) %>% ungroup() %>% 
      select(.data$covariateId, .data$treatment, .data$sum, 
             .data$mean, .data$sd) %>% collect()
    cohortMethodData$w <- NULL
  }
  else {
    cohortCounts <- cohorts %>% group_by(.data$treatment) %>% 
      count()
    result <- cohortMethodData$covariates %>% inner_join(select(cohorts, 
                                                                .data$rowId, .data$treatment), by = "rowId") %>% 
      group_by(.data$covariateId, .data$treatment) %>% 
      summarise(sum = sum(.data$covariateValue, na.rm = TRUE), 
                sumSqr = sum(.data$covariateValue^2, na.rm = TRUE)) %>% 
      inner_join(cohortCounts, by = "treatment") %>% mutate(sd = sqrt((.data$sumSqr - 
                                                                         (.data$sum^2/.data$n))/.data$n), mean = .data$sum/.data$n) %>% 
      ungroup() %>% select(.data$covariateId, .data$treatment, 
                           .data$sum, .data$mean, .data$sd) %>% collect()
  }
  target <- result %>% filter(.data$treatment == 1) %>% select(.data$covariateId, 
                                                               sumTarget = .data$sum, meanTarget = .data$mean, sdTarget = .data$sd)
  comparator <- result %>% filter(.data$treatment == 0) %>% 
    select(.data$covariateId, sumComparator = .data$sum, 
           meanComparator = .data$mean, sdComparator = .data$sd)
  result <- target %>% full_join(comparator, by = "covariateId") %>% 
    mutate(sd = sqrt((.data$sdTarget^2 + .data$sdComparator^2)/2)) %>% 
    select(!c(.data$sdTarget, .data$sdComparator))
  return(result)
}


computeCovariateBalanceRe <- function (population, studyPopulation, cohortMethodData, subgroupCovariateId = NULL) {
  
  cohortMethodData$tempCohortsBeforeMatching <- studyPopulation %>%
    select(.data$rowId, .data$treatment)
  cohortMethodData$tempCohortsAfterMatching <- population %>%
    select(.data$rowId, .data$treatment, .data$stratumId)
  
  on.exit(cohortMethodData$tempCohortsBeforeMatching <- NULL)
  on.exit(cohortMethodData$tempCohortsAfterMatching <- NULL,
          add = TRUE)
  
  beforeMatching <- computeMeansPerGroup(cohortMethodData$tempCohortsBeforeMatching,
                                         cohortMethodData)
  afterMatching <- computeMeansPerGroup(cohortMethodData$tempCohortsAfterMatching,
                                        cohortMethodData)
  
  colnames(beforeMatching)[colnames(beforeMatching) == "meanTarget"] <- "beforeMatchingMeanTarget"
  colnames(beforeMatching)[colnames(beforeMatching) == "meanComparator"] <- "beforeMatchingMeanComparator"
  colnames(beforeMatching)[colnames(beforeMatching) == "sumTarget"] <- "beforeMatchingSumTarget"
  colnames(beforeMatching)[colnames(beforeMatching) == "sumComparator"] <- "beforeMatchingSumComparator"
  colnames(beforeMatching)[colnames(beforeMatching) == "sd"] <- "beforeMatchingSd"
  colnames(afterMatching)[colnames(afterMatching) == "meanTarget"] <- "afterMatchingMeanTarget"
  colnames(afterMatching)[colnames(afterMatching) == "meanComparator"] <- "afterMatchingMeanComparator"
  colnames(afterMatching)[colnames(afterMatching) == "sumTarget"] <- "afterMatchingSumTarget"
  colnames(afterMatching)[colnames(afterMatching) == "sumComparator"] <- "afterMatchingSumComparator"
  colnames(afterMatching)[colnames(afterMatching) == "sd"] <- "afterMatchingSd"
  
  balance <- beforeMatching %>% full_join(afterMatching, by = "covariateId") %>%
    inner_join(collect(cohortMethodData$covariateRef), by = "covariateId") %>%
    mutate(beforeMatchingStdDiff = (.data$beforeMatchingMeanTarget - .data$beforeMatchingMeanComparator) / .data$beforeMatchingSd,
           afterMatchingStdDiff = (.data$afterMatchingMeanTarget - .data$afterMatchingMeanComparator) / .data$afterMatchingSd
    )
  
  balance$beforeMatchingStdDiff[balance$beforeMatchingSd == 0] <- 0
  balance$afterMatchingStdDiff[balance$beforeMatchingSd == 0] <- 0
  balance <- balance[order(-abs(balance$beforeMatchingStdDiff)),]
  
  return(balance)
}


computeCovariateBalanceRevise <- function(row, cmOutputFolder, balanceFolder) {
  outputFileName <- file.path(balanceFolder,
                              sprintf("bal_t%s_c%s_o%s_a%s.rds", row$targetId, row$comparatorId, row$outcomeId, row$analysisId))
  if (!file.exists(outputFileName)) {
    ParallelLogger::logTrace("Creating covariate balance file ", outputFileName)
    cohortMethodDataFile <- file.path(cmOutputFolder, row$cohortMethodDataFile)
    cohortMethodData <- CohortMethod::loadCohortMethodData(cohortMethodDataFile)
    strataFile <- file.path(cmOutputFolder, row$strataFile)
    strata <- readRDS(strataFile)
    studyPopFile <- file.path(cmOutputFolder, row$studyPopFile)
    studyPopulation <- readRDS(studyPopFile)
    balance <- computeCovariateBalanceRe(population = strata, studyPopulation = studyPopulation, cohortMethodData = cohortMethodData)
    saveRDS(balance, outputFileName)
  }
}

enforceMinCellValue <- function(data, fieldName, minValues, silent = FALSE) {
  toCensor <- !is.na(pull(data, fieldName)) & pull(data, fieldName) < minValues & pull(data, fieldName) != 0
  if (!silent) {
    percent <- round(100 * sum(toCensor)/nrow(data), 1)
    ParallelLogger::logInfo("   censoring ",
                            sum(toCensor),
                            " values (",
                            percent,
                            "%) from ",
                            fieldName,
                            " because value below minimum")
  }
  if (length(minValues) == 1) {
    data[toCensor, fieldName] <- -minValues
  } else {
    data[toCensor, fieldName] <- -minValues[toCensor]
  }
  return(data)
}