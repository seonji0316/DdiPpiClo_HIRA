analysisResult <- function(target_id,
                           comparator_id,
                           outcome_id,
                           analysis_id,
                           database_id,
                           outputFolder,
                           connectionDetails,
                           cdmDatabaseSchema,
                           cohortDatabaseSchema,
                           cohortTable) {
  
  connection <- DatabaseConnector::connect(connectionDetails) 
  cmOutput <- file.path(outputFolder, "cmOutput")
  om <- readRDS(file.path(cmOutput, "outcomeModelReference.rds")) 
  omIr <- om %>% filter(targetId == target_id,
                        comparatorId == comparator_id,
                        outcomeId == outcome_id,
                        analysisId == analysis_id)
  
  cmdata <- CohortMethod::loadCohortMethodData(file.path(cmOutput,omIr$cohortMethodDataFile)) 
  cohort <- cmdata$cohorts %>% as.data.frame()
  cohort$personId <- as.numeric(cohort$personId)
  
  stratPop <- readRDS(file.path(cmOutput, omIr$strataFile))
  stratPop <- inner_join(stratPop, cohort [,c("rowId", "personId")], by = "rowId")
  
  if (connectionDetails$dbms == "oracle") {
    sqlFolder <- file.path(getwd(), "inst/sql/oracle")
  } else {
    sqlFolder <- file.path(getwd(), "inst/sql/sql_server")
  }
  
  sql <- SqlRender::readSql(file.path(sqlFolder, "GetClopidogrelRecord.sql"))
  
  clopidogrelRecord <- DatabaseConnector::renderTranslateQuerySql(connection,
                                                                  sql,
                                                                  cdm_database_schema = cdmDatabaseSchema,
                                                                  target_database_schema = cohortDatabaseSchema, 
                                                                  target_cohort_table = cohortTable,
                                                                  cohort_ids = c(target_id, comparator_id), 
                                                                  target_id = target_id,
                                                                  comparator_id = comparator_id)
  colnames(clopidogrelRecord) <- SqlRender::snakeCaseToCamelCase(colnames(clopidogrelRecord))
  
  
  stratPop <- inner_join(stratPop, clopidogrelRecord, by = c("personId", "cohortStartDate", "treatment"))
  
  stratPop <- stratPop %>%
    mutate(cloDays = as.numeric(difftime(stratPop$clopidogrelEndDate, stratPop$cohortStartDate, units = "days"))) %>% 
    filter(cloDays > 0) %>%
    mutate(combDays = pmin(daysToCohortEnd, cloDays, na.rm=T),
           daysToEvent = ifelse(daysToEvent > combDays, NA, daysToEvent), 
           outcomeCount = ifelse(is.na(daysToEvent), 0, 1),
           survivalTime = pmin(combDays, daysToEvent, na.rm = T))
  
  outcomeModel <- CohortMethod::fitOutcomeModel(population = stratPop,
                                                modelType = "cox",
                                                stratified = FALSE)
                                                 
                                                
  coefficient <- as.vector(coef(outcomeModel))
  ci <- confint(outcomeModel)
  
  if(is.null(coefficient)) {
    p <- NA
  } else {
    p <- EmpiricalCalibration::computeTraditionalP(logRr = coefficient,
                                                   seLogRr = outcomeModel$outcomeModelTreatmentEstimate$seLogRr)

  }
  
  result <- data.frame(targetId = target_id,
                       comparatorId = comparator_id,
                       outcomeId = outcome_id,
                       analysisId = analysis_id,
                       rr = if (is.null(coefficient)) NA else exp(coefficient), 
                       ci95Lb = if (is.null(coefficient)) NA else exp(ci[1]),
                       ci95Ub = if (is.null(coefficient)) NA else exp(ci[2]),
                       p = !!p,
                       i_2 = NA,
                       logRr = if (is.null(coefficient)) NA else coefficient,
                       seLogRr = if (is.null(coefficient)) NA else outcomeModel$outcomeModelTreatmentEstimate$seLogRr,
                       targetSubjects = outcomeModel$populationCounts $targetPersons,
                       comparatorSubjects = outcomeModel$populationCounts$comparatorPersons,
                       targetDays = outcomeModel$timeAtRisk$targetDays,
                       comparatorDays = outcomeModel$timeAtRisk$comparatorDays,
                       targetOutcomes = outcomeModel$outcomeCounts$targetOutcomes,
                       comparatorOutcomes = outcomeModel$outcomeCounts$comparatorOutcomes,
                       calibratedP = NA, 
                       calibratedRr = NA,
                       calibratedCi95Lb = NA,
                       calibratedC195ub = NA,
                       calibratedLogRr = NA,
                       calibratedSeLogRr = NA,
                       databaseId = database_id 
                       )
  
  file_name <- sprintf("concomitant_cohort_method_result_t%s_c%s_o%s_a%s.csv", target_id, comparator_id, outcome_id, analysis_id) 
  write.csv(result, file.path(outputFolder, "export", file_name), row.names = F)
  
  print("Done!")
  
}