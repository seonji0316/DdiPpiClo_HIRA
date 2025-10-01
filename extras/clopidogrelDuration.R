clopidogrelDuration <- function(target_id, 
                                comparator_id, 
                                outcome_id, 
                                analysis_id, 
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
  
  stratPop <- readRDS(file.path(cmOutput,omIr$strataFile))
  stratPop <- inner_join(stratPop, cohort[,c("rowId","personId")], by = "rowId")
  
  sqlFolder <- file.path(getwd(), "inst/sql/sql_server")
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
    mutate(diffdays = as.numeric(difftime(stratPop$clopidogrelEndDate, stratPop$cohortStartDate, units = "days")),
           cloDuration = ifelse(diffdays >= timeAtRisk, timeAtRisk, diffdays),
           outcomeDuringExposure = ifelse(cloDuration >= daysToEvent, 1, 0),
           gapToEvent = daysToEvent-cloDuration)
  
  cloDurationResult <- stratPop %>%
    group_by(cohortDefinitionId) %>%
    summarise(numOfPatients = n(),
              min = min(cloDuration, na.rm = T),
              max = max(cloDuration, na.rm = T),
              q1 = quantile(cloDuration, 0.25, na.rm = TRUE),
              median =  quantile(cloDuration, 0.5, na.rm = TRUE),
              q3 = quantile(cloDuration, 0.75, na.rm = TRUE),
              mean = mean(cloDuration, na.rm = T),
              sd = sd(cloDuration, na.rm = T)) %>%
    mutate(outcomeDuringExposure = NA,
           analysisType = "duration")
  
  outcomeResult <- stratPop %>%
    filter(outcomeCount >= 1) %>%
    group_by(cohortDefinitionId, outcomeDuringExposure) %>%
    summarise(numOfPatients = n(),
              min = min(gapToEvent, na.rm = T),
              max = max(gapToEvent, na.rm = T),
              q1 = quantile(gapToEvent, 0.25, na.rm = TRUE),
              median =  quantile(gapToEvent, 0.5, na.rm = TRUE),
              q3 = quantile(gapToEvent, 0.75, na.rm = TRUE),
              mean = mean(gapToEvent, na.rm = T),
              sd = sd(gapToEvent, na.rm = T)) %>%
    mutate(analysisType = "outcome")
  
  final <- bind_rows(cloDurationResult, outcomeResult) %>% 
    mutate(targetId = target_id, 
           comparatorId = comparator_id, 
           outcomeId = outcome_id,
           analysisId = analysis_id
           ) %>%
    select(targetId, comparatorId, outcomeId, analysisId, analysisType, cohortDefinitionId, outcomeDuringExposure, numOfPatients, min, max, q1, median, q3, mean, sd)
  
  
  
  file_name <- sprintf("clopidogrel_duration_t%s_c%s_o%s_a%s.csv", target_id, comparator_id, outcome_id, analysis_id)
  write.csv(final, file.path(outputFolder, "export", file_name), row.names = F)
  
  DatabaseConnector::disconnect(connection)
  print("Done!")
}