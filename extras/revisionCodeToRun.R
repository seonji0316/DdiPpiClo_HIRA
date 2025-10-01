library(DdiPpiClo)
library(dplyr)

options(encoding = "UTF-8")

# Optional: specify where the temporary files (used by the Andromeda package) will be created:
options(andromedaTempFolder = "")

# Maximum number of cores to be used:
maxCores <- parallel::detectCores()

# The folder where the study intermediate and result files will be written:
outputFolder <- ""

# Details for connecting to the server:
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "",
                                                                server="",
                                                                user="",
                                                                port = "",
                                                                password="",
                                                                pathToDriver = "")


# The name of the database schema where the CDM data can be found:
cdmDatabaseSchema <- ""

# The name of the database schema and table where the study-specific cohorts will be instantiated:
cohortDatabaseSchema <-  ""
cohortTable <- "" # You must use the cohortTable previously created

# Some meta-information that will be used by the export function:
databaseId <- ""
databaseName <- ""
databaseDescription <- "Drug-drug interaction of PPI and clopidogrel"

# For some database platforms (e.g. Oracle): define a schema that can be used to emulate temp tables:
options(sqlRenderTempEmulationSchema = NULL)

#### Make subgroup cohorts ####

######################### CAUTION: It must be run only once !!!! #################################

connection <- DatabaseConnector::connect(connectionDetails) 

sqlFolder <- file.path(getwd(), "inst/sql/sql_server")
sql <- SqlRender::readSql(file.path(sqlFolder, "CombinePPICohort.sql"))

DatabaseConnector::renderTranslateExecuteSql(connection, 
                                             sql, 
                                             cdm_database_schema = cdmDatabaseSchema,
                                             target_database_schema = cohortDatabaseSchema,
                                             target_cohort_table = cohortTable)

fileList <- list.files(file.path(getwd(), "inst/sql/sql_server/subgroup"))

for (i in 1:length(fileList)) {
  
  print(paste0("Subgroup: ", sub("\\.sql$", "", fileList[i])))
  sql <- SqlRender::readSql(file.path(sqlFolder, "subgroup", fileList[i]))
  
  DatabaseConnector::renderTranslateExecuteSql(connection,
                                               sql, 
                                               cdm_database_schema = cdmDatabaseSchema,
                                               target_database_schema = cohortDatabaseSchema,
                                               target_cohort_table = cohortTable)
  if (i==length(fileList)) {
    print("Done!")
  } 
} 

DatabaseConnector::disconnect(connection)

#######################################################################################################

execute(connectionDetails = connectionDetails,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cohortDatabaseSchema = cohortDatabaseSchema,
        cohortTable = cohortTable,
        outputFolder = outputFolder,
        databaseId = databaseId,
        databaseName = databaseName,
        databaseDescription = databaseDescription,
        verifyDependencies = FALSE,
        createCohorts = FALSE,
        synthesizePositiveControls = FALSE,
        runAnalyses = TRUE,
        packageResults = TRUE,
        maxCores = maxCores)

#### Ome vs Esome ####
connection <- DatabaseConnector::connect(connectionDetails)
cmOutput <- file.path(outputFolder, "cmOutput")
om <- readRDS(file.path(cmOutput, "outcomeModelReference.rds"))

omIr <- om %>% filter(analysisId == 2, targetId == 289, comparatorId == 290, outcomeId == 70)
cmdata <- CohortMethod::loadCohortMethodData(file.path(cmOutput,omIr$cohortMethodDataFile))
cohort <- cmdata$cohorts %>% as.data.frame()
cohort$personId <- as.numeric(cohort$personId)

stratPop <- readRDS(file.path(cmOutput,omIr$strataFile))
stratPop <- inner_join(stratPop, cohort[,c("rowId","personId")], by = "rowId")

stratPop <- stratPop %>% filter(treatment == 1)

sqlFolder <- file.path(getwd(), "inst/sql/sql_server")
sql <- SqlRender::readSql(file.path(sqlFolder, "TargetDrugClassification.sql"))

drugRecord <- DatabaseConnector::renderTranslateQuerySql(connection,
                                                       sql, 
                                                       cdm_database_schema = cdmDatabaseSchema,
                                                       person_id = as.numeric(stratPop$personId))

colnames(drugRecord) <- SqlRender::snakeCaseToCamelCase(colnames(drugRecord))

stratPop <- inner_join(stratPop, drugRecord, by = c("personId" = "personId", "cohortStartDate" = "drugEraStartDate"))
drugClassification <- stratPop %>% group_by(drugGroup) %>% summarise(numOfRecord = n())
drugClassification <- drugClassification %>% mutate(analysisId = 2, targetId = 289, comparatorId = 290, outcomeId = 70) %>%
  select(targetId, comparatorId, outcomeId, analysisId, drugGroup, numOfRecord)

file_name <- sprintf("targetDrugClassification_t%i_c%i_o%i_a%i.csv", 289, 290, 70, 2)

write.csv(drugClassification, file.path(outputFolder, "export", file_name), row.names = F)

#### Clopidogrel duration ####
source("./extras/clopidogrelDuration.R")
clopidogrelDuration(target_id = 289, 
                    comparator_id = 290, 
                    outcome_id = 70, 
                    analysis_id = 2, 
                    outputFolder = outputFolder,
                    connectionDetails = connectionDetails,
                    cdmDatabaseSchema = cdmDatabaseSchema,
                    cohortDatabaseSchema = cohortDatabaseSchema,
                    cohortTable = cohortTable) 

clopidogrelDuration(target_id = 289, 
                    comparator_id = 290, 
                    outcome_id = 70, 
                    analysis_id = 1, 
                    outputFolder = outputFolder,
                    connectionDetails = connectionDetails,
                    cdmDatabaseSchema = cdmDatabaseSchema,
                    cohortDatabaseSchema = cohortDatabaseSchema,
                    cohortTable = cohortTable) 