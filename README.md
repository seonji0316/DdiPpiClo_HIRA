DdiPpiClo
==============================


Requirements
============

- A database in [Common Data Model version 5](https://ohdsi.github.io/CommonDataModel/) in one of these platforms: SQL Server, Oracle, PostgreSQL, IBM Netezza, Apache Impala, Amazon RedShift, Google BigQuery, Spark, or Microsoft APS.
- R version 4.0.0 or newer
- On Windows: [RTools](http://cran.r-project.org/bin/windows/Rtools/)
- [Java](http://java.com)
- 25 GB of free disk space

How to run
==========
1. Open your study package in RStudio. Use the following code to install all the dependencies:

    ```r
    renv::deactivate()
    ```

2. In RStudio, select 'Build' then 'Install and Restart' to build the package.

   Once installed, you can execute the study by modifying and using the code below. For your convenience, this code is also provided under `extras/CodeToRun.R`:

    ```r
    library(DdiPpiClo)

    # Optional: specify where the temporary files (used by the Andromeda package) will be created:
    options(andromedaTempFolder = "D:/andromedaTemp")
	
    # Maximum number of cores to be used:
    maxCores <- parallel::detectCores()
	
    # Minimum cell count when exporting data:
    minCellCount <- 5
	
    # The folder where the study intermediate and result files will be written:
    outputFolder <- "D:/DdiPpiClo"
	
    # Details for connecting to the server:
    # See ?DatabaseConnector::createConnectionDetails for help
    connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "sql server",
                                                                server = Sys.getenv("server"),
                                                                user = Sys.getenv("user"),
                                                                password = Sys.getenv("password"),
                                                                pathToDriver = "D:/pathToDriver")
    # The name of the database schema where the CDM data can be found:
    cdmDatabaseSchema <- "cdm_db"

    # The name of the database schema and table where the study-specific cohorts will be instantiated:
    cohortDatabaseSchema <- "scratch.dbo"
    cohortTable <- "ddippicloCohortStudy"

    # Some meta-information that will be used by the export function:
    databaseId <- "DdiPpiClo"
    databaseName <- "DdiPpiClo"
    databaseDescription <- "Drug-drug interaction of PPI and clopidogrel"

    # For some database platforms (e.g. Oracle): define a schema that can be used to emulate temp tables:
    options(sqlRenderTempEmulationSchema = NULL)

    execute(connectionDetails = connectionDetails,
            cdmDatabaseSchema = cdmDatabaseSchema,
            cohortDatabaseSchema = cohortDatabaseSchema,
            cohortTable = cohortTable,
            outputFolder = outputFolder,
            databaseId = databaseId,
            databaseName = databaseName,
            databaseDescription = databaseDescription,
            verifyDependencies = FALSE,
            createCohorts = TRUE,
            synthesizePositiveControls = TRUE,
            runAnalyses = TRUE,
            packageResults = TRUE,
            maxCores = maxCores)
    ```

3. Share the file export/Results_<DatabaseId>.zip in the output folder to the study coordinator

License
=======
The DdiPpiClo package is licensed under Apache License 2.0

Development
===========
DdiPpiClo was developed in ATLAS and R Studio.

### Development status

Unknown
