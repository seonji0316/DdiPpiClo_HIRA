INSERT INTO @target_database_schema.@target_cohort_table (
    cohort_definition_id,
    subject_id,
    cohort_start_date,
    cohort_end_date
)
SELECT 5289 AS cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
FROM @target_database_schema.@target_cohort_table a
WHERE a.cohort_definition_id = 289
AND EXISTS (
	SELECT 1
	FROM @cdm_database_schema.PERSON p 
	WHERE a.subject_id = p.person_id
	AND YEAR(a.cohort_start_date) - p.year_of_birth < 65
)