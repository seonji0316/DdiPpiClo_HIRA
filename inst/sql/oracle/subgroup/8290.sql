INSERT INTO @target_database_schema.@target_cohort_table (
    cohort_definition_id,
    subject_id,
    cohort_start_date,
    cohort_end_date
)
SELECT 8290 AS cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
FROM @target_database_schema.@target_cohort_table a
WHERE a.cohort_definition_id = 290
AND EXISTS (
	SELECT 1
	FROM @cdm_database_schema.PERSON p 
	WHERE a.subject_id = p.person_id
	AND  p.gender_concept_id = 8532
)