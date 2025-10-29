INSERT INTO @target_database_schema.@target_cohort_table (
    cohort_definition_id,
    subject_id,
    cohort_start_date,
    cohort_end_date
)
WITH diabetes_concept as (select distinct c.concept_id
  from @cdm_database_schema.CONCEPT c
  join @cdm_database_schema.CONCEPT_ANCESTOR ca on c.concept_id = ca.descendant_concept_id
  and ca.ancestor_concept_id in (201820)
  and c.invalid_reason is null
),
target_cohort as (select * 
 FROM @target_database_schema.@target_cohort_table a 
 where a.cohort_definition_id IN (289)
),
diabetes_cohort as
(select * from @cdm_database_schema.CONDITION_OCCURRENCE b
WHERE b.condition_concept_id IN (select concept_id from diabetes_concept)
AND b.person_id IN (select distinct subject_id FROM target_cohort)
)
SELECT 1289 AS cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
FROM target_cohort a
WHERE EXISTS (
  SELECT 1
  FROM diabetes_cohort b
  WHERE a.subject_id = b.person_id
  AND b.condition_start_date <= a.cohort_start_date 
  AND b.condition_start_date >= DATEADD(day, -365, a.cohort_start_date)
)