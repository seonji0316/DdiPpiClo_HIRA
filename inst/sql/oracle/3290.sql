INSERT INTO @target_database_schema.@target_cohort_table (
    cohort_definition_id,
    subject_id,
    cohort_start_date,
    cohort_end_date
)
WITH aspirin_concept as (select distinct c.concept_id
  from @cdm_database_schema.CONCEPT c
  join @cdm_database_schema.CONCEPT_ANCESTOR ca on c.concept_id = ca.descendant_concept_id
  and ca.ancestor_concept_id in (1112807)
  and c.invalid_reason is null
),
target_cohort as (select * 
 FROM @target_database_schema.@target_cohort_table a 
 where a.cohort_definition_id IN (290)
),
aspirin_cohort as
(select * from @cdm_database_schema.DRUG_EXPOSURE b
WHERE b.drug_concept_id IN (select concept_id from aspirin_concept)
AND b.person_id IN (select distinct subject_id FROM target_cohort)
)
SELECT 3290 AS cohort_definition_id, subject_id, cohort_start_date, cohort_end_date
FROM target_cohort a
WHERE EXISTS (
  SELECT 1 
  FROM aspirin_cohort b
  WHERE a.subject_id = b.person_id
  AND b.drug_exposure_start_date <= a.cohort_start_date 
  AND b.drug_exposure_start_date >= DATEADD(day, -30, a.cohort_start_date)
)