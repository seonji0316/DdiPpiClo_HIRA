WITH target_cohort as (
	SELECT * 
	FROM @target_database_schema.@target_cohort_table
	WHERE cohort_definition_id IN (@cohort_ids)
),
clopidogrel_concept as (select distinct c.concept_id
  from @cdm_database_schema.CONCEPT c
  join @cdm_database_schema.CONCEPT_ANCESTOR ca on c.concept_id = ca.descendant_concept_id
  and ca.ancestor_concept_id in (21600989, 1322184)
  and c.invalid_reason is null
),
drug_record as (
    select de.*
    from @cdm_database_schema.DRUG_ERA de
    join clopidogrel_concept c 
    on de.drug_concept_id = c.concept_id
)
SELECT
  de.person_id,
  de.drug_era_start_date AS clopidogrel_start_date,
  de.drug_era_end_date   AS clopidogrel_end_date,
  t.cohort_definition_id,
  t.cohort_start_date,
  t.cohort_end_date,
  CASE
    WHEN t.cohort_definition_id = @target_id     THEN 1
    WHEN t.cohort_definition_id = @comparator_id THEN 0
  END AS treatment
FROM drug_record de
JOIN target_cohort t
  ON de.person_id = t.subject_id
WHERE de.drug_era_start_date <= DATEADD(day,-30,t.cohort_start_date)
  AND de.drug_era_end_date   >= t.cohort_start_date