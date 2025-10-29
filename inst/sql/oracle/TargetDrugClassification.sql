With omeprazole as (select distinct c.concept_id
  from @cdm_database_schema.CONCEPT c
  join @cdm_database_schema.CONCEPT_ANCESTOR ca on c.concept_id = ca.descendant_concept_id
  and ca.ancestor_concept_id in (21600096, 923645)
  and c.invalid_reason is null
),
esomeprazole as (select distinct c.concept_id
  from @cdm_database_schema.CONCEPT c
  join @cdm_database_schema.CONCEPT_ANCESTOR ca on c.concept_id = ca.descendant_concept_id
  and ca.ancestor_concept_id in (21600100, 904453)
  and c.invalid_reason is null
),
drug_record as (select * 
	from @cdm_database_schema.DRUG_ERA
	where person_id IN (@person_id)
    AND drug_concept_id IN (
      SELECT concept_id FROM omeprazole
      UNION
      SELECT concept_id FROM esomeprazole
    )
)
SELECT 
    dr.person_id,
    dr.drug_era_id,
    dr.drug_concept_id,
    CASE 
        WHEN dr.drug_concept_id IN (SELECT concept_id FROM omeprazole) THEN 'omeprazole'
        WHEN dr.drug_concept_id IN (SELECT concept_id FROM esomeprazole) THEN 'esomeprazole'
        ELSE 'other'
    END AS drug_group,
    dr.drug_era_start_date,
    dr.drug_era_end_date
FROM drug_record dr;