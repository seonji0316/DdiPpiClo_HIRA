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
drug_record as (
    select de.*, 
           combined_drugs.drug_group
    from @cdm_database_schema.DRUG_ERA de
    join (
        select concept_id, 'omeprazole' as drug_group 
        from omeprazole
        union all
        select concept_id, 'esomeprazole' as drug_group 
        from esomeprazole
    ) combined_drugs 
    on de.drug_concept_id = combined_drugs.concept_id
),
target_record as (select * 
	from @target_database_schema.@target_cohort_table
	where cohort_definition_id = @target_id
)
SELECT 
    dr.person_id,
    dr.drug_era_id,
    dr.drug_concept_id,
    dr.drug_era_start_date,
    dr.drug_era_end_date,
    dr.drug_group
FROM target_record tr
JOIN drug_record dr
    ON tr.subject_id = dr.person_id
    AND tr.cohort_start_date = dr.drug_era_start_date