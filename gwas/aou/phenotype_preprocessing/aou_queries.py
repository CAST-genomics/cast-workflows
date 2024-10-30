"""
SQL queries for AoU phenotypes
"""

import os

########## Phenotype-specific info ###########
BLOOD_UNITS = ["IU/L", "No matching concept", "international unit per liter", \
            "no value", "unit per liter"]

############################################################
# Demographics query - general for all phenotypes
demographics_sql = """
    SELECT
        person.person_id,
        person.birth_datetime as date_of_birth,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth 
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.person` person   
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` p_sex_at_birth_concept
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
    WHERE
        person.PERSON_ID IN (
            SELECT
                distinct person_id  
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
            WHERE
                cb_search_person.person_id IN (
                    SELECT
                        person_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
                    WHERE
                        has_whole_genome_variant = 1 
                ) 
            )"""

def GetUnits(units):
    if units == "blood":
        return BLOOD_UNITS
    else: return [item.strip() for item in units.split(",")]

def GetPhenotypeRange(range):
    minval, maxval = range.strip().split(",")
    return float(minval), float(maxval)

def ConstructDrugExposureSQL(concept_id):
    return """
    SELECT
        d_exposure.person_id,
        d_exposure.drug_concept_id,
        d_standard_concept.concept_name as standard_concept_name,
        d_standard_concept.concept_code as standard_concept_code,
        d_standard_concept.vocabulary_id as standard_vocabulary,
        d_exposure.drug_exposure_start_datetime,
        d_exposure.drug_exposure_end_datetime,
        d_exposure.verbatim_end_date,
        d_exposure.drug_type_concept_id,
        d_type.concept_name as drug_type_concept_name,
        d_exposure.stop_reason,
        d_exposure.refills,
        d_exposure.quantity,
        d_exposure.days_supply,
        d_exposure.sig,
        d_exposure.route_concept_id,
        d_route.concept_name as route_concept_name,
        d_exposure.lot_number,
        d_exposure.visit_occurrence_id,
        d_visit.concept_name as visit_occurrence_concept_name,
        d_exposure.drug_source_value,
        d_exposure.drug_source_concept_id,
        d_source_concept.concept_name as source_concept_name,
        d_source_concept.concept_code as source_concept_code,
        d_source_concept.vocabulary_id as source_vocabulary,
        d_exposure.route_source_value,
        d_exposure.dose_unit_source_value 
    FROM
        ( SELECT
            * 
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.drug_exposure` d_exposure 
        WHERE
            (
                drug_concept_id IN  (
                    SELECT
                        DISTINCT ca.descendant_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria_ancestor` ca 
                    JOIN
                        (
                            select
                                distinct c.concept_id 
                            FROM
                                `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                            JOIN
                                (
                                    select
                                        cast(cr.id as string) as id 
                                    FROM
                                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr 
                                    WHERE
                                        concept_id IN (""" + str(concept_id) + """
                                        ) 
                                        AND full_text LIKE '%_rank1]%'
                                ) a 
                                    ON (
                                        c.path LIKE CONCAT('%.',
                                    a.id,
                                    '.%') 
                                    OR c.path LIKE CONCAT('%.',
                                    a.id) 
                                    OR c.path LIKE CONCAT(a.id,
                                    '.%') 
                                    OR c.path = a.id) 
                                WHERE
                                    is_standard = 1 
                                    AND is_selectable = 1
                                ) b 
                                    ON (
                                        ca.ancestor_id = b.concept_id
                                    )
                            )
                        )  
                        AND (
                            d_exposure.PERSON_ID IN (
                                SELECT
                                    distinct person_id  
                            FROM
                                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
                            WHERE
                                cb_search_person.person_id IN (
                                    SELECT
                                        person_id 
                                    FROM
                                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
                                    WHERE
                                        has_whole_genome_variant = 1 
                                ) 
                            )
                    )
                ) d_exposure 
            LEFT JOIN
                `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_standard_concept 
                    ON d_exposure.drug_concept_id = d_standard_concept.concept_id 
            LEFT JOIN
                `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_type 
                    ON d_exposure.drug_type_concept_id = d_type.concept_id 
            LEFT JOIN
                `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_route 
                    ON d_exposure.route_concept_id = d_route.concept_id 
            LEFT JOIN
                `""" + os.environ["WORKSPACE_CDR"] + """.visit_occurrence` v 
                    ON d_exposure.visit_occurrence_id = v.visit_occurrence_id 
            LEFT JOIN
                `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_visit 
                    ON v.visit_concept_id = d_visit.concept_id 
            LEFT JOIN
                `""" + os.environ["WORKSPACE_CDR"] + """.concept` d_source_concept 
                    ON d_exposure.drug_source_concept_id = d_source_concept.concept_id"""

def ConstructTraitSQL(concept_id, ppi):
    # PPI and LOINC queries are different in two places.
    # Here we set the values based on PPI or LOINC (default).
    if ppi:
        measurement_concept_id = "measurement_source_concept_id"
        is_standard = "0"
    else:
        measurement_concept_id = "measurement_concept_id"
        is_standard = "1"
    return """
    SELECT
        measurement.person_id,
        measurement.measurement_concept_id,
        m_standard_concept.concept_name as standard_concept_name,
        m_standard_concept.concept_code as standard_concept_code,
        m_standard_concept.vocabulary_id as standard_vocabulary,
        measurement.measurement_datetime,
        measurement.measurement_type_concept_id,
        m_type.concept_name as measurement_type_concept_name,
        measurement.operator_concept_id,
        m_operator.concept_name as operator_concept_name,
        measurement.value_as_number,
        measurement.value_as_concept_id,
        m_value.concept_name as value_as_concept_name,
        measurement.unit_concept_id,
        m_unit.concept_name as unit_concept_name,
        measurement.range_low,
        measurement.range_high,
        measurement.visit_occurrence_id,
        m_visit.concept_name as visit_occurrence_concept_name,
        measurement.measurement_source_value,
        measurement.measurement_source_concept_id,
        m_source_concept.concept_name as source_concept_name,
        m_source_concept.concept_code as source_concept_code,
        m_source_concept.vocabulary_id as source_vocabulary,
        measurement.unit_source_value,
        measurement.value_source_value 
    FROM
        ( SELECT
            * 
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.measurement` measurement 
        WHERE
            (
                """ + measurement_concept_id + """ IN  (
                    SELECT
                        DISTINCT c.concept_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                    JOIN
                        (
                            select
                                cast(cr.id as string) as id 
                            FROM
                                `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr 
                            WHERE
                                concept_id IN (""" + str(concept_id) + """
                                ) 
                                AND full_text LIKE '%_rank1]%'
                        ) a 
                            ON (
                                c.path LIKE CONCAT('%.',
                            a.id,
                            '.%') 
                            OR c.path LIKE CONCAT('%.',
                            a.id) 
                            OR c.path LIKE CONCAT(a.id,
                            '.%') 
                            OR c.path = a.id) 
                        WHERE
                            is_standard = """ + is_standard + """
                            AND is_selectable = 1
                        )
                )  
                AND (
                    measurement.PERSON_ID IN (
                        SELECT
                            distinct person_id  
                        FROM
                            `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
                        WHERE
                            cb_search_person.person_id IN (
                                SELECT
                                    person_id 
                                FROM
                                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
                                WHERE
                                    has_whole_genome_variant = 1 
                            ) 
                        )
                )
            ) measurement 
        LEFT JOIN
            `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_standard_concept 
                ON measurement.measurement_concept_id = m_standard_concept.concept_id 
        LEFT JOIN
            `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_type 
                ON measurement.measurement_type_concept_id = m_type.concept_id 
        LEFT JOIN
            `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_operator 
                ON measurement.operator_concept_id = m_operator.concept_id 
        LEFT JOIN
            `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_value 
                ON measurement.value_as_concept_id = m_value.concept_id 
        LEFT JOIN
            `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_unit 
                ON measurement.unit_concept_id = m_unit.concept_id 
        LEFT JOIn
            `""" + os.environ["WORKSPACE_CDR"] + """.visit_occurrence` v 
                ON measurement.visit_occurrence_id = v.visit_occurrence_id 
        LEFT JOIN
            `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_visit 
                ON v.visit_concept_id = m_visit.concept_id 
        LEFT JOIN
            `""" + os.environ["WORKSPACE_CDR"] + """.concept` m_source_concept 
                ON measurement.measurement_source_concept_id = m_source_concept.concept_id"""

def ConstructSnomedSQL(concept_id):
    return """
    SELECT
        c_occurrence.person_id,
        c_occurrence.condition_concept_id,
        c_standard_concept.concept_name as standard_concept_name,
        c_standard_concept.concept_code as standard_concept_code,
        c_standard_concept.vocabulary_id as standard_vocabulary,
        c_occurrence.condition_start_datetime,
        c_occurrence.condition_end_datetime,
        c_occurrence.condition_type_concept_id,
        c_type.concept_name as condition_type_concept_name,
        c_occurrence.stop_reason,
        c_occurrence.visit_occurrence_id,
        visit.concept_name as visit_occurrence_concept_name,
        c_occurrence.condition_source_value,
        c_occurrence.condition_source_concept_id,
        c_source_concept.concept_name as source_concept_name,
        c_source_concept.concept_code as source_concept_code,
        c_source_concept.vocabulary_id as source_vocabulary,
        c_occurrence.condition_status_source_value,
        c_occurrence.condition_status_concept_id,
        c_status.concept_name as condition_status_concept_name 
    FROM
        ( SELECT
            * 
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.condition_occurrence` c_occurrence 
        WHERE
            (
                condition_concept_id IN (SELECT
                    DISTINCT c.concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id       
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr       
                    WHERE
                        concept_id IN (""" + concept_id + """)       
                        AND full_text LIKE '%_rank1]%'      ) a 
                        ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                        OR c.path LIKE CONCAT('%.', a.id) 
                        OR c.path LIKE CONCAT(a.id, '.%') 
                        OR c.path = a.id) 
                WHERE
                    is_standard = 1 
                    AND is_selectable = 1)
            )) c_occurrence 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_standard_concept 
            ON c_occurrence.condition_concept_id = c_standard_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_type 
            ON c_occurrence.condition_type_concept_id = c_type.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.visit_occurrence` v 
            ON c_occurrence.visit_occurrence_id = v.visit_occurrence_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` visit 
            ON v.visit_concept_id = visit.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_source_concept 
            ON c_occurrence.condition_source_concept_id = c_source_concept.concept_id 
    LEFT JOIN
        `""" + os.environ["WORKSPACE_CDR"] + """.concept` c_status 
            ON c_occurrence.condition_status_concept_id = c_status.concept_id"""

