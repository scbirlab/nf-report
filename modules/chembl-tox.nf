process fetch_chembl_tox {

    tag "v${chembl_version}:${chembl_url}:${cell_ids.join(',')}"

    publishDir( 
        "${params.outputs}/toxicity", 
        mode: 'copy',
    )
    
    errorStrategy { if ( "${chembl_db}" == 'placeholder' ) { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' } else { return 'terminate' } }
    maxRetries 5

    input:
    val chembl_url
    val cell_ids
    val chembl_version
    path chembl_db

    output:
    path "chembl_tox.tsv.gz", emit: main

    script:
    if ( "${chembl_db}" == 'placeholder' ) {
    """
    set -x
    # == config 
    SLEEP_TIME=0.3

    parse_assay () (
        jq -r '.assays[] 
        | [
            .cell_chembl_id, 
            .assay_cell_type, 
            .assay_chembl_id, 
            .assay_type
        ] | @tsv'
    )

    parse_activity () (
        jq -r '
            .activities[] 
            | [
                .assay_chembl_id, 
                .molecule_chembl_id,
                .canonical_smiles,
                .standard_type,
                .standard_value,
                .standard_units
            ] | @tsv
        '
    )

    parse_mechansisms () (
        jq -r '
            .mechanisms[] 
            | [
                .molecule_chembl_id, 
                .target_chembl_id, 
                .mechanism_of_action, 
                .max_phase
            ] | @tsv' 
    )

    get_col_number () (
        head -n1 | tr \$'\\t' \$'\\n' | grep -nFw "\$1" | cut -d: -f1
    )

    # == Get all cell line assays
    root_url="${chembl_url}/chembl/api/data/assay.json"

    query="assay_cell_type__in=${cell_ids.join(',')}&assay_type__in=F,T&cell_chembl_id__isnull=False"
    init_url="\${root_url}?\${query}&limit=0"

    header=(cell_chembl_id cell_type assay_chembl_id assay_type)
    
    printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
    > chembl_assays.tsv

    curl -s "\${init_url}" > response.json

    parse_assay < response.json >> chembl_assays.tsv
    jq -r '.page_meta.next' < response.json > next_page.txt

    while [ "\$(cat next_page.txt)" != "null" ]
    do  
        sleep \$SLEEP_TIME
        curl -s  -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
            "${chembl_url}\$(cat next_page.txt)" > response.json
        parse_assay < response.json >> chembl_assays.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt
    done    

    head -n1 chembl_assays.tsv \
    | cat - <(tail -n+2 chembl_assays.tsv | sort -u | sort -k1 ) \
    > chembl_assays-sorted.tsv \
    && mv chembl_assays-sorted.tsv chembl_assays.tsv


    # == Get all IC50 activities
    root_url="${chembl_url}/chembl/api/data/activity.json"
    base_query="standard_type__in=IC50,CC50&potential_duplicate=0&standard_value__gte=0&standard_units=nM"

    header=(assay_chembl_id molecule_chembl_id molecule_smiles assay_measurement_type assay_standard_value assay_units)
    
    printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
    > chembl_ic50.tsv
    
    assay_id_col=\$(get_col_number assay_chembl_id < chembl_assays.tsv)
    tail -n+2 chembl_assays.tsv | cut -f"\$assay_id_col" | sort -u |split -l 20 - 'ids_'
    for id_file in ids_*
    do
        sleep \$SLEEP_TIME
        ids=\$(tr \$'\\n' , < "\$id_file")
        query="assay_chembl_id__in=\${ids}"
        init_url="\${root_url}?\${base_query}&\${query}&limit=0"

        curl -s  -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
            "\${init_url}" > response.json

        parse_activity < response.json >> chembl_ic50.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt

        while [ "\$(cat next_page.txt)" != "null" ]
        do  
            sleep \$SLEEP_TIME
            curl -s "${chembl_url}\$(cat next_page.txt)" > response.json
            parse_activity < response.json >> chembl_ic50.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt
        done    
    done

    # == Get all molecule targets
    root_url="${chembl_url}/chembl/api/data/mechanism.json"
    base_query="direct_interaction=1&molecular_mechanism=1&action_type__in=ANTAGONIST,INHIBITOR"

    header=(molecule_chembl_id target_chembl_id molecule_mechanism moleculae_max_phase)
    
    printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
    > targets.tsv
    
    mol_col=\$(get_col_number molecule_chembl_id < chembl_ic50.tsv)
    tail -n+2 chembl_ic50.tsv | cut -f"\$mol_col" | sort -u | split -l 20 - 'mols_'
    for id_file in mols_*
    do
        sleep \$SLEEP_TIME
        ids=\$(tr \$'\\n' , < "\$id_file")
        query="molecule_chembl_id__in=\${ids}"
        init_url="\${root_url}?\${base_query}&\${query}&limit=0"

        curl -s  -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
            "\${init_url}" > response.json

        parse_mechansisms < response.json >> targets.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt
        while [ "\$(cat next_page.txt)" != "null" ]
        do  
            sleep \$SLEEP_TIME
            curl -s "${chembl_url}\$(cat next_page.txt)" > response.json
            parse_mechansisms < response.json >> targets.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt
        done  
    done

    # == Get all molecule-target human IC50 or Ki
    parse_inhibition () (
        jq -r '.activities[] | [
            .molecule_chembl_id, 
            .target_chembl_id,
            .standard_type,
            .standard_value,
            .standard_units
        ] | @tsv'
    )
    root_url="${chembl_url}/chembl/api/data/activity.json"
    base_query="assay_type=B&standard_type__in=Ki,IC50&pchembl_value__gte=0&potential_duplicate=0&standard_units=nM"

    header=(molecule_chembl_id target_chembl_id molecule_target_measurement molecule_target_inhibition_value molecule_target_inhibition_units)
    
    printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
    > inhibition.tsv
    
    mol_col=\$(get_col_number molecule_chembl_id < targets.tsv)
    target_col=\$(get_col_number target_chembl_id < targets.tsv)
    tail -n+2 targets.tsv | while read target_line
    do
        sleep \$SLEEP_TIME
        mol_id=\$(echo "\$target_line" | cut -f"\$mol_col")
        target_id=\$(echo "\$target_line" | cut -f"\$target_col")
        query="molecule_chembl_id=\${mol_id}&target_chembl_id=\${target_id}"
        init_url="\${root_url}?\${base_query}&\${query}&limit=0"

        curl -s "\${init_url}" > response.json

        parse_inhibition < response.json >> inhibition.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt
        while [ "\$(cat next_page.txt)" != "null" ]
        do  
            sleep \$SLEEP_TIME
            curl -s -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
                "${chembl_url}\$(cat next_page.txt)" > response.json
            parse_inhibition < response.json >> inhibition.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt
        done  
    done

    python -c '
    import pandas as pd

    (
        pd.merge(
            pd.read_csv("chembl_assays.tsv", sep="\\t"),
            pd.read_csv("chembl_ic50.tsv", sep="\\t"),
        )
        .merge(
            pd.read_csv("targets.tsv", sep="\\t", how="left")
        )
        .merge(
            pd.read_csv("inhibition.tsv", sep="\\t", how="left")
        )
        .drop_duplicates()
        .to_csv("chembl_tox.tsv.gz", sep="\\t", index=False)
    )
    
    '

    """
    }
    else {
        """
        set -euo pipefail

        echo "Using local ChEMBL SQLite DB: ${chembl_db}" >&2

        # DuckDB + SQLite extension query to reproduce the same columns as the API
        mkdir duckdb
        duckdb << 'EOF'
        SET home_directory='duckdb';
        INSTALL sqlite;
        LOAD sqlite;

        ATTACH '${chembl_db}' AS chembl (TYPE sqlite, READ_ONLY);
        USE chembl;

        -- 1) All cell line assays matching the requested cell_ids
        CREATE TEMP TABLE chembl_assays AS
        SELECT DISTINCT
            a.cell_id          AS cell_chembl_id,
            a.assay_cell_type  AS cell_type,
            a.chembl_id        AS assay_chembl_id,
            a.assay_type       AS assay_type
        FROM assays a
        LEFT JOIN cell_dictionary cd
            ON a.cell_id = cd.cell_id
        WHERE
            a.assay_cell_type IN (
                ${cell_ids.collect { "'${it}'" }.join(',')}
            )
        AND a.assay_type IN ('F','T','A')
        AND a.cell_id IS NOT NULL;

        -- 2) All IC50/CC50 activities for those assays
        CREATE TEMP TABLE chembl_ic50 AS
        SELECT DISTINCT
            a.chembl_id         AS assay_chembl_id,
            a.assay_type        AS assay_type,
            a.confidence_score  AS assay_target_confidence_score,
            md.chembl_id         AS molecule_chembl_id,
            md.pref_name         AS molecule_name,
            cs.canonical_smiles  AS molecule_smiles,
            act.standard_type    AS assay_measurement_type,
            act.standard_value   AS assay_standard_value,
            act.standard_units   AS assay_units,
            act.pchembl_value   AS molecule_target_pchembl
        FROM activities act
        JOIN assays a               ON act.assay_id = a.assay_id
        JOIN chembl_assays ca       ON a.chembl_id  = ca.assay_chembl_id
        JOIN molecule_dictionary md ON act.molregno = md.molregno
        LEFT JOIN compound_structures cs
            ON act.molregno = cs.molregno
        WHERE
            act.standard_type IN ('IC50','CC50')
            AND COALESCE(act.potential_duplicate, 0) = 0
            AND act.standard_value >= 0
            AND act.standard_units = 'nM';

        -- 3) Mechanisms for molecules in those IC50 assays
        CREATE TEMP TABLE tox_targets AS
        WITH mech_raw AS (
            SELECT DISTINCT
                md.chembl_id             AS molecule_chembl_id,
                td.chembl_id             AS mech_target_chembl_id,
                td.target_type,
                td.tid                   AS mech_tid,
                mech.mechanism_of_action AS molecule_mechanism,
                di.max_phase_for_ind     AS molecule_max_phase
            FROM drug_mechanism mech
            JOIN molecule_dictionary md  ON mech.molregno = md.molregno
            LEFT JOIN drug_indication di      ON mech.molregno = di.molregno
            JOIN target_dictionary td    ON mech.tid = td.tid
            JOIN chembl_ic50 ci          ON md.chembl_id = ci.molecule_chembl_id
            WHERE mech.direct_interaction = 1
            AND mech.molecular_mechanism = 1
        )
        SELECT DISTINCT
            molecule_chembl_id,
            mech_target_chembl_id AS target_chembl_id,
            molecule_mechanism,
            molecule_max_phase
        FROM mech_raw
        WHERE target_type = 'SINGLE PROTEIN'

        UNION ALL

        SELECT DISTINCT
            mr.molecule_chembl_id,
            t_child.chembl_id     AS target_chembl_id,
            mr.molecule_mechanism,
            mr.molecule_max_phase
        FROM mech_raw mr
        JOIN target_relations tr     ON mr.mech_tid = tr.tid
        JOIN target_dictionary t_child
            ON tr.related_tid = t_child.tid
            AND t_child.target_type = 'SINGLE PROTEIN'
        WHERE mr.target_type IN ('PROTEIN COMPLEX', 'PROTEIN FAMILY',
                                'PROTEIN COMPLEX GROUP')
        AND tr.relationship = 'SUPERSET OF';

        -- 4) Ki / IC50 inhibition data for those molecule–target pairs
        CREATE TEMP TABLE inhibition AS
        SELECT DISTINCT
            md.chembl_id       AS molecule_chembl_id,
            t.chembl_id        AS target_chembl_id,
            act.standard_type  AS molecule_target_measurement,
            act.standard_value AS molecule_target_inhibition_value,
            act.standard_units AS molecule_target_inhibition_units
        FROM activities act
        JOIN molecule_dictionary md 
            ON act.molregno = md.molregno
        JOIN assays a            
            ON act.assay_id = a.assay_id
        JOIN target_dictionary t 
            ON a.tid = t.tid
        JOIN tox_targets tg
            ON md.chembl_id = tg.molecule_chembl_id
            AND t.chembl_id = tg.target_chembl_id
        WHERE
            a.assay_type = 'B'
            AND act.standard_type IN ('Ki','IC50')
            AND act.pchembl_value >= 0
            AND COALESCE(act.potential_duplicate, 0) = 0
            AND act.standard_units = 'nM';

        COPY (
            SELECT DISTINCT *
            FROM chembl_assays
            JOIN chembl_ic50 USING (assay_chembl_id)
            LEFT JOIN tox_targets USING (molecule_chembl_id)
            LEFT JOIN inhibition  USING (molecule_chembl_id, target_chembl_id)
        ) TO 'chembl_tox.tsv' (HEADER, DELIMITER '\\t');
        EOF

        gzip --best chembl_tox.tsv

        """
    }

}