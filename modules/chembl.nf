process chembl_status {

    tag "${date}"

    errorStrategy "retry"  // sometimes doesn't respond
    maxRetries 4

    publishDir( 
        "${params.outputs}/info", 
        mode: 'copy'
    )

    input:
    val date

    output:
    env 'chembl_version', emit: version
    path 'chembl_version.json', emit: info

    script:
    """
    set +x
    curl -v \
        -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
        -H 'Accept: application/json' \
        https://www.ebi.ac.uk/chembl/api/data/status.json \
    > chembl_version.json
    chembl_version=\$(jq -r '.chembl_db_version' < chembl_version.json)
    echo \$chembl_version

    """

}


process fetch_chembl_target_sequences {

    tag "v${chembl_version}"
    stageInMode 'link'

    publishDir( 
        "${params.outputs}/sequences", 
        mode: 'copy',
        saveAs: { "chembl_targets-v${chembl_version}.fasta.gz" }
    )

    input:
    val chembl_url
    val chembl_version

    output:
    tuple val( chembl_version ), path( "chembl_targets.fasta.gz" )

    script:
    """
    curl -s "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_${chembl_version}/chembl_${chembl_version}.fa.gz" \
    | zcat \
    | sed 's/^> />/' \
    | gzip --best \
    > chembl_targets.fasta.gz
    
    """

}


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
    base_query="target_tax_id=9606&assay_type=B&standard_type__in=Ki,IC50&pchembl_value__gte=0&potential_duplicate=0&standard_units=nM"

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
            pd.read_csv("targets.tsv", sep="\\t")
        )
        .merge(
            pd.read_csv("inhibition.tsv", sep="\\t")
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
        AND a.assay_type IN ('F','T')
        AND a.cell_id IS NOT NULL;

        -- 2) All IC50/CC50 activities for those assays
        CREATE TEMP TABLE chembl_ic50 AS
        SELECT DISTINCT
            a.chembl_id          AS assay_chembl_id,
            md.chembl_id         AS molecule_chembl_id,
            cs.canonical_smiles  AS molecule_smiles,
            act.standard_type    AS assay_measurement_type,
            act.standard_value   AS assay_standard_value,
            act.standard_units   AS assay_units
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
        SELECT DISTINCT
            md.chembl_id             AS molecule_chembl_id,
            td.chembl_id             AS target_chembl_id,
            mech.mechanism_of_action AS molecule_mechanism,
            di.max_phase_for_ind     AS molecule_max_phase
        FROM drug_mechanism mech
        JOIN molecule_dictionary md 
            ON mech.molregno = md.molregno
        JOIN drug_indication di 
            ON mech.molregno = di.molregno
        JOIN target_dictionary td
            ON mech.tid = td.tid
        JOIN chembl_ic50 ci
            ON md.chembl_id = ci.molecule_chembl_id
        WHERE
            mech.direct_interaction = 1
            AND mech.molecular_mechanism = 1
            AND mech.action_type IN ('ANTAGONIST','INHIBITOR');

        -- 4) Human (tax_id=9606) Ki / IC50 inhibition data for those molecule–target pairs
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
            t.tax_id = 9606
            AND a.assay_type = 'B'
            AND act.standard_type IN ('Ki','IC50')
            AND act.pchembl_value >= 0
            AND COALESCE(act.potential_duplicate, 0) = 0
            AND act.standard_units = 'nM';

        COPY (
            SELECT DISTINCT *
            FROM chembl_assays
            JOIN chembl_ic50 USING (assay_chembl_id)
            JOIN tox_targets USING (molecule_chembl_id)
            JOIN inhibition  USING (molecule_chembl_id, target_chembl_id)
        ) TO 'chembl_tox.tsv' (HEADER, DELIMITER '\\t');
        EOF

        gzip --best chembl_tox.tsv

        """
    }

}


process fetch_chembl_targets {

    tag "v${chembl_version}:${chembl_url}"

    publishDir( 
        "${params.outputs}/targets", 
        mode: 'copy',
        saveAs: { "${chembl_version}.${it}" },
    )
    
    input:
    val chembl_url
    val chembl_version
    path chembl_db

    output:
    tuple val( chembl_version), path( "chembl_targets.tsv.gz" )

    script:
    if ( "${chembl_db}" == 'placeholder' ) {
        """
        set -x
        parse_json () (
            tr \$'\\t' '\\t' \
            | jq -r '
                .targets[] 
                | select( .species_group_flag? | not )
                | [
                    .tax_id, 
                    .organism, 
                    (
                        (
                            .target_components[0]
                            .target_component_synonyms
                            // empty
                        )
                        | map(select( .syn_type == "GENE_SYMBOL" ))
                        | first
                        | .component_synonym 
                        // "NA"
                    ), 
                    (
                        (
                            .target_components[0]
                            .target_component_synonyms
                            // empty
                        )
                        | map(select( .syn_type == "EC_NUMBER" ))
                        | first
                        | .component_synonym 
                        // "NA"
                    ), 
                    (
                        (.target_components[0].target_component_xrefs // empty) 
                        | map(select( .xref_src_db == "GoProcess" ))  
                        | map(.xref_id) | join(";") 
                        // "NA"
                    ),
                    (
                        (.target_components[0].target_component_xrefs // empty) 
                        | map(select( .xref_src_db == "GoProcess" ))  
                        | map(.xref_name) | join("; ") 
                        // "NA"
                    ),
                    .target_chembl_id, 
                    .target_components[0].accession, 
                    .pref_name
                ] 
                | @tsv
            '
        )

        root_url="${chembl_url}/chembl/api/data/target.json"
        query="target_type=SINGLE%20PROTEIN"
        init_url="\${root_url}?\${query}&limit=0"
        header=(target_taxon_id target_organism_name target_gene_symbol target_ec_number target_go_process_id target_go_process_name target_chembl_id target_uniprot_id target_name)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > chembl_targets.tsv

        curl -s "\${init_url}" > response.json

        parse_json < response.json >> chembl_targets.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt

        while [ "\$(cat next_page.txt)" != "null" ]
        do  
            sleep 0.3
            curl -s -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
                "${chembl_url}\$(cat next_page.txt)" > response.json
            parse_json < response.json >> chembl_targets.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt
        done    

        head -n1 chembl_targets.tsv \
        | cat - <(tail -n+2 chembl_targets.tsv | sort -u | sort -k1 ) \
        > chembl_targets-sorted.tsv
        
        mv chembl_targets-sorted.tsv chembl_targets.tsv
        gzip --best chembl_targets.tsv

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

        /*
         * Output columns:
         *  target_taxon_id
         *  target_organism_name
         *  target_gene_symbol
         *  target_ec_number
         *  target_go_process_id
         *  target_go_process_name
         *  target_chembl_id
         *  target_uniprot_id
         *  target_name
         */

        COPY (
          SELECT
            td.tax_id      AS target_taxon_id,
            td.organism    AS target_organism_name,

            -- first / representative gene symbol synonym
            max(CASE WHEN cs.syn_type = 'GENE_SYMBOL'
                     THEN cs.component_synonym END) AS target_gene_symbol,

            -- first / representative EC number synonym
            max(CASE WHEN cs.syn_type = 'EC_NUMBER'
                     THEN cs.component_synonym END) AS target_ec_number,

            -- GO “process” terms aggregated as in your JSON:
            -- IDs joined with ';', names joined with '; '
            string_agg(DISTINCT CASE
                                  WHEN cg.aspect = 'P'
                                  THEN cg.pref_name
                                END,
                       ';')     AS target_go_process_id,

            string_agg(DISTINCT CASE
                                  WHEN cg.aspect = 'P'
                                  THEN cg.pref_name
                                END,
                       '; ')    AS target_go_process_name,

            td.chembl_id   AS target_chembl_id,

            -- Uniprot accession (component_sequences)
            max(cseq.accession) AS target_uniprot_id,

            td.pref_name   AS target_name

          FROM target_dictionary      AS td
          JOIN target_type            AS tt   ON td.target_type = tt.target_type
          JOIN target_components      AS tc   ON td.tid = tc.tid
          LEFT JOIN component_sequences AS cseq
                 ON tc.component_id = cseq.component_id
          LEFT JOIN component_synonyms  AS cs
                 ON tc.component_id = cs.component_id
          LEFT JOIN (
              SELECT * 
              FROM component_go
              LEFT JOIN go_classification        AS gc
                 ON component_go.go_id = gc.go_id
            )
              AS cg
                 ON tc.component_id = cg.component_id
          

          WHERE
                td.target_type = 'SINGLE PROTEIN'
            AND tt.parent_type = 'PROTEIN'
            AND coalesce(td.species_group_flag, 0) = 0

          GROUP BY
            td.tax_id,
            td.organism,
            td.chembl_id,
            td.pref_name

          ORDER BY
            td.tax_id,
            td.chembl_id
        ) TO 'chembl_targets.tsv' (HEADER, DELIMITER '\\t');
        EOF

        # Optional: mimic your old sort-by-tax-id + uniq behaviour explicitly
        head -n1 chembl_targets.tsv \
          | cat - <(tail -n+2 chembl_targets.tsv | sort -u | sort -k1) \
          > chembl_targets-sorted.tsv

        mv chembl_targets-sorted.tsv chembl_targets.tsv
        gzip --best chembl_targets.tsv
        
        """
    }
    

}


process fetch_target_taxonomy {

    tag "v${chembl_version}"

    publishDir( 
        "${params.outputs}/taxonomy", 
        mode: 'copy',
        saveAs: { "${chembl_version}-${it}" }
    )
    
    input:
    val chembl_url
    val chembl_version
    path chembl_db

    output:
    path "taxon.tsv"

    script:
    if ( "${chembl_db}" == 'placeholder' ) {
        """
        set -x
        parse_json () (
            jq -r '.organisms[] | [.tax_id, .l1, .l2, .l3] | @tsv'
        )

        OUTFILE=taxon.tsv

        root_url="${chembl_url}/chembl/api/data/organism.json"
        url="\$root_url"'?limit=0'
        curl -s "\$url" > init_response.json

        header=(target_taxon_id target_taxon_l1 target_taxon_l2 target_taxon_l3)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > "\$OUTFILE"

        parse_json < init_response.json >> "\$OUTFILE"
        jq -r '.page_meta.next' < init_response.json > next_page.txt

        while [ "\$(cat next_page.txt)" != "null" ]
        do  
            sleep 0.3
            curl -s -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
                "${chembl_url}\$(cat next_page.txt)" > new_response.json
            parse_json < new_response.json >> "\$OUTFILE"
            jq -r '.page_meta.next' < new_response.json > next_page.txt
        done

        TEMP=\$(basename "\$OUTFILE" .tsv)-sorted.tsv
        head -n1 "\$OUTFILE" | cat - <(tail -n+2 "\$OUTFILE" | sort -u) > "\$TEMP" \
        && mv "\$TEMP" "\$OUTFILE"

        """
    }
    else {
        """
        set -euo pipefail

        OUTFILE=taxon.tsv
        echo "Using local ChEMBL SQLite DB: ${chembl_db}" >&2

        mkdir duckdb
        duckdb <<EOF
        SET home_directory='duckdb';

        INSTALL sqlite;
        LOAD sqlite;

        ATTACH '${chembl_db}' AS chembl (TYPE sqlite, READ_ONLY);
        USE chembl;

        COPY (
            SELECT DISTINCT
                tax_id AS target_taxon_id,
                l1     AS target_taxon_l1,
                l2     AS target_taxon_l2,
                l3     AS target_taxon_l3
            FROM organism_class
            ORDER BY tax_id
        ) TO 'taxon.tsv' (HEADER, DELIMITER '\\t');
        EOF

        TEMP=\$(basename "\$OUTFILE" .tsv)-sorted.tsv
        head -n1 "\$OUTFILE" | cat - <(tail -n+2 "\$OUTFILE" | sort -u) > "\$TEMP" \
        && mv "\$TEMP" "\$OUTFILE"

        """

    }

}


process fetch_chembl_inhibitors {

    tag "v${chembl_version}:${target_ids[0]}...${target_ids[-1]}: pChEMBL ≥ ${min_pchembl}"
    stageInMode 'link'
    // maxForks 2
    
    errorStrategy { if ( "${chembl_db}" == 'placeholder' ) { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' } else { return 'terminate' } }
    // errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    publishDir( 
        "${params.outputs}/inhibitors/by-target",
        mode: 'copy',
        saveAs: { "${target_ids[0]}-${target_ids[-1]}.${it}" },
    )
    
    input:
    val target_ids
    val chembl_url
    val chembl_version
    path chembl_db
    val min_pchembl

    output:
    tuple val( target_ids ), path( "inhibitors.tsv.gz" )

    script:
    if ( "${chembl_db}" == 'placeholder' ) {
        """
        set -euox pipefail

        fetch_json() (
            local url="\$1"
            local tries=8
            local delay=2
            for i in \$(seq 1 \$tries)
            do
                if curl -sS --fail-with-body \
                    -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
                    -H 'Accept: application/json' \
                    -D headers.txt \
                    -o response.json \
                    --connect-timeout 10 --max-time 120 \
                    "\$url"
                then
                    head -c1 response.json | grep -q '[\\{[]' && return 0
                fi
                echo "WARN: attempt \$i failed for \$url" >&2
                sleep "\$((delay ** i))"
            done
            echo "ERROR: giving up on \$url" >&2
            return 1
        )

        parse_assay () (
            jq -r '
                .assays[] | [
                    .assay_chembl_id, 
                    .assay_type,
                    .target_chembl_id,
                    .confidence_score,
                ] | @tsv
            '
        )

        parse_activity () (
            jq -r '
                .activities[] | [
                    .target_tax_id, 
                    .target_organism, 
                    .target_chembl_id, 
                    .target_pref_name, 
                    .molecule_chembl_id, 
                    .molecule_pref_name, 
                    .canonical_smiles
                ] | @tsv' \
            | sort -u
        )

        get_col_number () (
            head -n1 | tr \$'\\t' \$'\\n' | grep -nFw "\$1" | cut -d: -f1
        )

        # == Get assays
        root_url="${chembl_url}/chembl/api/data/assay.json"

        query="assay_type__in=F,B&confidence_score__gte=6"
        init_url="\${root_url}?\${query}&limit=0"

        header=(assay_chembl_id assay_type target_chembl_id assay_target_confidence_score)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > chembl_assays.tsv

        fetch_json "\$init_url"
        jq -r '.page_meta.next' < response.json > next_page.txt
        parse_assay < response.json >> chembl_assays.tsv

        while [ "\$(cat next_page.txt)" != "null" ]
        do  
            fetch_json "${chembl_url}\$(cat next_page.txt)"
            parse_assay < response.json >> chembl_assays.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt
        done    

        head -n1 chembl_assays.tsv \
        | cat - <(tail -n+2 chembl_assays.tsv | sort -u | sort -k1 ) \
        > chembl_assays-sorted.tsv \
        && mv chembl_assays-sorted.tsv chembl_assays.tsv


        # == Get activities
        root_url="${chembl_url}/chembl/api/data/activity.json"
        base_query="pchembl_value__gte=${min_pchembl}&potential_duplicate=0"

        header=(target_taxon_id target_organism target_chembl_id target_name molecule_chembl_id molecule_name molecule_smiles)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > chembl_ic50.tsv
        
        assay_id_col=\$(get_col_number assay_chembl_id < chembl_assays.tsv)
        tail -n+2 chembl_assays.tsv | cut -f"\$assay_id_col" | sort -u |split -l 20 - 'ids_'
        for id_file in ids_*
        do
            ids=\$(tr \$'\\n' , < "\$id_file")
            query="assay_chembl_id__in=\${ids}"
            init_url="\${root_url}?\${base_query}&\${query}&limit=0"

            fetch_json "\${init_url}"
            parse_activity < response.json >> chembl_ic50.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt

            while [ "\$(cat next_page.txt)" != "null" ]
            do  
                fetch_json "${chembl_url}\$(cat next_page.txt)"
                parse_activity < response.json >> chembl_ic50.tsv
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
            .drop_duplicates()
            .to_csv("inhibitors.tsv.gz", sep="\\t", index=False)
        )
        
        '

        """
    }
    else {
        """
        set -euo pipefail

        echo "Using local ChEMBL SQLite DB: ${chembl_db}" >&2

        mkdir duckdb
        duckdb <<EOF
        SET home_directory='duckdb';

        INSTALL sqlite;
        LOAD sqlite;

        ATTACH '${chembl_db}' AS chembl (TYPE sqlite, READ_ONLY);
        USE chembl;
        COPY (
            SELECT DISTINCT
                t.tax_id            AS target_taxon_id,
                t.organism          AS target_organism_name,
                t.chembl_id         AS target_chembl_id,
                t.pref_name         AS target_name,
                a.chembl_id         AS assay_chembl_id,
                a.assay_type        AS assay_type,
                a.confidence_score  AS assay_target_confidence_score,
                md.chembl_id        AS molecule_chembl_id,
                md.pref_name        AS molecule_name,
                cs.canonical_smiles AS molecule_smiles,
                act.pchembl_value   AS molecule_target_pchembl
            FROM activities            AS act
            JOIN assays                AS a   ON act.assay_id  = a.assay_id
            JOIN target_dictionary     AS t   ON a.tid         = t.tid
            JOIN molecule_dictionary   AS md  ON act.molregno  = md.molregno
            LEFT JOIN compound_structures AS cs ON md.molregno = cs.molregno
            WHERE
                    a.confidence_score >= 6
                AND act.pchembl_value  >= ${min_pchembl}
                AND COALESCE(act.potential_duplicate, 0) = 0
                AND t.chembl_id IN (${target_ids.collect { "'${it}'" }.join(',')})
        ) TO 'inhibitors.tsv' (HEADER, DELIMITER '\\t');
        EOF

        head -n1 inhibitors.tsv | cat - <(tail -n+2 inhibitors.tsv | sort -u) > inhibitors-sorted.tsv \
        && mv inhibitors-sorted.tsv inhibitors.tsv
        gzip --best inhibitors.tsv
        
        """

    }

}


process fetch_chembl_compound_mechanisms {

    tag "v${chembl_version}:${target_ids[0]}...${target_ids[-1]}"
    stageInMode 'link'
    // maxForks 2
    
    errorStrategy { if ( "${chembl_db}" == 'placeholder' ) { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' } else { return 'terminate' } }
    // errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    publishDir( 
        "${params.outputs}/mechanism/by-target",
        mode: 'copy',
        saveAs: { "${target_ids[0]}-${target_ids[-1]}.${it}" },
    )
    
    input:
    val target_ids
    val chembl_url
    val chembl_version
    path chembl_db

    output:
    tuple val( target_ids ), path( "compounds.tsv.gz" )

    script:
    if ( "${chembl_db}" == 'placeholder' ) {
        """
        set -euox pipefail
        UA='scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)'

        fetch_json() (
            local url="\$1"
            local tries=8
            local delay=2
            for i in \$(seq 1 \$tries)
            do
                if curl -sS --fail-with-body \
                    -A "\$UA" \
                    -H 'Accept: application/json' \
                    -D headers.txt \
                    -o response.json \
                    --connect-timeout 10 --max-time 120 \
                    "\$url"
                then
                    head -c1 response.json | grep -q '[\\{[]' && return 0
                fi
                echo "WARN: attempt \$i failed for \$url" >&2
                sleep "\$((delay ** i))"
            done
            echo "ERROR: giving up on \$url" >&2
            return 1
        )

        parse_json () (
            jq -r '
                .mechanisms[] | [
                    .target_chembl_id, 
                    .molecule_chembl_id,
                    .action_type,
                    .mechanism_of_action,
                    .max_phase,
                ] | @tsv' \
            | sort -u
        )

        root_url="${chembl_url}/chembl/api/data/mechanism.json"
        header=(target_chembl_id molecule_chembl_id parent_molecule_chembl_id action_type mechanism_of_action max_phase)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > mech.tsv

        query="target_chembl_id__in=${target_ids.join(",")}&direct_interaction=1&molecular_mechanism=1"
        init_url="\$root_url"'?limit=1000&'"\$query"
        
        fetch_json "\$init_url"
        jq -r '.page_meta.next' < response.json > next_page.txt
        parse_json < response.json >> inhibitors.tsv

        np=\$(cat next_page.txt)
        while [ "\$np" != "null" ]
        do  
            sleep 0.3
            fetch_json "${chembl_url}\$np"
            parse_json < response.json >> inhibitors.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt
            np=\$(cat next_page.txt)
        done

        head -n1 mech.tsv | cat - <(tail -n+2 mech.tsv | sort -u) > mech-sorted.tsv \
        && mv mech-sorted.tsv mech.tsv
        gzip --best mech.tsv

        # == Get all compounds from mech
        parse_mol () (
            jq -r '.molecules[] | [
                .molecule_chembl_id, 
                .pref_name,
                .molecule_structures.standard_inchi_key,
                .molecule_structures.smiles,
            ] | @tsv'
        )
        root_url="${chembl_url}/chembl/api/data/molecule.json"
        base_query=""

        header=(molecule_chembl_id molecule_name molecule_inchikey molecule_smiles)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > mol.tsv
        
        zcat mech.tsv.gz | tail -n+2 | while read target_line
        do
            sleep \$SLEEP_TIME
            target_id=\$(echo "\$target_line" | cut -f1)
            query="target_chembl_id=\${target_id}"
            init_url="\${root_url}?\${base_query}&\${query}&limit=0"

            curl -s "\${init_url}" > response.json

            parse_mol < response.json >> mol.tsv
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
                pd.read_csv("mech.tsv.gz", sep="\\t"),
                pd.read_csv("mol.tsv", sep="\\t"),
            )
            .drop_duplicates()
            .to_csv("compounds.tsv.gz", sep="\\t", index=False)
        )
        '

        """
    }
    else {
        """
        set -euo pipefail

        echo "Using local ChEMBL SQLite DB: ${chembl_db}" >&2

        mkdir duckdb
        duckdb <<EOF
        SET home_directory='duckdb';

        INSTALL sqlite;
        LOAD sqlite;

        ATTACH '${chembl_db}' AS chembl (TYPE sqlite, READ_ONLY);
        USE chembl;
        COPY (
            SELECT DISTINCT
                t.chembl_id              AS target_chembl_id,
                md.chembl_id             AS molecule_chembl_id,
                md.pref_name             AS molecule_name,
                cs.standard_inchi_key    AS molecule_inchikey,
                cs.canonical_smiles      AS molecule_smiles,
                mech.action_type         AS action_type,
                mech.mechanism_of_action AS mechanism_of_action,
                md.max_phase             AS max_phase
            FROM drug_mechanism        AS mech
            JOIN target_dictionary     AS t   ON mech.tid       = t.tid
            JOIN molecule_dictionary   AS md  ON mech.molregno  = md.molregno
            LEFT JOIN compound_structures cs
                ON md.molregno = cs.molregno
            WHERE
                    mech.direct_interaction = 1
                AND mech.molecular_mechanism = 1
                AND t.chembl_id IN (${target_ids.collect { "'${it}'" }.join(',')})
        ) TO 'compounds.tsv' (HEADER, DELIMITER '\\t');
        EOF

        head -n1 compounds.tsv | cat - <(tail -n+2 compounds.tsv | sort -u) > compounds-sorted.tsv \
        && mv compounds-sorted.tsv compounds.tsv
        gzip --best compounds.tsv
        
        """

    }

}


process fetch_pubchem_id {

    tag "${chembl_id[0]}...${chembl_id[-1]}:v${chembl_version}"
    stageInMode 'link'
    
    input:
    val chembl_id
    val chembl_url
    val chembl_version
    path chembl_db

    output:
    tuple val( chembl_id ), path( "inhibitors.tsv.gz" )

    script:
    if ( "${chembl_db}" == 'placeholder' ) {
        """
        set -euox pipefail
        
        parse_json () (
            jq -r '
                .molecules[] | [
                    .molecule_chembl_id, 
                    .pref_name, 
                    .molecule_structures.canonical_smiles, 
                    .molecule_structures.standard_inchi_key,
                    .oral, 
                    .topical, 
                    .parenteral, 
                    .orphan, 
                    .natural_product, 
                    .chemical_probe, 
                    .black_box_warning, 
                    .max_phase
                ] | @tsv' \
            | sort -u
        )

        root_url="${chembl_url}/chembl/api/data/molecule.json"
        query="molecule_chembl_id__in=${chembl_id.join(",")}"
        init_url="\$root_url"'?limit=0&'"\$query"
        
        header=(molecule_chembl_id molecule_name molecule_smiles molecule_inchikey is_oral is_topical is_parenteral is_orphan is_natural_product is_chemcial_probe has_black_box max_phase)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > inhibitors.tsv

        curl -s "\$init_url" \
            -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
        > init_response.json
        jq -r '.page_meta.next' < init_response.json > next_page.txt
        parse_json < init_response.json >> inhibitors.tsv

        while [ "\$(cat next_page.txt)" != "null" ]
        do  
            sleep 0.3
            curl -s -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
                "${chembl_url}\$(cat next_page.txt)" \
                > new_response.json
            parse_json < new_response.json >> inhibitors.tsv
            jq -r '.page_meta.next' < new_response.json > next_page.txt
        done

        head -n1 inhibitors.tsv | cat - <(tail -n+2 inhibitors.tsv | sort -u) > inhibitors-sorted.tsv
        mv inhibitors-sorted.tsv inhibitors.tsv
        gzip --best inhibitors.tsv
        
        """
    } 
    else {
        """
        set -euox pipefail

        echo "Using local ChEMBL SQLite DB: ${chembl_db}" >&2

        # Pull molecule metadata from local ChEMBL DB
        mkdir duckdb
        duckdb :memory: << EOF
        SET home_directory='duckdb';

        INSTALL sqlite;
        LOAD sqlite;

        ATTACH '${chembl_db}' AS chembl (TYPE sqlite, READ_ONLY);
        USE chembl;

        COPY (
            SELECT
                md.chembl_id          AS molecule_chembl_id,
                md.pref_name          AS molecule_name,
                cs.canonical_smiles   AS molecule_smiles,
                cs.standard_inchi_key AS molecule_inchikey,
                md.oral               AS is_oral,
                md.topical            AS is_topical,
                md.parenteral         AS is_parenteral,
                md.orphan             AS is_orphan,
                md.natural_product    AS is_natural_product,
                md.chemical_probe     AS is_chemcial_probe,
                md.black_box_warning  AS has_black_box,
                md.max_phase          AS max_phase
            FROM molecule_dictionary md
            LEFT JOIN compound_structures cs
                ON md.molregno = cs.molregno
            WHERE md.chembl_id IN (
                ${chembl_id.collect { "'${it}'" }.join(',')}
            )
        ) TO 'inhibitors.tsv' (HEADER, DELIMITER '\\t');
        EOF

        head -n1 inhibitors.tsv | cat - <(tail -n+2 inhibitors.tsv | sort -u) > inhibitors-sorted.tsv
        mv inhibitors-sorted.tsv inhibitors.tsv
        gzip --best inhibitors.tsv

        """
    }


}


process fetch_vendors {

    tag "${chembl_id[0]}...${chembl_id[-1]}"
    stageInMode 'link'

    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5
    
    input:
    tuple val( chembl_id ), path( table )

    output:
    tuple val( chembl_id ), path( "purchasable.tsv.gz" )

    script:

    """
    set -eux

    parse_json_unichem () (
        jq -r '
        # helper: pick first source with given shortName, or {} if none
        def pick(\$name):
            ((map(select(.shortName == \$name)) | first) // {})
            | (.compoundId // "NA", .url // "NA");

        # start from sources; if no compound or no sources, use [] so map() is safe
        (.compounds[0].sources? // []) as \$srcs
        | [
            (.compounds[0].standardInchiKey // "NA"),
            (\$srcs | pick("chembl")),
            (\$srcs | pick("pubchem")),
            (\$srcs | pick("drugbank")),
            (\$srcs | pick("zinc")),
            (\$srcs | pick("emolecules")),
            (\$srcs | pick("selleck")),
            (\$srcs | pick("mcule")),
            (\$srcs | pick("molport")),
            (\$srcs | pick("MedChemExpress"))
        ]
        | @tsv
        '
    )

    inchikey_col=\$(zcat "${table}" | head -n1 | tr \$'\\t' \$'\\n' | grep -n -Fw "molecule_inchikey" | cut -d: -f1)

    header=(molecule_inchikey molecule_chembl_id_2 chembl_url pubchem_id pubchem_url drugbank_id drugbank_url vendor_zinc_id zinc_url vendor_emolecules_id emolecules_url vendor_selleck selleck_url vendor_mcule mcule_url vendor_molport molport_url vendor_mce mce_url)
    
    printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
    > pubchem_ids.txt
    for key in \$(zcat "${table}" | tail -n+2 | cut -f"\$inchikey_col")
    do
        sleep 0.1
        curl -s --request POST \
            -H "accept: application/json" \
            -H "Content-Type: application/json" \
            -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
            --url https://www.ebi.ac.uk/unichem/api/v1/compounds \
            --data '{
                "type": "inchikey",
                "compound": "'"\$key"'"
            }' \
        > unichem-response.json

        parse_json_unichem < unichem-response.json \
        >> pubchem_ids.txt
    done

    python -c '
    import pandas as pd
    NA = "NA"
    (
        pd.read_csv("'"${table}"'", sep="\\t")
        .merge(
            pd.read_csv(
                "pubchem_ids.txt", 
                sep="\\t",
            )
            .query("not molecule_inchikey.isna() and molecule_inchikey != @NA"),
            how="left",
        )
        .to_csv(
            "purchasable.tsv", 
            sep="\\t", 
            index=False,
        )
    )
    '

    gzip --best purchasable.tsv
    
    """

}
