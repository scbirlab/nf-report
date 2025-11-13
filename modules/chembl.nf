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
    
    errorStrategy { sleep(Math.pow(2, task.attempt) * 60000 as long); return 'retry' }
    maxRetries 5

    input:
    val chembl_url
    val cell_ids
    val chembl_version
    path chembl_db

    output:
    path "chembl_tox.tsv", emit: main
    path "*.tsv", emit: tables

    script:
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
        .to_csv("chembl_tox.tsv", sep="\\t", index=False)
    )
    
    '

    """

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
    tuple val( chembl_version), path( "chembl_targets.tsv" )

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

    tag "v${chembl_version}:${id}:${target_ids[0]}...${target_ids[-1]}: pChEMBL ≥ ${min_pchembl}"
    stageInMode 'link'
    // maxForks 2
    
    errorStrategy { if ( "${chembl_db}" == 'placeholder' ) { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' } else { return 'terminate' } }
    // errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    publishDir( 
        "${params.outputs}/inhibitors/by-target", 
        mode: 'copy',
        saveAs: { "${id}.${target_ids[0]}-${target_ids[-1]}.${it}" },
    )
    
    input:
    tuple val( id ), val( target_ids )
    val chembl_url
    val chembl_version
    path chembl_db
    val min_pchembl

    output:
    tuple val( id ), path( "inhibitors.tsv" )

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

        root_url="${chembl_url}/chembl/api/data/activity.json"
        header=(target_taxon_id target_organism_name target_chembl_id target_name molecule_chembl_id molecule_name molecule_smiles)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > inhibitors.tsv

        query="target_chembl_id__in=${target_ids.join(",")}&confidence_score__gte=6&pchembl_value__gte=${min_pchembl}&potential_duplicate=0"
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

        head -n1 inhibitors.tsv | cat - <(tail -n+2 inhibitors.tsv | sort -u) > inhibitors-sorted.tsv \
        && mv inhibitors-sorted.tsv inhibitors.tsv

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
                md.chembl_id        AS molecule_chembl_id,
                md.pref_name        AS molecule_name,
                cs.canonical_smiles AS molecule_smiles
            FROM activities            AS act
            JOIN assays                AS a   ON act.assay_id  = a.assay_id
            JOIN target_dictionary     AS t   ON a.tid         = t.tid
            JOIN molecule_dictionary   AS md  ON act.molregno  = md.molregno
            LEFT JOIN compound_structures AS cs ON md.molregno = cs.molregno
            WHERE
                    a.confidence_score        >= 6
                AND act.pchembl_value         >= ${min_pchembl}
                AND COALESCE(act.potential_duplicate, 0) = 0
                AND t.chembl_id IN (${target_ids.collect { "'${it}'" }.join(',')})
        ) TO 'inhibitors.tsv' (HEADER, DELIMITER '\\t');
        EOF

        head -n1 inhibitors.tsv | cat - <(tail -n+2 inhibitors.tsv | sort -u) > inhibitors-sorted.tsv \
        && mv inhibitors-sorted.tsv inhibitors.tsv
        
        """

    }

}


process fetch_chembl_inhibitor_activities {

    tag "v${chembl_version}:${id}:${target_id}: pChEMBL ≥ ${min_pchembl}"
    stageInMode 'link'

    errorStrategy 'retry'
    maxRetries 2
    
    input:
    tuple val( id ), val( target_id )
    val chembl_url
    val chembl_version
    path chembl_db
    val min_pchembl

    output:
    tuple val( id ), path( "*.tsv" )

    script:
    """
    set -x
    parse_json () (
        jq -r '.activities[] | [
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

    root_url="${chembl_url}/chembl/api/data/activity.json"
    base_query='confidence_score__gte=6&potential_duplicate=0'
    query="target_chembl_id=${target_id}&pchembl_value__gte=${min_pchembl}"
    init_url="\$root_url"'?limit=0&'"\$base_query"'&'"\$query"

    header=(target_taxon_id target_organism_name target_chembl_id target_name molecule_chembl_id molecule_name molecule_smiles)
    
    printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
    > inhibitors.tsv

    curl -s "\$init_url" > response.json
    jq -r '.page_meta.next' < response.json > next_page.txt
    parse_json < response.json >> inhibitors.tsv

    np=\$(cat next_page.txt)
    while [ "\$np" != "null" ]
    do  
        sleep 0.3
        curl -s -A 'scbirlab-nf-report/0.4 (+https://scbirlab.org; contact: eachan.johnson@crick.ac.uk)' \
            "${chembl_url}\$np" > response.json
        parse_json < response.json >> inhibitors.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt
        np=\$(cat next_page.txt)
    done

    head -n1 inhibitors.tsv | cat - <(tail -n+2 inhibitors.tsv | sort -u) > inhibitors-sorted.tsv \
    && mv inhibitors-sorted.tsv inhibitors.tsv

    """

}


process fetch_pubchem_id {

    tag "${chembl_id}:v${chembl_version}"
    stageInMode 'link'

    errorStrategy { "${chembl_db}" == 'placeholder' ? 'retry' : 'terminate' }
    maxRetries 2
    
    input:
    val chembl_id
    val chembl_url
    val chembl_version
    path chembl_db

    output:
    tuple val( chembl_id ), path( "inhibitors.tsv" )

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

        parse_json_unichem () (
            jq -r '
                .compounds[0].sources 
                | [
                    (map(select( .shortName == "chembl" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "pubchem" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "drugbank" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "zinc" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "emolecules" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "selleck" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "mcule" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "molport" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "MedChemExpress" )) | first | (.compoundId // "NA", .url // "NA"))
                ] | @tsv'
        )

        root_url="${chembl_url}/chembl/api/data/molecule.json"
        query="molecule_chembl_id=${chembl_id}"
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

        inchikey_col=\$(head -n1 inhibitors.tsv | tr \$'\\t' \$'\\n' | grep -n -Fw molecule_inchikey | cut -d: -f1)

        header=(molecule_chembl_id molecule_chembl_url pubchem_id pubchem_url drugbank_id drugbank_url vendor_zinc_id zinc_url vendor_emolecules_id emolecules_url vendor_selleck selleck_url vendor_mcule mcule_url vendor_molport molport_url vendor_mce mce_url)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > pubchem_ids.txt
        for key in \$(tail -n+2 inhibitors.tsv | cut -f"\$inchikey_col")
        do
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

        paste inhibitors.tsv pubchem_ids.txt > inhibitors-ids.tsv
        head -n1 inhibitors-ids.tsv | cat - <(tail -n+2 inhibitors-ids.tsv | sort -u) > inhibitors-sorted.tsv
        mv inhibitors-sorted.tsv inhibitors.tsv
        
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
            WHERE md.chembl_id = '${chembl_id}'
        ) TO 'inhibitors.tsv' (HEADER, DELIMITER '\\t');
        EOF

        parse_json_unichem () (
            jq -r '
                .compounds[0].sources 
                | [
                    (map(select( .shortName == "chembl" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "pubchem" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "drugbank" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "zinc" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "emolecules" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "selleck" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "mcule" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "molport" )) | first | (.compoundId // "NA", .url // "NA")),
                    (map(select( .shortName == "MedChemExpress" )) | first | (.compoundId // "NA", .url // "NA"))
                ] | @tsv'
        )

        inchikey_col=\$(head -n1 inhibitors.tsv | tr \$'\\t' \$'\\n' | grep -n -Fw molecule_inchikey | cut -d: -f1)

        header=(molecule_chembl_id molecule_chembl_url pubchem_id pubchem_url drugbank_id drugbank_url vendor_zinc_id zinc_url vendor_emolecules_id emolecules_url vendor_selleck selleck_url vendor_mcule mcule_url vendor_molport molport_url vendor_mce mce_url)
        
        printf "\$(IFS=\$'\\t'; echo "\${header[*]}")\\n" \
        > pubchem_ids.txt
        for key in \$(tail -n+2 inhibitors.tsv | cut -f"\$inchikey_col")
        do
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

        paste inhibitors.tsv pubchem_ids.txt > inhibitors-ids.tsv
        head -n1 inhibitors-ids.tsv | cat - <(tail -n+2 inhibitors-ids.tsv | sort -u) > inhibitors-sorted.tsv
        mv inhibitors-sorted.tsv inhibitors.tsv

        """
    }


}
