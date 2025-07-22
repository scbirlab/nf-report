process chembl_status {

    tag "${date}"

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
    curl -s https://www.ebi.ac.uk/chembl/api/data/status.json > chembl_version.json
    chembl_version=\$(jq -r '.chembl_db_version' < chembl_version.json)
    echo \$chembl_version

    """

}



process fetch_chembl_target_sequences {

    tag "v${chembl_version}"

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

    output:
    path "chembl_tox.tsv"

    script:
    """
    set -x
    # == config 
    SLEEP_TIME=0.3

    parse_assay () (
        jq -r '.assays[] | [.cell_chembl_id, .assay_cell_type, .assay_chembl_id, .assay_type] | @tsv'
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
        jq -r '.mechanisms[] | [.molecule_chembl_id, .target_chembl_id, .mechanism_of_action, .max_phase] | @tsv' 
    )

    get_col_number () (
        head -n1 | tr \$'\\t' \$'\\n' | grep -nFw "\$1" | cut -d: -f1
    )

    # == Get all cell line assays
    root_url="${chembl_url}/chembl/api/data/assay.json"

    query="assay_cell_type__in=${cell_ids.join(',')}&assay_type__in=F,T&cell_chembl_id__isnull=False"
    init_url="\${root_url}?\${query}&limit=0"

    printf "cell_chembl_id\\tcell_type\\tassay_chembl_id\\tassay_type\\n" \
    > chembl_assays.tsv

    curl -s "\${init_url}" > response.json

    parse_assay < response.json >> chembl_assays.tsv
    jq -r '.page_meta.next' < response.json > next_page.txt

    while [ "\$(cat next_page.txt)" != "null" ]
    do  
        sleep \$SLEEP_TIME
        curl -s "${chembl_url}\$(cat next_page.txt)" > response.json
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

    printf "assay_chembl_id\\tmolecule_chembl_id\\tsmiles\\tmeasurement_type\\tstandard_value\\tunits\\n" \
    > chembl_ic50.tsv
    
    assay_id_col=\$(get_col_number assay_chembl_id < chembl_assays.tsv)
    tail -n+2 chembl_assays.tsv | cut -f"\$assay_id_col" | sort -u |split -l 20 - 'ids_'
    for id_file in ids_*
    do
        sleep \$SLEEP_TIME
        ids=\$(tr \$'\\n' , < "\$id_file")
        query="assay_chembl_id__in=\${ids}"
        init_url="\${root_url}?\${base_query}&\${query}&limit=0"

        curl -s "\${init_url}" > response.json

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

    printf "molecule_chembl_id\\ttarget_chembl_id\\tmechanism\\tmax_phase\\n" \
    > targets.tsv
    
    mol_col=\$(get_col_number molecule_chembl_id < chembl_ic50.tsv)
    tail -n+2 chembl_ic50.tsv | cut -f"\$mol_col" | sort -u | split -l 20 - 'mols_'
    for id_file in mols_*
    do
        sleep \$SLEEP_TIME
        ids=\$(tr \$'\\n' , < "\$id_file")
        query="molecule_chembl_id__in=\${ids}"
        init_url="\${root_url}?\${base_query}&\${query}&limit=0"

        curl -s "\${init_url}" > response.json

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
    printf "molecule_chembl_id\\ttarget_chembl_id\\ttarget_measurement\\ttarget_inhibition_value\\ttarget_inhibition_units\\n" \
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
            curl -s "${chembl_url}\$(cat next_page.txt)" > response.json
            parse_inhibition < response.json >> inhibition.tsv
            jq -r '.page_meta.next' < response.json > next_page.txt
        done  
    done

    python -c '
    import pandas as pd

    pd.merge(
        pd.read_csv("chembl_assays.tsv", sep="\\t"),
        pd.read_csv("chembl_ic50.tsv", sep="\\t"),
    ).merge(
        pd.read_csv("targets.tsv", sep="\\t")
    ).merge(
        pd.read_csv("inhibition.tsv", sep="\\t")
    ).drop_duplicates().to_csv("chembl_tox.tsv", sep="\\t", index=False)
    
    '
    """

}


process fetch_chembl_targets {

    tag "v${chembl_version}:${chembl_url}"

    publishDir( 
        "${params.outputs}/targets", 
        mode: 'copy',
        saveAs: { "${it}" }
    )
    
    input:
    val chembl_url
    val chembl_version

    output:
    path "chembl_targets.tsv"

    script:
    """
    set -x
    parse_json () (
        tr \$'\\t' '\\t' \
        | jq -r '
            .targets[] 
            | [
                .tax_id, 
                .organism, 
                (
                    .target_components[0]
                    .target_component_synonyms
                    | map(select( .syn_type == "GENE_SYMBOL" ))
                    | first
                    | .component_synonym 
                    // "NA"
                ), 
                (
                    .target_components[0]
                    .target_component_synonyms
                    | map(select( .syn_type == "EC_NUMBER" ))
                    | first
                    | .component_synonym 
                    // "NA"
                ), 
                .target_chembl_id, 
                .target_components[0].accession, 
                .target_components[0].component_description
            ] 
            | @tsv
        '
    )

    root_url="${chembl_url}/chembl/api/data/target.json"
    query="target_type=SINGLE%20PROTEIN&confidence_score__gte=6&pchembl_value__gte=5&potential_duplicate=0"
    init_url="\${root_url}?\${query}&limit=0"

    printf "taxon_id\\torganism_name\\tgene\\tec_number\\ttarget_chembl_id\\ttarget_uniprot_id\\ttarget_name\\n" \
    > chembl_targets.tsv

    curl -s "\${init_url}" > response.json

    parse_json < response.json >> chembl_targets.tsv
    jq -r '.page_meta.next' < response.json > next_page.txt

    while [ "\$(cat next_page.txt)" != "null" ]
    do  
        sleep 0.3
        curl -s "${chembl_url}\$(cat next_page.txt)" > response.json
        parse_json < response.json >> chembl_targets.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt
    done    

    head -n1 chembl_targets.tsv \
    | cat - <(tail -n+2 chembl_targets.tsv | sort -u | sort -k1 ) \
    > chembl_targets-sorted.tsv \
    && mv chembl_targets-sorted.tsv chembl_targets.tsv

    """

}


process fetch_included_chembl_taxon {

    tag "v${chembl_version}:${id}:${organism_id}"

    publishDir( 
        "${params.outputs}/targets", 
        mode: 'copy',
        saveAs: { "${organism_id}-${it}" }
    )
    
    input:
    tuple val( id ), val( organism_id )
    val chembl_url
    val chembl_version
    val default_kingdom

    output:
    tuple val( id ), path( "included_taxid.tsv" )

    script:
    """
    set -x
    parse_json () (
        jq -r '.organisms[] | [.tax_id, .l1, .l2, .l3] | @tsv'
    )

    root_url="${chembl_url}/chembl/api/data/organism.json"
    query="tax_id=${organism_id}"
    init_url="\$root_url"?"\$query"
    curl -s "\$init_url" | jq -r '.organisms[0].l1 // empty' > kingdom.txt
    if [ ! -s kingdom.txt ]
    then
        echo "${default_kingdom}" > kingdom.txt
    fi

    #kingdom=\$(cat kingdom.txt)
    #url="\$root_url"'?limit=0&l1__regex=^(?!'\$kingdom'\$).*'
    url="\$root_url"'?limit=0'
    curl -s "\$url" > init_response.json
    printf "taxon_id\\ttaxon_l1\\ttaxon_l2\\ttaxon_l3\\n" > included_taxid.tsv

    parse_json < init_response.json >> included_taxid.tsv
    jq -r '.page_meta.next' < init_response.json > next_page.txt

    while [ "\$(cat next_page.txt)" != "null" ]
    do  
        sleep 0.3
        curl -s "${chembl_url}\$(cat next_page.txt)" > new_response.json
        parse_json < new_response.json >> included_taxid.tsv
        jq -r '.page_meta.next' < new_response.json > next_page.txt
    done

    head -n1 included_taxid.tsv | cat - <(tail -n+2 included_taxid.tsv | sort -u) > included_taxid-sorted.tsv \
    && mv included_taxid-sorted.tsv included_taxid.tsv

    """

}


process filter_targets_by_taxon {

    tag "${id}"

    publishDir( 
        "${params.outputs}/targets", 
        mode: 'copy',
        saveAs: { "${id}-${it}" }
    )
    
    input:
    tuple val( id ), path( included_taxid ), path( all_targets ) 

    output:
    tuple val( id ), path( "filtered_targets.tsv" )

    script:
    """
    set -x
    #awk -v OFS='\\t' 'NR > 1 { print "^"\$0,"" }' "${included_taxid}" > search.txt
    #grep -f search.txt "${all_targets}" > table-tail.tsv
    #head -n1 "${all_targets}" \
    #| cat - table-tail.tsv \
    #> filtered_targets.tsv
    # head -n1 "${all_targets}" | cat - <(tail -n+2 "${all_targets}" | sort -u | sort -k1 -n ) \
    # > sorted.tsv
    # join -t\$'\t' --nocheck-order --header -j1 \
    #     "${included_taxid}" sorted.tsv \
    # > filtered_targets.tsv

    python -c '
    import pandas as pd

    pd.merge(
        pd.read_csv("${included_taxid}", sep="\\t"),
        pd.read_csv("${all_targets}", sep="\\t"),
    ).to_csv("filtered_targets.tsv", sep="\\t", index=False)
    '

    """

}


process fetch_chembl_inhibitors {

    tag "v${chembl_version}:${id}:${target_id}"

    errorStrategy 'retry'
    maxRetries 2
    
    input:
    tuple val( id ), val( target_id )
    val chembl_url
    val chembl_version

    output:
    tuple val( id ), path( "*.tsv" )

    script:
    """
    set -x
    parse_json () (
        jq -r '.activities[] | [.target_chembl_id, .target_pref_name, .target_organism, .target_tax_id, .molecule_chembl_id, .molecule_pref_name, .canonical_smiles] | @tsv' \
        | sort -u
    )

    root_url="${chembl_url}/chembl/api/data/activity.json"
    query="target_chembl_id=${target_id}&confidence_score__gte=6&pchembl_value__gte=5&potential_duplicate=0"
    init_url="\$root_url"'?limit=0&'"\$query"
    
    printf "target_chembl_id\\ttarget_name\\torganism_name\\ttaxon_id\\tmolecule_chembl_id\\tmolecule_name\\tsmiles\\n" \
    > inhibitors.tsv

    curl -s "\$init_url" > response.json
    jq -r '.page_meta.next' < response.json > next_page.txt
    parse_json < response.json >> inhibitors.tsv

    np=\$(cat next_page.txt)
    while [ "\$np" != "null" ]
    do  
        sleep 0.3
        curl -s "${chembl_url}\$np" > response.json
        parse_json < response.json >> inhibitors.tsv
        jq -r '.page_meta.next' < response.json > next_page.txt
        np=\$(cat next_page.txt)
    done

    head -n1 inhibitors.tsv | cat - <(tail -n+2 inhibitors.tsv | sort -u) > inhibitors-sorted.tsv \
    && mv inhibitors-sorted.tsv inhibitors.tsv

    """

}


process fetch_chembl_inhibitor_activities {

    tag "v${chembl_version}:${id}:${target_id}: pChEMBL ≥ ${min_pchembl}"

    errorStrategy 'retry'
    maxRetries 2
    
    input:
    tuple val( id ), val( target_id )
    val chembl_url
    val chembl_version
    val min_pchembl

    output:
    tuple val( id ), path( "*.tsv" )

    script:
    """
    set -x
    parse_json () (
        jq -r '.activities[] | [.target_chembl_id, .target_pref_name, .target_organism, .target_tax_id, .molecule_chembl_id, .molecule_pref_name, .canonical_smiles] | @tsv' \
        | sort -u
    )

    root_url="${chembl_url}/chembl/api/data/activity.json"
    base_query='confidence_score__gte=6&potential_duplicate=0'
    query="target_chembl_id=${target_id}&pchembl_value__gte=${min_pchembl}"
    init_url="\$root_url"'?limit=0&'"\$base_query"'&'"\$query"
    
    printf "target_chembl_id\ttarget_name\torganism_name\ttaxon_id\tmolecule_chembl_id\tmolecule_name\tsmiles\n" \
    > inhibitors.tsv

    curl -s "\$init_url" > response.json
    jq -r '.page_meta.next' < response.json > next_page.txt
    parse_json < response.json >> inhibitors.tsv

    np=\$(cat next_page.txt)
    while [ "\$np" != "null" ]
    do  
        sleep 0.3
        curl -s "${chembl_url}\$np" > response.json
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

    errorStrategy 'retry'
    maxRetries 2
    
    input:
    val chembl_id
    val chembl_url
    val chembl_version

    output:
    tuple val( chembl_id ), path( "inhibitors.tsv" )

    script:
    """
    set -x
    parse_json () (
        jq -r '.molecules[] | [.molecule_chembl_id, .pref_name, .oral, .topical, .parenteral, .orphan, .natural_product, .chemical_probe, .black_box_warning, .max_phase, .molecule_structures.canonical_smiles, .molecule_structures.standard_inchi_key] | @tsv' \
        | sort -u
    )

    root_url="${chembl_url}/chembl/api/data/molecule.json"
    query="molecule_chembl_id=${chembl_id}"
    init_url="\$root_url"'?limit=0&'"\$query"
    
    printf "molecule_chembl_id\tmolecule_name\tis_oral\tis_topical\tis_parenteral\tis_orphan\tis_natural_product\tis_chemcial_probe\thas_black_box\tmax_phase\tsmiles\tinchikey\n" \
    > inhibitors.tsv

    curl -s "\$init_url" > init_response.json
    jq -r '.page_meta.next' < init_response.json > next_page.txt
    parse_json < init_response.json >> inhibitors.tsv

    while [ "\$(cat next_page.txt)" != "null" ]
    do  
        sleep 0.3
        curl -s "${chembl_url}\$(cat next_page.txt)" > new_response.json
        parse_json < new_response.json >> inhibitors.tsv
        jq -r '.page_meta.next' < new_response.json > next_page.txt
    done

    inchikey_col=\$(head -n1 inhibitors.tsv | tr \$'\\t' \$'\\n' | grep -n -Fw inchikey | cut -d: -f1)
    printf "pubchem_id\\tpubchem_url\\tmolecule_chembl_id\\tmolecule_chembl_url\\tdrugbank_id\\tdrugbank_url\\tvendor_zinc_id\\tzinc_url\\tvendor_emolecules_id\\temolecules_url\\tvendor_selleck\\tselleck_url\\tvendor_mcule\\tmcule_url\\tvendor_molport\\tmolport_url\\tvendor_mce\\tmce_url\\n" > pubchem_ids.txt
    for key in \$(tail -n+2 inhibitors.tsv | cut -f\$inchikey_col)
    do
        curl -s --request POST \
            -H "accept: application/json" -H "Content-Type: application/json" \
            --url https://www.ebi.ac.uk/unichem/api/v1/compounds \
            --data '{
                "type": "inchikey",
                "compound": "'"\$key"'"
            }' > unichem-response.json
        jq -r '
            .compounds[0].sources 
            | [
                (map(select( .shortName == "pubchem" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "chembl" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "drugbank" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "zinc" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "emolecules" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "selleck" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "mcule" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "molport" )) | first | (.compoundId // "NA", .url // "NA")),
                (map(select( .shortName == "MedChemExpress" )) | first | (.compoundId // "NA", .url // "NA"))
            ] | @tsv
        ' \
        < unichem-response.json \
        >> pubchem_ids.txt
    done

    paste inhibitors.tsv pubchem_ids.txt > inhibitors-ids.tsv
    head -n1 inhibitors-ids.tsv | cat - <(tail -n+2 inhibitors-ids.tsv | sort -u) > inhibitors-sorted.tsv
    mv inhibitors-sorted.tsv inhibitors.tsv
    
    """

}