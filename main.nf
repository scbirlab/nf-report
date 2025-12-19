#!/usr/bin/env nextflow

/*
========================================================================================
   Pipeline to identify orthologs that have purchasable inhibitors
========================================================================================
   Github   : https://github.com/scbirlab/nf-reclaim
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

pipeline_title = """\
   R E P O R T   P I P E L I N E
   =========================================================================
   Nextflow pipeline to identify orthologs that have purchasable inhibitors.
   
   """
   .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println pipeline_title + """\
         Command-line usage:
            nextflow run scbirlab/nf-report --organism_id <taxon ID>
         Config/sample sheet usage:
            nextflow run scbirlab/nf-report -c <config-file>

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """
   .stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   if ( !params.organism_id ) {
      throw new Exception("!!! PARAMETER MISSING: Please provide a sample sheet or at least --organism_id.")
   }
}

log.info pipeline_title + """\
   test mode               : ${params.test}
   inputs
      sample sheet         : ${params.sample_sheet}
      Taxon ID             : ${params.organism_id}
      Chembl URL           : ${params.chembl_url}
      TaxonKit DB          : ${params.taxonkit_db_url}
      gNOMAD version       : ${params.gnomad_version}
   parameters
      min. pChEMBL         : ${params.min_pchembl}
      min. LOEUF           : ${params.min_loeuf}
      min. ID              : ${params.min_identity}
      min. coverage        : ${params.min_coverage}
      batch size           : ${params.batch_size}
      fetch inhibitors?    : ${params.inhibitors}
   output                  : ${params.outputs}
      make plots?          : ${params.plots}
   """
   .stripIndent()


/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

// load modules
include { 
   chembl_status;
   fetch_chembl_inhibitors;
   fetch_chembl_targets;
   fetch_chembl_tox;
   fetch_pubchem_id;
   fetch_target_taxonomy;
   fetch_vendors;
} from './modules/chembl.nf'
include { 
   make_diamond_db;
   diamond_blastp;
} from './modules/diamond.nf'
include { 
   fetch_gnomad_constraints;
} from './modules/gnomad.nf'
include { 
   factorise_nmf;
} from './modules/nmf.nf'
include { 
   describe;
   find_coverage_cutoff;
} from './modules/orthology.nf'
include { 
   describe_provenance;
} from './modules/provenance.nf'
include { 
   silhouette;
} from './modules/silhouette.nf'
include { 
   fetch_taxonkit_db;
   fetch_taxonomic_ranks;
} from './modules/taxonkit.nf'
include { 
   umaps_of_rbh_matrix;
} from './modules/umap.nf'
include { 
   fetch_fastas_from_organism_id;
   fetch_fastas_from_uniprot_ids;
   fetch_species_gene_names;
} from './modules/uniprot.nf'
include { 
   stack_tables;
   stack_tables as stack_tables2;
   stack_tables as stack_tables3;
   stack_tables as stack_tables4;
   subset_table;
   filter_target_list;
   make_rbh_matrix;
   merge_tox_gnomad;
   merge_tables;
   merge_tables as merge_tables2;
   merge_tables as merge_tables3;
   merge_tables as merge_tables4;
} from './modules/utils.nf'



workflow {

   if ( !params.chembl_db ) {
      chembl_status(
         Channel.of( workflow.start )
      )
      chembl_status.out.version
         .first()
         .set { chembl_version }
   } else {
      Channel.value( file( params.chembl_db ).simpleName )
         .set { chembl_version }
   }

   Channel.value( 
      params.chembl_db 
      ? file( params.chembl_db, checkIfExists: true ) 
      : file( 'placeholder' ) 
   )
      .set { chembl_db }

   Channel.value( params.chembl_url )
      .set { chembl_url }

   if ( params.sample_sheet ) {

      Channel.fromPath( 
         params.sample_sheet,
         checkIfExists: true, 
      )
         .set { sample_sheet_ch }
         

   }

   else {

      Channel.of( "organism_id", "${params.organism_id}" )
         .collectFile( 
            name: "sample-sheet.csv", 
            keepHeader: true, 
            newLine: true, 
            storeDir: "${params.outputs}/sample-sheet",
         )
         .set { sample_sheet_ch }

   }

   sample_sheet_ch
      .splitCsv( header: true )
      .set { sample_rows }

   fetch_taxonkit_db(
      Channel.of( params.taxonkit_db_url ),
   )

   fetch_taxonomic_ranks(
      sample_sheet_ch.combine( fetch_taxonkit_db.out ),
      Channel.value( "organism_id" ),
   )

   fetch_target_taxonomy(
      chembl_url,
      chembl_version,
      chembl_db,
   )
   fetch_chembl_targets(
      chembl_url,
      chembl_version,
      chembl_db,
   )  
   fetch_gnomad_constraints(
      Channel.of( params.gnomad_version ),
   )

   if ( !params.test && params.fetch_tox ) {

      fetch_chembl_tox(
         chembl_url,
         Channel.value( params.tox_cell_lines ),
         chembl_version,
         chembl_db,
      )

      merge_tox_gnomad(
         fetch_chembl_targets.out
            .map { v -> v[-1] },
         fetch_chembl_tox.out.main,
         fetch_gnomad_constraints.out,
         Channel.value( "inner" ),
      )
      
   }

   fetch_chembl_targets.out
      .map { it[-1] }
      .splitCsv( 
         header: true, 
         elem: 0, //1, 
         sep: '\t',
      )
      .map { v -> tuple( v.target_uniprot_id, v.target_chembl_id ) }
      .unique()
      .set { id_to_uniprot_to_chembl }
   
   ( params.test ? id_to_uniprot_to_chembl.take(100) : id_to_uniprot_to_chembl )
      .map { it[0] }
      .unique()
      .toSortedList()
      .flatten()
      .buffer( size: Math.min( params.batch_size, 100 ), remainder: true )
      | fetch_fastas_from_uniprot_ids
   
   fetch_fastas_from_uniprot_ids.out
      .collectFile( 
         name: "canonical-targets.fasta", 
         storeDir: "${params.outputs}/target-sequences",
      )
      .set { uniprot_fastas }

   sample_rows
      .map { v -> tuple( v.organism_id.toString(), v.organism_id ) }
      .unique()
      | fetch_fastas_from_organism_id  // Organism ID, FASTAs gz
      | make_diamond_db

   make_diamond_db.out
      .combine( uniprot_fastas )
      | diamond_blastp

   merge_tables(
      diamond_blastp.out.data
         .combine( 
            fetch_chembl_targets.out
               .map { v -> v[-1] }, 
         ),
      Channel.value( "inner" ),
      Channel.value( false ),
      Channel.value( false ),

   )

   merge_tables.out
      .combine( fetch_gnomad_constraints.out )
      .set { named_orthologs }

   merge_tables2(
      named_orthologs,
      Channel.value( "left" ),
      Channel.value( "targets/tables" ),
      Channel.value( "target_list.tsv" ),
   )
      | set { target_lists }
   
   fetch_species_gene_names(
      target_lists.combine( fetch_target_taxonomy.out ),
      Channel.value( "ortholog_uniprot_id" ),
   )
   stack_tables4(
      fetch_species_gene_names.out
         .map { tuple( "_all", it[1] ) }
         .groupTuple( by: 0 ),
      Channel.value( "targets" ),
      Channel.value( "tsv" ),
   ) | describe

   make_rbh_matrix(
      stack_tables4.out
         .combine( fetch_taxonomic_ranks.out ),
   )
   make_rbh_matrix.out.matrix
      .combine(make_rbh_matrix.out.row_data, by: 0)
      .combine(make_rbh_matrix.out.col_data, by: 0)
      .tap { rbh_matrix }
      | (
         umaps_of_rbh_matrix 
         & silhouette
      )
   factorise_nmf(
      make_rbh_matrix.out.matrix,
   )

   find_coverage_cutoff(
      stack_tables4.out
         .combine( make_rbh_matrix.out.table, by: 0 )
   )

   describe_provenance(
      stack_tables4.out
         .combine( find_coverage_cutoff.out.cutoff, by: 0 )
         .combine( fetch_taxonomic_ranks.out )
   )

   if ( params.inhibitors ) {

      filter_target_list(
         fetch_species_gene_names.out
            .combine( find_coverage_cutoff.out.cutoff.map { v -> v[1] } ),
         Channel.value( params.min_loeuf ),
         Channel.value( params.min_identity ),
      )   

      filter_target_list.out
         .splitCsv( header: true, sep: '\t', elem: 1 )
         .map { v -> tuple( v[1].target_accession.split("\\|")[1], v[0] ) }  // uniprot_id, ID
         .unique()
         .combine( id_to_uniprot_to_chembl, by: 0 )  // uniprot_id, ID, chembl_id
         .map { v -> tuple( v[1], v[-1] ) }  // ID, chembl_id
         .unique()
         .tap { org_id_to_target_chembl }
         .map { v -> v[-1] }
         .unique()
         .set { chembl_targets_conserved }
         // .groupTuple( by: 0, sort: true )
         // .map { 
         //    v -> tuple(
         //       v[0],
         //       v[1].withIndex().collect { 
         //          el, i -> Math.round(Math.floor(i / params.batch_size)) 
         //       },
         //       v[1],
         //    ) 
         // }  //  ID, [batch_i, ...], [chembl_id, ...],
         // .transpose()  //  ID, batch_i, chembl_id
         // .groupTuple( 
         //    by: [0, 1],
         //    sort: true,
         // )  //  ID, batch_i, [chembl_id, ...]
         // .filter { v -> v[-1].size() > 0 }  // filter out trivial (size-0) elements
         // .unique()
         // .map { v -> tuple( v[0], v[2] ) }
         
         // .map { v -> v[2] }
         // .transpose()
         // .unique()

      fetch_chembl_inhibitors(
         ( params.test ? chembl_targets_conserved.take(3) : chembl_targets_conserved )
            .buffer( 
               size: params.batch_size, 
               remainder: true,
         ),
         chembl_url,
         chembl_version,
         chembl_db,
         Channel.value( params.min_pchembl ),
      )

      stack_tables(
         fetch_chembl_inhibitors.out
            .map { v -> tuple( "all", v[-1] ) }
            .groupTuple( by: 0 ),
         Channel.value( false ),
         Channel.value( false ),
      )
      stack_tables.out
         .map { v -> v[-1] }
         .tap { inhibitor_table }
         .splitCsv( 
            header: true, 
            sep: '\t', 
            elem: 1, 
            by: 1, //Math.min( params.batch_size, ( params.test ? 10 : 1000 ) ),
         )
         .set { inhibitors }

      ( params.test ? inhibitors.take(1) : inhibitors )
         .map { v -> tuple( v.target_chembl_id, v.molecule_chembl_id ) }
         // .view()
         .set { inhibitors_by_target }

      fetch_pubchem_id(
         inhibitors_by_target
            .map { v -> v[-1] }
            .buffer( size: Math.min( params.batch_size, ( params.test ? 10 : 1000 ) ), remainder: true ),
            // .transpose()
            // .unique()
            // .toSortedList()
            // .flatten()
            // .buffer( size: Math.min( params.batch_size, 100 ), remainder: true ),
         chembl_url,
         chembl_version,
         chembl_db,
      )
         | fetch_vendors

      // subset_table(
      //    fetch_vendors.out.transpose(),
      //    Channel.value( "molecule_chembl_id" ),
      //    Channel.value( false ),
      //    Channel.value( false ),
      // )
      //    | set { purchasable_cmpds }

      stack_tables2(
         inhibitors_by_target
            .map { v -> tuple( v[1], v[0] ) }  // mol chembl id, target chembl id
            .combine( fetch_vendors.out.transpose(), by: 0 )  // mol chembl id, target chembl id, vendor table
            .map { v -> tuple( v[1], v[-1] ) } // target chembl id, vendor table
            .combine( 
               org_id_to_target_chembl
                  .map { v -> tuple( v[1], v[0] ) },
               by: 0,
            )  // target chembl id, vendor table, org id
            .map { v -> tuple( v[-1], v[1] ) }  // org id, vendor table
            .groupTuple( by: 0 ),
         Channel.value( false ),
         Channel.value( false ),
      )

      stack_tables2.out.view()
         .combine( inhibitor_table )
         .set { inhib_to_target }

      merge_tables3(
         inhib_to_target,
         Channel.value( "inner" ),
         Channel.value( false ),
         Channel.value( false ),
      )
      merge_tables4(
         fetch_species_gene_names.out
            .combine(
               merge_tables3.out,
               by: 0,
            ),
         Channel.value( "inner" ),
         Channel.value( "inhibitors-with-targets" ),
         Channel.value( "tsv" ),
      )

      stack_tables3(
         merge_tables3.out
         .groupTuple( by: 0 ),
         Channel.value( "inhibitors" ),
         Channel.value( "tsv" ),

      )
         | set { pubchem_ids }

   }
   

}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/
