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
            nextflow run scbirlab/nf-reclaim --organism_id <taxon ID>
         Config/sample sheet usage:
            nextflow run scbirlab/nf-reclaim -c <config-file>

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
      gNOMAD version       : ${params.gnomad_version}
   parameters
      min. pChEMBL         : ${params.min_pchembl}
      min. LOEUF           : ${params.min_loeuf}
      min. ID              : ${params.min_identity}
      min. coverage        : ${params.min_coverage}
      batch size           : ${params.batch_size}
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
   fetch_pubchem_id;
   // fetch_included_chembl_taxon;
   // filter_targets_by_taxon;
   fetch_target_taxonomy;
   fetch_chembl_tox;
} from './modules/chembl.nf'
include { 
   make_diamond_db;
   diamond_blastp;
} from './modules/diamond.nf'
include { 
   fetch_gnomad_constraints;
} from './modules/gnomad.nf'
include { 
   fetch_pubchem_vendors;
} from './modules/pubchem.nf'
include {
   fetch_rhea_database;
   match_uniprot_to_reactants;
} from './modules/rhea.nf'
include { 
   fetch_fastas_from_organism_id;
   fetch_fastas_from_uniprot_ids;
   fetch_species_gene_names;
} from './modules/uniprot.nf'
include { 
   concat_files;
   stack_tables;
   stack_tables as stack_tables2;
   stack_tables as stack_tables3;
   merge_tables;
   filter_target_list;
   merge_tox_gnomad;
   merge_tables as merge_tables2;
   merge_tables as merge_tables3;
   merge_tables as merge_tables4;
   merge_tables as merge_tables5;
} from './modules/utils.nf'

workflow {

   // Channel.of( params.rhea_url ).set { rhea_url }
   chembl_status(
      Channel.of( workflow.start )
   )
   chembl_status.out.version.first().set { chembl_version }


   if ( params.sample_sheet ) {

      Channel.fromPath( 
         params.sample_sheet,
         checkIfExists: true, 
      )
         .splitCsv( header: true )
         .set { sample_rows }

   }

   else {

      Channel.of( [
         organism_id: params.organism_id,
      ] )
      .set { sample_rows }

   }

   // fetch_included_chembl_taxon(
   //    sample_rows.map { tuple( it.organism_id.toString(), it.organism_id ) }.unique(),
   //    Channel.value( params.chembl_url ),
   //    chembl_version,
   //    Channel.value( params.default_kingdom )
   // )
   fetch_target_taxonomy(
      Channel.of( params.chembl_url ),
      chembl_version,
   )
   fetch_chembl_targets(
      Channel.of( params.chembl_url ),
      chembl_version,
   )  
   fetch_gnomad_constraints(
      Channel.of( params.gnomad_version ),
   )

   if ( !params.test ) {

      fetch_chembl_tox(
         Channel.of( params.chembl_url ),
         Channel.value( params.tox_cell_lines ),
         chembl_version,
      )

      merge_tox_gnomad(
         fetch_chembl_targets.out.map { it[-1] },
         fetch_chembl_tox.out.main,
         fetch_gnomad_constraints.out,
         Channel.value( "inner" ),
      )
      
   }

   // fetch_included_chembl_taxon.out
   //    .combine( fetch_chembl_targets.out )
   //    | filter_targets_by_taxon
   // filter_targets_by_taxon.out
   fetch_chembl_targets.out
      .map { it[-1] }
      .splitCsv( 
         header: true, 
         elem: 0, //1, 
         sep: '\t',
      )
      // .map { tuple( it[0].toString(), it[1].target_uniprot_id, it[1].target_chembl_id ) }
      .map { tuple( it.target_uniprot_id, it.target_chembl_id ) }
      .unique()
      .set { id_to_uniprot_to_chembl }
   
   ( params.test ? id_to_uniprot_to_chembl.take(100) : id_to_uniprot_to_chembl )
      .map { it[0] }
      .unique()
      .toSortedList()
      .flatten()
      .buffer( size: params.batch_size, remainder: true )
      | fetch_fastas_from_uniprot_ids
   
   fetch_fastas_from_uniprot_ids.out
      // .take(10)
      // .splitFasta( record: [id: true, text: true] )
      // .take(10)
      // .map { tuple( it.id, it.id.split("\\|")[1], it.text ) }
      // .collectFile { [ "${it[1]}.fasta", it[2] ] }
      // .map { tuple( it.getSimpleName(), it ) }
      // .combine( 
      //    // id_to_uniprot_to_chembl.map { it[1..0] }, 
      //    id_to_uniprot_to_chembl, //.map { it[1] }, 
      //    // by: 0,
      // )
      // .map { tuple( it[-1], it[1] ) }
      .collectFile( 
         name: "canonical-targets.fasta", 
         // sort: true, 
         storeDir: "${params.outputs}/target-sequences",
      )
      // .map { tuple( it.getSimpleName(), it ) }
      .view()
      .set { uniprot_fastas }

   sample_rows
      .map { tuple( it.organism_id.toString(), it.organism_id ) }
      .unique()
      | fetch_fastas_from_organism_id  // Organism ID, FASTAs gz
      | make_diamond_db

   make_diamond_db.out
      .combine( uniprot_fastas )
      | diamond_blastp

   merge_tables(
      // filter_targets_by_taxon.out
      // .map { [ it[0].toString() ] + it[1..-1] }
      diamond_blastp.out.data
         // .map { [ it[0].toString() ] + it[1..-1] }
         .combine( 
            fetch_chembl_targets.out.map { it[-1] }, 
            // by: 0,
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
      Channel.value( "targets" ),
      Channel.value( "target_list.tsv" ),
   )
      | set { target_lists }
   
   fetch_species_gene_names(
      target_lists.combine( fetch_target_taxonomy.out ),
      Channel.value( "ortholog_uniprot_id" ),
   )

   filter_target_list(
      fetch_species_gene_names.out,
      Channel.value( params.min_loeuf ),
      Channel.value( params.min_identity ),
      Channel.value( params.min_coverage ),
   )

   filter_target_list.out
      .splitCsv( header: true, sep: '\t', elem: 1 )
      .map { tuple( it[1].target_accession.split("\\|")[1], it[0]  ) }
      .unique()
      .combine( id_to_uniprot_to_chembl, by: 0 )
      .map { tuple( it[1], it[-1] ) }
      .unique()
      .set { chembl_targets_conserved }
   fetch_chembl_inhibitors(
      ( params.test ? chembl_targets_conserved.take(10) : chembl_targets_conserved ),
      Channel.value( params.chembl_url ),
      chembl_version,
      Channel.value( params.min_pchembl ),
   )

   stack_tables(
      fetch_chembl_inhibitors.out
         .groupTuple( by: 0 ),
      Channel.value( false ),
      Channel.value( false ),
   )
   // fetch_chembl_inhibitors.out
   //    .collectFile( 
   //       { [ "${it[0]}-inhibitors.tsv", it[1] ] }, 
   //       keepHeader: true, 
   //       skip: 1, 
   //       storeDir: "${params.outputs}/inhibitors",
   //    )
      // | map { tuple( it.getSimpleName().split(".inhibitors")[0], it ) }
   stack_tables.out
      .tap { inhibitor_table }
      .splitCsv( header: true, sep: '\t', elem: 1 )
      .set { inhibitors }

   ( params.test ? inhibitors.take(10) : inhibitors )
      .map { tuple( it[0].toString(), it[1].target_chembl_id, it[1].molecule_chembl_id ) }
      .set { inhibitors_by_target }

   fetch_pubchem_id(
      inhibitors_by_target
         .map { it[2] }
         .unique(),
      Channel.value( params.chembl_url ),
      chembl_version,
   )

   stack_tables2(
      inhibitors_by_target
         .map { tuple( it[2], it[0] ) }
         .combine( fetch_pubchem_id.out, by: 0 )
         .map { tuple( it[1], it[-1] ) }
         .groupTuple( by: 0 ),
      Channel.value( false ),
      Channel.value( false ),

   )
   // inhibitors_by_target
   //    .map { tuple( it[2], it[0] ) }
   //    .combine( fetch_pubchem_id.out, by: 0 )
   //    .map { tuple( it[1], it[-1] ) }
   //    .collectFile( 
   //       { [ "${it[0]}-pubchem.tsv", it[1] ] }, 
   //       keepHeader: true, 
   //       skip: 1, 
   //       storeDir: "${params.outputs}/inhibitors",
   //    )
   //    .map { tuple( it.getSimpleName().split("-pubchem")[0], it ) }
   stack_tables2.out
      .combine( inhibitor_table, by: 0 )
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
  
      // | stack_tables
      // .collectFile( 
      //    { [ "${it[0]}.conserved.tsv", it[1] ] }, 
      //    keepHeader: true, 
      //    skip: 1, 
      //    storeDir: "${params.outputs}/inhibitors",
      // )
      // .set { pubchem_ids }

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