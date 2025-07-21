#!/usr/bin/env nextflow

import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

// Define a formatter
DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss")

// Get the current date and time
LocalDateTime now = LocalDateTime.now()

// Format the current date and time
String formattedDateTime = now.format(formatter)

// Print the formatted date and time
println "Current date and time: $formattedDateTime"

nextflow.enable.dsl = 2

include { hash_seqs }                                           from './modules/hash_seqs.nf'
include { seq_qc }                                              from './modules/blast.nf'
include { blastn }                                              from './modules/blast.nf'
include { blastn_ncbi }                                         from './modules/blast.nf'
include { taxonkit_annotation as taxonkit_annotation_local }    from './modules/blast.nf'
include { taxonkit_annotation as taxonkit_annotation_ncbi }     from './modules/blast.nf'
include { filter_by_regex as filter_by_regex_local }            from './modules/blast.nf'
include { filter_by_regex as filter_by_regex_ncbi }             from './modules/blast.nf'
include { filter_best_bitscore as filter_best_bitscore_local }  from './modules/blast.nf'
include { filter_best_bitscore as filter_best_bitscore_ncbi }   from './modules/blast.nf'
include { build_report }                                        from './modules/blast.nf'
include { collect_provenance }                                  from './modules/provenance.nf'
include { pipeline_provenance }                                 from './modules/provenance.nf'


workflow {

  ch_pipeline_metadata = Channel.value([
    workflow.sessionId,
    workflow.runName,
    workflow.manifest.name,
    workflow.manifest.version,
    workflow.start,
  ])

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fasta = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['FILE']] }
  } else {
    ch_fasta = Channel.fromPath(params.fasta_search_path)
  }

  if (params.databases != 'NO_FILE') {
    ch_db = Channel.fromPath(params.databases).splitCsv(header: true).map{ it -> [it['ID'], it['DBNAME'], it['PATH']] }
  } else {
    ch_db = Channel.of()
  }

  ch_ncbi_db = Channel.fromPath(params.ncbi_db)

  ch_seqs = ch_fasta.splitFasta(record: [id: true, seqString: true])

  main:

    hash_seqs(ch_seqs)


    seq_qc(ch_seqs)
    ch_blast = blastn(ch_seqs.combine(ch_db)).blast_report
    ch_blast = taxonkit_annotation_local(ch_blast).blast_report

    ch_blast_ncbi = blastn_ncbi(ch_seqs.combine(ch_ncbi_db)).blast_report
    ch_blast_ncbi = taxonkit_annotation_ncbi(ch_blast_ncbi).blast_report

    if (params.filter_regexes != 'NO_FILE') {
      ch_regexes = Channel.fromPath(params.filter_regexes)
      ch_blast = filter_by_regex_local(ch_blast.combine(ch_regexes)).blast_filtered
      ch_blast_ncbi = filter_by_regex_ncbi(ch_blast_ncbi.combine(ch_regexes)).blast_filtered
    }

    ch_blast_collect = ch_blast.collectFile(it -> it[2], name: "collected_blast.csv", storeDir: params.outdir, keepHeader: true, skip: 1)
    
    ch_blast_ncbi_collect = ch_blast_ncbi.collectFile(it -> it[2], name: "collected_blast_ncbi.csv", storeDir: params.outdir, keepHeader: true, skip: 1)

    filter_best_bitscore_local(ch_blast)

    filter_best_bitscore_ncbi(ch_blast_ncbi)
    
    filter_best_bitscore_local.out.blast_best_bitscore_csv.collectFile(it -> it[1], name: "collected_blast_best_bitscore.csv", storeDir: params.outdir, keepHeader: true, skip: 1)

    filter_best_bitscore_ncbi.out.blast_best_bitscore_csv.collectFile(it -> it[1], name: "collected_blast_ncbi_best_bitscore.csv", storeDir: params.outdir, keepHeader: true, skip: 1)


    build_report(ch_blast_collect, ch_blast_ncbi_collect, Channel.fromPath(params.databases))

    // Build pipeline provenance 
    ch_pipeline_provenance = pipeline_provenance(ch_pipeline_metadata, build_report.out.provenance)
    
    //Pool Provenance data
    ch_provenance = hash_seqs.out.provenance
    ch_provenance = ch_provenance.join(seq_qc.out.provenance).map{ it -> [it[0], [it[1]] << it[2]] }
    ch_provenance = ch_provenance.join(blastn.out.provenance.groupTuple()).map{ it -> [it[0], (it[1] + it[2]).flatten() ] } 
    ch_provenance = ch_provenance.join(blastn_ncbi.out.provenance.groupTuple()).map{ it -> [it[0], (it[1] + it[2]).flatten() ] } 
    //ch_provenance = ch_provenance.join(filter_best_bitscore.out.provenance.groupTuple()).map{ it -> [it[0], (it[1] + it[2]).flatten()] }
    ch_provenance = ch_provenance.join(seq_qc.out.provenance.map{it -> it[0]}.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }
    collect_provenance(ch_provenance)
}
