#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process downloadFiles {
  container ='veupathdb/bowtiemapping:1.0.0'  

  input:
    val id

  output:
    tuple val(id), path("${id}**.fastq")

  script:
    template 'downloadFiles.bash'
}


process filterFastqs {
  container = 'veupathdb/markergeneanalysis16sdada2:1.0.0'
  
  input:
    tuple val(genomeName), path(fastqs)

  output:
    tuple val(genomeName), path('filtered/*.fast*')

  script:
    template 'filterFastqs.bash'
}


process buildErrors {
  container = 'veupathdb/markergeneanalysis16sdada2:1.0.0'

  input:
    tuple val(genomeName), path(fastasfiltered)

  output:
    tuple val(genomeName), path('err.rds'), path('filtered/*.fast*')

  script:
    template 'buildErrors.bash'
}


process fastqToAsv {
  container = 'veupathdb/markergeneanalysis16sdada2:1.0.0'

  input:
    tuple val(genomeName), path(errorFile), path(fastqsFiltered)

  output:
    tuple val(genomeName), path('featureTable.rds')

  script:
    template 'fastqToAsv.bash'
}


process mergeAsvsAndAssignToOtus {
  container = 'veupathdb/markergeneanalysis16sdada2:1.0.0'

  publishDir params.outputDir, mode: 'copy'

  input:
    tuple val(genomeName), path('featureTable.rds')

  output:
    path '*_output'
    path '*_output.bootstraps'
    path '*_output.full'

  script:
    template 'mergeAsvsAndAssignToOtus.bash'
}


workflow markerGeneAnalysis {
  take:
    accessions

  main:
    ids = Channel.fromList( accessions )
    downloadFiles( ids ) \
      | filterFastqs \
      | buildErrors \
      | fastqToAsv \
      | mergeAsvsAndAssignToOtus
}
