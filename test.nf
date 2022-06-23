import nextflow.splitter.CsvSplitter

nextflow.enable.dsl=2

def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

process prepASVPaired {
  input:
  tuple val(genomeName), path(genomeReads) 
  output:
  val $genomeName
  path '*.fasta'
  """
  gzip -d --force ${genomeReads[0]} 
  gzip -d --force ${genomeReads[1]} 
  """
}

process prepASVSingle {
  input:
  tuple val(genomeName), path(genomeReads) 
  output:
  tuple val(genomeName), path('*.fasta')
  """
  gzip -d --force $genomeReads 
  """  
}

process filterFastqs {
  input:
  val genomeName
  path fastqs
  output:
  tuple val(genomeName), path('./filtered/*')
  """
  Rscript /usr/bin/filterFastqs.R --fastqsInDir . --fastqsOutDir ./filtered --isPaired $params.isPaired --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --truncLenR $params.truncLenR --maxLen $params.maxLen --platform $params.platform
  rm *.fastq
  """
}

process buildErrors {
  input:
  tuple val(genomeName), path(fastasfiltered)
  output:
  tuple val(genomeName), path('err.rds'), path('*.fastq')
  """
  Rscript /usr/bin/buildErrorsN.R --fastqsInDir . --errorsOutDir . --errorsFileNameSuffix err.rds --isPaired $params.isPaired --platform $params.platform --nValue $params.nValue
  """
}

process fastqToAsv {
  input:
  tuple val(genomeName), path('err.rds'), path(fastqsFiltered)
  output:
  tuple val(genomeName), path('featureTable.rds')
  """
  Rscript /usr/bin/fastqToAsv.R  --fastqsInDir .  --errorsRdsPath ./err.rds --outRdsPath ./featureTable.rds --isPaired $params.isPaired --platform $params.platform --mergeTechReps $params.mergeTechReps
  """
}

process mergeAsvsAndAssignToOtus {
   errorStrategy 'ignore'
   publishDir params.outputDir, mode: 'copy'
   input:
   tuple val(genomeName), path('featureTable.rds')
   output:
   path '*_output'
   path '*_output.bootstraps'
   path '*_output.full'
   """
   Rscript /usr/bin/mergeAsvsAndAssignToOtus.R --asvRdsInDir . --assignTaxonomyRefPath $params.trainingSet --addSpeciesRefPath $params.speciesAssignment --outPath ./"$genomeName"_output
   """
}

workflow {
    accessions = fetchRunAccessions( params.studyIdFile )
    if(params.isPaired == "false") {
        channel.fromSRA( accessions, apiKey: params.apiKey, protocol: "http" ) | prepASVSingle | filterFastqs | buildErrors | fastqToAsv | mergeAsvsAndAssignToOtus
    }
    if(params.isPaired == "true"){
        channel.fromSRA( accessions, apiKey: params.apiKey, protocol: "http" ) | prepASVPaired | filterFastqs | buildErrors | fastqToAsv | mergeAsvsAndAssignToOtus
    } 
}