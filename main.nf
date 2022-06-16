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

process ASVPaired {
  
  publishDir params.outputDir, mode: 'copy'  
  input:
  
  tuple val(genomeName), path(genomeReads) 
  
  output:
  path '*_*' optional true
    
  """
  gzip -d --force ${genomeReads[0]} 
  gzip -d --force ${genomeReads[1]} 
  
  mkdir filtered errors feature
  touch ./errors/err.rds
  touch ./feature/featureTable.rds
  
  Rscript /usr/bin/filterFastqs.R --fastqsInDir . --fastqsOutDir ./filtered --isPaired $params.isPaired --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --truncLenR $params.truncLenR --maxLen $params.maxLen --platform $params.platform

  rm *.fastq
  
  Rscript /usr/bin/buildErrorsN.R --fastqsInDir ./filtered --errorsOutDir ./errors --errorsFileNameSuffix err.rds --isPaired $params.isPaired --platform $params.platform --nValue $params.nValue
  
  Rscript /usr/bin/fastqToAsv.R  --fastqsInDir ./filtered  --errorsRdsPath ./errors/err.rds --outRdsPath ./feature/featureTable.rds --isPaired $params.isPaired --platform $params.platform --mergeTechReps $params.mergeTechReps

  Rscript /usr/bin/mergeAsvsAndAssignToOtus.R --asvRdsInDir ./feature  --assignTaxonomyRefPath $params.trainingSet --addSpeciesRefPath $params.speciesAssignment --outPath ./"$genomeName"_$params.outputName
  """
}

process ASVSingle {
  errorStrategy 'ignore'
  publishDir params.outputDir, mode: 'copy'  
  input:
  
  tuple val(genomeName), path(genomeReads) 

  output:
  path '*_output'
  path '*_output.bootstraps'
  path '*_output.full'
  
  """
  gzip -d --force $genomeReads 
    
  mkdir filtered errors feature
  touch ./errors/err.rds
  touch ./feature/featureTable.rds
  
  Rscript /usr/bin/filterFastqs.R --fastqsInDir . --fastqsOutDir ./filtered --isPaired $params.isPaired --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --truncLenR $params.truncLenR --maxLen $params.maxLen --platform $params.platform

  rm *.fastq
  
  Rscript /usr/bin/buildErrorsMemSafe.R --fastqsInDir ./filtered --errorsOutDir ./errors --errorsFileNameSuffix err.rds --isPaired $params.isPaired --platform $params.platform
  
  Rscript /usr/bin/fastqToAsv.R  --fastqsInDir ./filtered  --errorsRdsPath ./errors/err.rds --outRdsPath ./feature/featureTable.rds --isPaired $params.isPaired --platform $params.platform --mergeTechReps $params.mergeTechReps

  Rscript /usr/bin/mergeAsvsAndAssignToOtus.R --asvRdsInDir ./feature  --assignTaxonomyRefPath $params.trainingSet --addSpeciesRefPath $params.speciesAssignment --outPath ./"$genomeName"_output
  """
}


workflow {

    accessions = fetchRunAccessions( params.sraStudyIdFile )
    if(params.isPaired == "false") {
        channel.fromSRA( accessions, apiKey: params.apiKey, protocol: "http" ) | ASVSingle
    }
    if(params.isPaired == "true"){
        channel.fromSRA( accessions, apiKey: params.apiKey, protocol: "http" ) | ASVPaired
    }
}