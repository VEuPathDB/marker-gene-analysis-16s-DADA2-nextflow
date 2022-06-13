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

process ASVTableTask {
  debug true
  input:
  
  tuple val(genomeName), path(genomeReads) 
  
  output:
  path 'featureTable.rds' 
  
  """
  echo "Unzipping FastQs"
  gzip -d --force ${genomeReads[0]} 
  gzip -d --force ${genomeReads[1]} 
  
  echo "Making Directories"
  mkdir ./filtered
  mkdir ./errors
  touch ./errors/err.rds
  touch featureTable.rds
  
  echo "Running Filter"
  Rscript /usr/bin/filterFastqs.R --fastqsInDir . --fastqsOutDir ./filtered --isPaired $params.isPaired --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --truncLenR $params.truncLenR --maxLen $params.maxLen --platform $params.platform
  
  echo "Building Error File"
  Rscript /usr/bin/buildErrorModels.R --fastqsInDir ./filtered --errorsOutDir ./errors --errorsFileNameSuffix err.rds --isPaired $params.isPaired --platform $params.platform
  
  echo "Running FastqtoASV"
  Rscript /usr/bin/fastqToAsv.R  --fastqsInDir ./filtered  --errorsRdsPath ./errors/err.rds --outRdsPath ./featureTable.rds --isPaired $params.isPaired --platform $params.platform --mergeTechReps $params.mergeTechReps
  """
}


workflow {

    accessions = fetchRunAccessions( params.sraStudyIdFile )

    fastq_qch = channel
        .fromSRA( accessions, apiKey: params.apiKey )
        .view()
	
    results = ASVTableTask(fastq_qch)

    results = results.collectFile(storeDir: params.outputDir)
}