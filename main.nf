nextflow.enable.dsl=1

process test {
  input:
  
  set genomeName, path(genomeReads) from Channel.fromSRA("$params.sraStudyId", apiKey: "$params.apiKey")
  
  output:
  path 'featureTable.rds' into output_ch
  
  """
  gzip -d --force ${genomeReads[0]} 
  gzip -d --force ${genomeReads[1]} 

  mkdir ./filtered
  mkdir ./errors
  touch ./errors/err.rds
  touch featureTable.rds

  Rscript /usr/bin/filterFastqs.R --fastqsInDir . --fastqsOutDir ./filtered --isPaired $params.isPaired --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --truncLenR $params.truncLenR --maxLen $params.maxLen --platform $params.platform

  Rscript /usr/bin/buildErrorModels.R --fastqsInDir ./filtered --errorsOutDir ./errors --errorsFileNameSuffix err.rds --isPaired $params.isPaired --platform $params.platform

  Rscript /usr/bin/fastqToAsv.R  --fastqsInDir ./filtered  --errorsRdsPath ./errors/err.rds --outRdsPath ./featureTable.rds --isPaired $params.isPaired --platform $params.platform --mergeTechReps $params.mergeTechReps
  """
}

results = output_ch.collectFile(storeDir: params.outputDir)