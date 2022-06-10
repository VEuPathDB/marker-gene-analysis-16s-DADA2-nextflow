nextflow.enable.dsl=1

process test {
  input:
  set genomeName, path(genomeReads) from Channel.fromSRA("$params.sraStudyId", apiKey: "$params.apiKey")
  output:
  path 'featureTable.rds' into output_ch
  """
  mkdir ./fastqsInDir
  gzip -d --force ${genomeReads[0]} 
  gzip -d --force ${genomeReads[1]} 
  mv *.fastq* ./fastqsInDir
  perl /usr/bin/newASVTask.pl --platform $params.platform --isPaired $params.isPaired --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --truncLenR $params.truncLenR  --maxLen $params.maxLen --mergeTechReps $params.mergeTechReps --trainingSetFile $params.trainingSetFile --speciesAssignmentFile $params.speciesAssignmentFile  
  """
}

results = output_ch.collectFile(storeDir: params.outputDir)