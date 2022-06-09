nextflow.enable.dsl=1

fileInput_qch = Channel.fromPath("$params.dataDir/*")

process ASVTableTask {
    input:
    path inputFile from fileInput_qch 
    output:
    path 'output' into results_qch
    """
    newASVTask.pl --sraStudyId $params.sraStudyId --sraSampleAndRunIdsPath $params.sraSampleAndRunIdsPath --platform $params.platform --isPaired $params.isPaired --samplesInfoFile $params.samplesInfoFile --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --truncLenR $params.truncLenR  --maxLen $params.maxLen --mergeTechReps $params.mergeTechReps --trainingSetFile $params.trainingSetFile --speciesAssignmentFile $params.speciesAssignmentFile 
    """
}

//results = results_qch.collectFile(storeDir: params.outputDir, name: params.resultsFile)

