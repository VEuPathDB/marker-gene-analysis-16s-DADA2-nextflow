nextflow.enable.dsl=1

process ASVTableTask {
    output:
    path '' into fastqInDir_vch
    """
    newASVTask.pl --sraStudyId $params.sraStudyId --sraSampleAndRunIdsPath $params.sraSampleAndRunIdsPath --dataDir $params.dataDir --platform $params.platform --samplesInfoFile $params.samplesInfoFile --trimLeft $params.trimLeft --trimLeftR $params.trimLeftR --truncLen $params.truncLen --maxLen $params.maxLen --mergeTechReps $params.mergeTechReps --trainingSetFile $params.trainingSetFile --speciesAssignmentFile $params.speciesAssignmentFile --resultFile $params.resultFile --isPaired $params.isPaired --samplesInfo $params.samplesInfo
    """
}

samResults = sam_qch.collectFile(storeDir: params.outputDir, name: params.samFile)

