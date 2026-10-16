nextflow.enable.dsl=2


def fetchRunAccessions( tsv ) {
    def splitter = new nextflow.splitter.CsvSplitter().options( header:true, sep:'\t' )
    List<Map> rows = splitter.target( file( tsv ) ).list()
    return rows.collect { row -> row['run_accession'] }
}


//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { markerGeneAnalysis } from './modules/markerGeneAnalysis.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------


workflow {
  
  if(params.studyIdFile) {
    accessions = fetchRunAccessions( params.studyIdFile )
  }
  else {
    error "Missing params.studyIdFile"
  }

  markerGeneAnalysis(accessions)

}