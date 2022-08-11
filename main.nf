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


//--------------------------------------------------------------------------
// Param Checking
//--------------------------------------------------------------------------

if(params.studyIdFile) {
  accessions = fetchRunAccessions( params.studyIdFile )
}
else {
  throw new Exception("Missing params.studyIdFile")
}

//--------------------------------------------------------------------------
// Includes
//--------------------------------------------------------------------------

include { markerGeneAnalysis } from './modules/markerGeneAnalysis.nf'

//--------------------------------------------------------------------------
// Main Workflow
//--------------------------------------------------------------------------


workflow {
  
  markerGeneAnalysis(accessions)

}