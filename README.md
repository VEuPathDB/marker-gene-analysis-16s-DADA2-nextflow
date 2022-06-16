# Nextflow Conversion of ASVTableTask.pm

In order to run the script, you will need to supply to the nextflow.config file:
1. An apikey to be used to collect the fastq files from NCBI
2. The location of your training set file (a short example file has been supplied in /data)
3. The location of your species assignment file (a short example file has been supplied in /data)

*The parameter names for these have already been specified, you only need to enter them into the empty quotes.*

### Get Started
  * Install Nextflow
    
    `curl https://get.nextflow.io | bash`
  
  * Run the script
    
    `nextflow run VEuPathDB/ASVTableTask -with-trace -c  <config_file> -r main`

Running this locally can commonly result in an OOME. This results from the derepFastq function in buildErrorModels.R not having n specified.

`drpsF[[i]] <- derepFastq(inputFilesF[[i]], verbose = verbose)`

The default n is 1 million reads. If an OOME occurs, a lower number can be entered. 

Ex: `drpsF[[i]] <- derepFastq(inputFilesF[[i]],n = 1e+03, verbose = verbose)`

*This would then involve a rebuilding of the docker container. If you do not wish to do this, I will supply in the container a file named buildErrorsMemSafe.R. In the main.nf, this could be substituted for buildErrorModels.R. It is the same tool, but has n set to 1000.*

**Might add functionality to specify 'n' in the nextflow.config file**

Documentation for derepFastq and other dada2 functionality can be found [here](https://rdrr.io/github/benjjneb/dada2/man/derepFastq.html).
