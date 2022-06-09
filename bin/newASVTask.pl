#!/usr/bin/perl

use Utils;
use strict;
use Getopt::Long;

my($sraStudyId,$sraSampleAndRunIdsPath,$dataDir,$workingDir,$platform,$samplesInfoFile,$trimLeft,$trimLeftR,$truncLen,$truncLenR,$maxLen,$mergeTechReps,$trainingSetFile,$speciesAssignmentFile,$isPaired);
$workingDir = ".";

###############Collecting Inputs#################################
&GetOptions( 
            "sraStudyId=s"=> \$sraStudyId,
            "sraSampleAndRunIdsPath=s" => \$sraSampleAndRunIdsPath,
            "workingDir=s" => \$workingDir,
            "platform=s" => \$platform,
            "samplesInfoFile=s" => \$samplesInfoFile,
            "trimLeft=s" => \$trimLeft,
            "trimLeftR=s" => \$trimLeftR,
            "truncLen=s" => \$truncLen,
            "truncLenR=s" => \$truncLenR,
            "maxLen=s" => \$maxLen,
            "mergeTechReps" => \$mergeTechReps,
            "trainingSetFile=s" => \$trainingSetFile,
            "speciesAssignmentFile=s" => \$speciesAssignmentFile,
            "isPaired=s" => \$isPaired,
            );

$dataDir = ".";
my $fastqsInDir; 
###################Checking Inputs###############################
die "Must provide either --sraStudyId or --dataDir or --sraSampleAndRunIdsPath"
  unless ($dataDir ne "none" || $sraStudyId ne "none" || $sraSampleAndRunIdsPath ne "none");
if ($platform eq "illumina") {
  if ($isPaired ne 'true' && $isPaired ne 'false') {
    die "argument to 'isPaired' must be 'true' or 'false'";
  }
} elsif ($platform eq '454') {
   $isPaired = 'false'; 
} else {
    die "currently only 'illumina' and '454' are supported platforms.";
} 

if ($maxLen eq "none" and $platform eq '454') {
    die "must provide maximum read length as parameter 'maxLen' for 454 data.";
} 

#####################Generating fastqsInDir##############################
my $samplesDir = $dataDir;
# If dataDir was specified, will unzip and rename all files for fastq
if ($samplesDir ne "none"){
  my @zipFiles = glob("$samplesDir/*.gz");
  foreach my $zipFile (@zipFiles) {
    print STDERR "unzipping $zipFile\n";
    `gunzip $zipFile`;
    }
  
  if (glob("$samplesDir/*.fq")) {
    print STDERR "changing samples file ext from .fq to .fastq\n";
    `rename .fq .fastq $samplesDir/*.fq`;
  }
  die "$samplesDir contains no fastqs"
    unless glob("$samplesDir/*fastq");
     
  if (glob("$samplesDir/*_R1_001.fastq")){
    `rename _R1_001.fastq _1.fastq $samplesDir/*_R1_001.fastq`;
  }
  if (glob("$samplesDir/*_R2_001.fastq")){
    `rename _R2_001.fastq _2.fastq $samplesDir/*_R2_001.fastq`;
  }

  if($isPaired eq 'true'){
     die "No forward fastqs" unless glob("$samplesDir/*_1.fastq");
     die "No reverse fastqs" unless glob("$samplesDir/*_2.fastq");
  }
  $fastqsInDir = $samplesDir;
  # If sraStudyID was specified
}elsif ($sraStudyId) {
  `mkdir $workingDir/fastqsFromSra`;
  &runCmd("getFastqFromSra.pl --workingDir $workingDir/fastqsFromSra --studyId '$sraStudyId' --pairs $isPaired"); 
  $fastqsInDir = "$workingDir/fastqsFromSra";
}elsif($sraSampleAndRunIdsPath){
  `mkdir $workingDir/fastqsFromSra`;
  &runCmd("getFastqFromSra.pl --workingDir $workingDir/fastqsFromSra --sampleAndRunIdsPath '$sraSampleAndRunIdsPath' --pairs $isPaired");
  $fastqsInDir = "$workingDir/fastqsFromSra";
}

##################Running filterFastqs.R on fastqsInDir######################
my $cmd = <<"EOF";
Rscript /usr/bin/filterFastqs.R 
--fastqsInDir $fastqsInDir
--fastqsOutDir $workingDir/filtered
--isPaired $isPaired 
--maxLen $maxLen 
--platform $platform 
EOF
    $cmd =~ s/\n/ /g;
    $cmd .= " --samplesInfo $samplesInfoFile" if $samplesInfoFile ne "false";
    $cmd .= " --trimLeft $trimLeft" if $trimLeft ne "false";
    $cmd .= " --trimLeftR $trimLeftR" if $trimLeftR ne "false";
    $cmd .= " --truncLen $truncLen" if $truncLen ne "false";
    $cmd .= " --truncLenR $truncLenR" if $truncLenR ne "false";

print STDERR "Running: $cmd";
&runCmd($cmd);
die "Failed: no fastqs in output filtering dir $workingDir/filtered"
  unless glob("$workingDir/filtered/*.fastq");

######################Running buildErrorModels.R#############################
my $cmd = <<"EOF";
Rscript /usr/bin/buildErrorModels.R 
--fastqsInDir $workingDir/filtered
--errorsOutDir $workingDir/errors 
--errorsFileNameSuffix err.rds
--isPaired $isPaired 
--platform $platform 
EOF
$cmd =~ s/\n/ /g;
$cmd .= " --samplesInfo $samplesInfoFile" if $samplesInfoFile;
print STDERR "Running: $cmd";
&runCmd($cmd);

  # These naughty files triggered this bug, fixed in Jan 2020: https://github.com/Bioconductor/ShortRead/pull/2
  # If you're reading this in year 2022 verify Bioconductor has released sometime in 2020 or 2021, and that you have the good ShortRead locally
  # Then delete this code and the message
  if ($sraStudyId eq 'ERP005806'){
    unlink("$workingDir/filtered/ERS454101.ERR503975_filt.fastq");
    unlink("$workingDir/filtered/ERS453453.ERR503327_filt.fastq");
  }

######################Running fastqToAsv.R###############################
my $cmd = <<EOF;
Rscript /usr/bin/fastqToAsv.R
--fastqsInDir $workingDir/filtered
--errorsRdsPath $workingDir/errors
--outRdsPath $workingDir/featureTable.rds
--isPaired $isPaired
--platform $platform
--mergeTechReps $mergeTechReps
EOF
$cmd =~ s/\n/ /g;
print STDERR "Running: $cmd";
&runCmd($cmd);

my $cmd = <<EOF; 
Rscript /usr/bin/mergeAsvsAndAssignToOtus.R
--asvRdsInDir $workingDir/featureTable.rds
--assignTaxonomyRefPath $trainingSetFile
--addSpeciesRefPath $speciesAssignmentFile
--outPath "output"
EOF
$cmd =~ s/\n/ /g;
print STDERR "Running: $cmd";
&runCmd($cmd);







