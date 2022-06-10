#!/usr/bin/perl

use strict;
use Getopt::Long;

my($platform,$samplesInfoFile,$trimLeft,$trimLeftR,$truncLen,$truncLenR,$maxLen,$mergeTechReps,$trainingSetFile,$speciesAssignmentFile,$isPaired);

&GetOptions( 
            "platform=s" => \$platform,
            #"samplesInfoFile=s" => \$samplesInfoFile,
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

system("mkdir ./filtered");
system("mkdir ./errors");
system("touch ./errors/err.rds");
system("touch featureTable.rds");

system("Rscript /usr/bin/filterFastqs.R --fastqsInDir ./fastqsInDir --fastqsOutDir ./filtered --isPaired $isPaired --trimLeft $trimLeft --trimLeftR $trimLeftR --truncLen $truncLen --truncLenR $truncLenR --maxLen $maxLen --platform $platform");

system("Rscript /usr/bin/buildErrorModels.R --fastqsInDir ./filtered --errorsOutDir ./errors --errorsFileNameSuffix err.rds --isPaired $isPaired --platform $platform");


system("Rscript /usr/bin/fastqToAsv.R  --fastqsInDir ./filtered  --errorsRdsPath ./errors/err.rds --outRdsPath ./featureTable.rds --isPaired true --platform illumina --mergeTechReps true");

#my $cmd = <<EOF; 
#Rscript /usr/bin/mergeAsvsAndAssignToOtus.R
#--asvRdsInDir $workingDir/featureTable.rds
#--assignTaxonomyRefPath $trainingSetFile
#--addSpeciesRefPath $speciesAssignmentFile
#--outPath "output"
#EOF
#$cmd =~ s/\n/ /g;
#print STDERR "Running: $cmd";
#&runCmd($cmd);
