package Sra;
#use autodie qw/system/; # requires 'IPC::System::Simple';
use LWP::Simple; 
use XML::Simple; 
use Data::Dumper;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(getRunIdsFromSraSampleIds getFastqForSraRunId getFastqForSampleIds getCsForSampleIds getRunIdsFromSraStudyId getFastqForStudyId runEfetch);
use strict;
use warnings;

sub getFastqForStudyId {
  my ($studyId, @args) = @_;
  print STDERR "fetching runs for study: $studyId\n";
  my @runs = &getRunsForStudy($studyId);
  print STDERR sprintf("fetched %s runs for study $studyId\n", scalar @runs);
  return getFastqForStudyRuns(\@runs, @args);
}
sub getFastqForStudyRuns {
  my ($runs, $isPairedEndDerp, $dontDownload, $deflineVars) = @_;

  my %sampleCount;

  FASTQ_FOR_RUN:
  for my $run (@{$runs}) {
    my ($sampleId, $runId, $libraryLayout) = @$run;

    # Derp the per-study $isPairedEnd, and meanwhile use it for a warning if $libraryLayout disagrees with it
    if (defined $isPairedEndDerp && (($isPairedEndDerp eq 'false' and $libraryLayout eq 'PAIRED') or ($isPairedEndDerp ne 'false' and $libraryLayout eq 'SINGLE'))){
      warn "getFastqForStudyId $sampleId has \$isPairedEnd=$isPairedEndDerp but $runId has library layout $libraryLayout";
    }

    next if $dontDownload;

    if ($libraryLayout eq 'PAIRED') {
      my $isPairedEnd = 1;
      next FASTQ_FOR_RUN if -f "$sampleId.${runId}_1.fastq" and -f "$sampleId.${runId}_2.fastq";
      getFastqForSraRunId($runId, $isPairedEnd , $deflineVars);
      rename("${runId}_1.fastq", "$sampleId.${runId}_1.fastq");
      rename("${runId}_2.fastq", "$sampleId.${runId}_2.fastq");
    } elsif ($libraryLayout eq 'SINGLE') {
      my $isPairedEnd = 0;
      next FASTQ_FOR_RUN if -f "$sampleId.${runId}.fastq";
      getFastqForSraRunId($runId, $isPairedEnd , $deflineVars);
      rename("${runId}_1.fastq", "$sampleId.${runId}.fastq"); 
    } else {
      die "Unknown library layout: $libraryLayout";
    }
  }
}

# The SRA payload has 15 MB and all the fields, doesn't let you slice it, and XML::Simple is too simple for it
# Locality doesn't matter when we're not downloading fastqs so use ENA instead
sub getRunsForStudy {
  my ($study_id) = @_;
  my $url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$study_id&result=read_run&fields=sample_accession,secondary_sample_accession,run_accession,library_layout";
  my $text = get ($url);
  die "No data for URL: $url" unless $text;
  open (my $fh, "<", \$text) or die $!;
  die "Bad header: $url" unless (my $header = <$fh>) eq "sample_accession\tsecondary_sample_accession\trun_accession\tlibrary_layout\n";
  my @result;
  while(<$fh>){
     chomp;
     my ($theOtherSampleId, $sampleId, $runId, $libraryLayout) = split "\t";
     die "Bad line: $_" unless ($theOtherSampleId =~ /^SAMN\d+$/ and $sampleId =~ /^[SED]RS\d+$/ and $runId =~ /^[SED]RR\d+$/ and $libraryLayout =~ /^SINGLE|PAIRED$/);
     push @result, [$sampleId, $runId, $libraryLayout];
  }
  return @result;
}

sub getFastqForSampleIds {
  my($sids,$fileoutone,$fileouttwo,$dontdownload,$isPairedEnd) = @_;

  $fileoutone = $fileoutone ? $fileoutone : "reads_1.fastq";
  $fileouttwo = $fileouttwo ? $fileouttwo : "reads_2.fastq";
  my @rids;
  foreach my $sid (@{$sids}){
    push(@rids,&getRunIdsFromSraSampleId($sid));
  }
  my @out;
  my $readCount = 0;
  my $fid = $rids[0]->[1];
  my $singleEnd = 0;
  my $doubleEnd = 0;
  my %tsid;
  my %done;
  foreach my $a (@rids){
    $tsid{$a->[0]} = 1;
    if($done{$a->[1]}){
      print STDERR "ERROR: already retrieved '$a->[1]' .. skipping\n";
      next;
    }
    push(@out,"$a->[0]:$a->[1]:$a->[2]:".($a->[2] ? $a->[3] / $a->[2] : 'undef').":$a->[4]");
    $readCount += $a->[2];
    my $id = $a->[1];
    $done{$id} = 1;
    next if $dontdownload;
    print STDERR " paired end is $isPairedEnd\n";
    &getFastqForSraRunId($id,$isPairedEnd);
    ##if single end will have single fastq file labeled _1.fastq .. otherwise two labeled _1 and _2.fastq
    my $foundFile = 0;
    if(-e "$id\_1.fastq"){
      $foundFile++;
      if($fid ne $id){
        system("cat $id\_1.fastq >> tmpReads_1.fastq");
        unlink("$id\_1.fastq");
      }else{
        rename("$id\_1.fastq","tmpReads_1.fastq");
      }
    }
    if(-e "$id\_3.fastq"){
      $doubleEnd = 1;
      die "ERROR: this sample '$a->[0]' contains both single and double end reads\n" if $singleEnd && $doubleEnd;
      $foundFile++;
      if($fid ne $id){
        system("cat $id\_3.fastq >> tmpReads_2.fastq");
        unlink("$id\_3.fastq");
      }else{
        rename("$id\_3.fastq","tmpReads_2.fastq");
      }
      unlink("$id\_2.fastq");  ##this one is the barcode
    }elsif(-e "$id\_2.fastq"){
      $doubleEnd = 1;
      die "ERROR: this sample '$a->[0]' contains both single and double end reads\n" if $singleEnd && $doubleEnd;
      $foundFile++;
      if($fid ne $id){
        system("cat $id\_2.fastq >> tmpReads_2.fastq");
        unlink("$id\_2.fastq");
      }else{
        rename("$id\_2.fastq","tmpReads_2.fastq");
      }
    }else{
      ##is single end only
      $singleEnd = 1;
      die "ERROR: this sample '$a->[0]' contains both single and double end reads\n" if $singleEnd && $doubleEnd;
    }
  } 
  print "input: (",join(", ",@{$sids}),") ", scalar(keys%tsid), " samples, ", scalar(@rids) , " runs, $readCount spots: " , (scalar(@rids) == 0 ? "ERROR: unable to retrieve runIds\n" : "(",join(", ",@out),")\n");
  ##now mv the files to a final filename ...
#  rename("$fid.fastq","reads.fastq") if (-e "$fid.fastq");
  rename("tmpReads_1.fastq","$fileoutone") if (-e "tmpReads_1.fastq");
  rename("tmpReads_2.fastq","$fileouttwo") if (-e "tmpReads_2.fastq");
}

sub getCsForSampleIds {
  my($sids,$fileoutone,$fileouttwo,$dontdownload) = @_;
  $fileoutone = $fileoutone ? $fileoutone : "reads_1.csfasta";
  $fileouttwo = $fileouttwo ? $fileouttwo : "reads_2.csfasta";
  my @rids;
  foreach my $sid (@{$sids}){
    push(@rids,&getRunIdsFromSraSampleId($sid));
  }
  my @out;
  my $readCount = 0;
  my $fid = $rids[0]->[1];
  my $singleEnd = 0;
  my $doubleEnd = 0;
  my %tsid;
  my %done;
  foreach my $a (@rids){
    $tsid{$a->[0]} = 1;
    if($done{$a->[1]}){
      print STDERR "ERROR: already retrieved '$a->[1]' .. skipping\n";
      next;
    }
    push(@out,"$a->[0]:$a->[1]:$a->[2]:".($a->[2] ? $a->[3] / $a->[2] : 'undef').":$a->[4]");
    $readCount += $a->[2];
    my $id = $a->[1];
    $done{$id} = 1;
    next if $dontdownload;
    &getCsForSraRunId($id);
    ##if single end will have two files labeled _F3.csfasta and _F3_QV.qual .. otherwise four labeled _F3.csfasta and _F3_QV.qual, _R3.csfasta and _R3_QV.qual
    my $foundFile = 0;
    if(-e "$id\_F3.csfasta" && -e "$id\_F3_QV.qual"){
      $foundFile++;
      if($fid ne $id){
        system("cat $id\_F3.csfasta >> tmpReads_1.csfasta");
        system("cat $id\_F3_QV.qual >> tmpReads_1.csfasta.qual");
        unlink("$id\_F3.csfasta");
        unlink("$id\_F3.csfasta");
      }else{
        rename("$id\_F3.csfasta","tmpReads_1.csfasta");
        rename("$id\_F3_QV.qual","tmpReads_1.csfasta.qual");
      }
    }
    if(-e "$id\_R3.csfasta" && -e "$id\_R3_QV.qual"){
      $doubleEnd = 1;
      die "ERROR: this sample '$a->[0]' contains both single and double end reads\n" if $singleEnd && $doubleEnd;
      $foundFile++;
      if($fid ne $id){
        system("cat $id\_R3.csfasta >> tmpReads_2.csfasta");
        system("cat $id\_R3_QV.qual >> tmpReads_2.csfasta.qual");
        unlink("$id\_R3.csfasta");
        unlink("$id\_R3.csfasta");
      }else{
        rename("$id\_R3.csfasta","tmpReads_2.csfasta");
        rename("$id\_R3_QV.qual","tmpReads_2.csfasta.qual");
      }
    }else{
      ##is single end only
      $singleEnd = 1;
      die "ERROR: this sample '$a->[0]' contains both single and double end reads\n" if $singleEnd && $doubleEnd;
    }
  } 
  print "input: (",join(", ",@{$sids}),") ", scalar(keys%tsid), " samples, ", scalar(@rids) , " runs, $readCount spots: " , (scalar(@rids) == 0 ? "ERROR: unable to retrieve runIds\n" : "(",join(", ",@out),")\n");
  ##now mv the files to a final filename ...

  rename("tmpReads_1.csfasta","$fileoutone") if (-e "tmpReads_1.csfasta");
  rename("tmpReads_2.csfasta","$fileouttwo") if (-e "tmpReads_2.csfasta");
  rename("tmpReads_1.csfasta.qual","$fileoutone.qual") if (-e "tmpReads_1.csfasta.qual");
  rename("tmpReads_2.csfasta.qual","$fileouttwo.qual") if (-e "tmpReads_2.csfasta.qual");
}

sub getRunIdsFromSraStudyId {
    warn __PACKAGE__ . " getRunIdsFromSraStudyId is derped, check out getRunsForStudy";
    my ($studyId) = @_;

    my $utils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";

    my $esearch = "$utils/esearch.fcgi?api_key=f2006d7a9fa4e92b2931d964bb75ada85a08&db=sra&retmax=1&usehistory=y&term=$studyId";

    my $esearch_result;
    my $counter=0;
    while (! defined $esearch_result) {
	if ($counter==10) {
	    die "esearch error. Could not access $esearch\n";
	}
	sleep 1;
	$esearch_result = get($esearch);
	$counter++;
    }

    $esearch_result =~ /<Count>(\d+)<\/Count>/s;
    my $Count = $1; 

    $esearch_result =~ /QueryKey>(\d+)<\/QueryKey>/s;
    my $QueryKey = $1; 

    $esearch_result =~ /<WebEnv>(\S+)<\/WebEnv>/s;
    my $WebEnv = $1; 

    my @ids = &runEfetch($Count, $QueryKey, $WebEnv);

    if ($studyId =~ /^[SED]RX/ ) {
	return @ids;
    } else {
	my @single_id=();
	foreach my $a (@ids) {
	    if ($a->[1] eq $studyId) {
		push(@single_id,$a);
	    }
	}
	return @single_id;
    }
}

sub runEfetch {
  my ($Count, $QueryKey, $WebEnv) = @_;

  my $utils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";

  my $efetch = "$utils/efetch.fcgi?api_key=f2006d7a9fa4e92b2931d964bb75ada85a08&rettype=xml&retmode=text&retmax=$Count&db=sra&query_key=$QueryKey&WebEnv=$WebEnv";

  my $efetch_result;
  my $counter=0;
  while (! defined $efetch_result) {
      if ($counter==10) {
	  die "runEfetch error. Could not access $efetch.\n";
      }
      sleep 1;
      $efetch_result = get($efetch);
      $counter++;
  }

  my $root = XMLin($efetch_result);

  my @ids;
  my @expPa;

  my $ep = $root->{'EXPERIMENT_PACKAGE'};

  if(ref($ep) eq 'ARRAY') {
    foreach my $a (@{$ep}) {
      push(@expPa,$a);
    }
  } else {
    push(@expPa,$ep);
  }

  foreach my $e (@expPa) {
    my %sids;
    if(ref($e->{SAMPLE}) eq 'ARRAY') {
      foreach my $a (@{$e->{SAMPLE}}) {
        $sids{$a->{accession}}++;
      }
    } else {
      $sids{$e->{SAMPLE}->{accession}}++;
    }

    my @runsets;
    if(ref($e->{RUN_SET}) eq 'ARRAY') {
      foreach my $a (@{$e->{RUN_SET}}) {
        push(@runsets,$a);
      }
    } else {
      push(@runsets,$e->{RUN_SET});
    }

    foreach my $rs (@runsets){
      if(ref($rs->{RUN}) eq 'ARRAY') {
        foreach my $r (@{$rs->{RUN}}) {
          push(@ids, [join(",",keys%sids),$r->{accession},$r->{total_spots},$r->{total_bases},($e->{EXPERIMENT}->{DESIGN}->{LIBRARY_DESCRIPTOR}->{LIBRARY_LAYOUT}->{PAIRED}) ? "PAIRED" : "SINGLE"]);
        }
      } else {
        push(@ids, [join(",",keys%sids),$rs->{RUN}->{accession},$rs->{RUN}->{total_spots},$rs->{RUN}->{total_bases},($e->{EXPERIMENT}->{DESIGN}->{LIBRARY_DESCRIPTOR}->{LIBRARY_LAYOUT}->{PAIRED}) ? "PAIRED" : "SINGLE" ]);
      }
    }

   #($e->{EXPERIMENT}->{DESIGN}->{LIBRARY_DESCRIPTOR}->{LIBRARY_LAYOUT}->{SINGLE}) ? print STDERR "ssingle\n" : print STDERR "ppaired\n"; 

  }

  return @ids;
}


sub getRunIdsFromSraSampleId { 
 my ($sid) = @_; 

 my $utils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"; 

 my $esearch = "$utils/esearch.fcgi?api_key=f2006d7a9fa4e92b2931d964bb75ada85a08&db=sra&retmax=1&usehistory=y&term=$sid"; 

 my $esearch_result;
 my $counter=0;
 while (! defined $esearch_result) {
     if ($counter==10) {
	 die "esearch error. Could not access $esearch\n";
     }
     sleep 1;
     $esearch_result = get($esearch);
     $counter++;
 }

 $esearch_result =~ /<Count>(\d+)<\/Count>/s;
 my $Count = $1; 

 $esearch_result =~ /QueryKey>(\d+)<\/QueryKey>/s;
 my $QueryKey = $1; 

 $esearch_result =~ /<WebEnv>(\S+)<\/WebEnv>/s;
 my $WebEnv = $1; 

 my @ids = &runEfetch($Count, $QueryKey, $WebEnv);

 if ($sid =~ /^[SED]RX/ ) {
     return @ids;
 } else {
     my @single_id=();
     foreach my $a (@ids) {
	 if ($a->[1] eq $sid) {
	     push(@single_id,$a);
	 }
     }
     return @single_id;
 }
}

##do wget, then split files with fastq-dump, delete .sra when complete and also the barcode if relevant
## wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR340/SRR340131/SRR340131.sra
## https://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR016/SRR016080/SRR016080.sra
## looks like a problem to construct .. look for full path in sample.xml
## fastq-dump --split-files SRR340224.sra

sub getFastqForSraRunId {
    my($runId,$isPairedEnd,$deflineVars) = @_;

    my $file = "$runId.sra";
    print STDERR "getFastqForSraRunId run $runId paired? $isPairedEnd\n";

    unlink("$runId.sra") if -e "$runId.sra";
    ###  replacing wget with prefetch       my $cmd = "wget https://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/".substr($runId,0,3)."/".substr($runId,0,6)."/$runId/$file";
    my $cmd = "prefetch -X 9999999999999 -O . $runId";
    #TODO test with something paired end, since not 100% sure how those files are named
    if (! -f "${runId}_1.fastq") {
      system($cmd);
    }
    if($?){
	die "ERROR ($?): Unable to fetch sra file for $runId\n";
    }
    if ($isPairedEnd) {
        if (defined $deflineVars) {
          system("fastq-dump -B -I --defline-seq '$deflineVars' --split-files ./$file");
        } else {
	  system("fastq-dump -B -I --split-files ./$file");
        }
    }
    else {
        if (defined $deflineVars) {
          system("fastq-dump --defline-seq '$deflineVars' --split-spot --skip-technical ./$file");
        } else {
  	  system("fastq-dump --split-spot --skip-technical ./$file");
        }
	my $first = $runId.".fastq";
	my $second = $runId."_1.fastq";
#  if (-e $first) {
	#     next;
	# }
	#else {
	system("mv $first $second");
    }
    my @files = glob("$runId*.fastq");
    print STDERR "DONE: ".scalar(@files)." files (".join(", ",@files).")\n";
    unlink($file);
}

sub getCsForSraRunId {
  my($runId,$pe) = @_;
  my $file = "$runId.sra";
  unlink("$runId.sra") if -e "$runId.sra";
  #####  replacing wget with prefetch      my $cmd = "wget https://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/".substr($runId,0,3)."/".substr($runId,0,6)."/$runId/$file";
  my $cmd = "prefetch -X 9999999999999 -O . $runId";
  print STDERR "retrieving $runId with $cmd\n";
  system($cmd);
  if($?){
    die "ERROR ($?): Unable to fetch sra file for $runId\n";
  }
  print STDERR "extracting fastq file(s)...";
  system("abi-dump ./$file");
  my @files = glob("$runId*.*");
  print STDERR "DONE: ".scalar(@files)." files (".join(", ",@files).")\n";
  unlink($file);
}

1;
