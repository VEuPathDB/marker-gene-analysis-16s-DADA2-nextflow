#!/usr/bin/perl

package Utils;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(commafy runCmd timestampLog appendToLogFile timeformat mail getProgramDir);

use strict;
use File::Basename;

sub commafy {
  my $Rv = shift;

  1 while $Rv =~ s/^(-?\d+)(\d{3})/$1,$2/;

  return $Rv;
}

sub runCmd {
    my ($cmd) = @_;

    my $output = `$cmd`;
    my $status = $? >> 8;
    die "Failed with status $status running '$cmd'\n" if $status;
    return $output;
}

sub runCmdTag {
  my $Tag = shift;
  my $Cmd = join(' ', @_);

  print STDERR join("\t", scalar(localtime), $Tag, $Cmd), "\n";

  my $output = `$Cmd`;
  my $status = $? >> 8;
  die "Failed with status $status running '$Cmd'\nOutput was $output" if $status;
  return $output;
}

sub appendToLogFile {
    my ($logFile, $msg) = @_;
    open(FILE, ">>$logFile") || die "couldn't open logFile $logFile";
    print FILE $msg;
    close(FILE);
}

sub timestampLog {
  print STDERR join("\t", scalar(localtime), $$, @_), "\n";
}

# convert secs to hr:min:sec
sub timeformat {
    my($tot_sec) = @_;

    my $secs = $tot_sec % 60;
    my $mins = int($tot_sec / 60) % 60;
    my $hours = int($tot_sec / 3600);

    $secs = "0$secs" if $secs < 10;
    $mins = "0$mins" if $mins < 10;
    $hours = "0$hours" if $hours < 10;
    return "$hours:$mins:$secs";
}

sub mail {
  my ($to, $from, $subject, $body, $cc, $bcc) = @_;

  unless($to && $from && $subject && $body) {
    die "mail usage:  CBIL::Util::Utils::mail(to, from, subject, body, cc, bcc)";
  }

  my $ccString = $cc ? "-c $cc" : "";
  my $bccString = $bcc ? "-b $bcc" : "";

  my $cmd = "export EMAIL=$from; echo '$body'|mutt -s '$subject' $ccString $bccString $to";

  my $result = system($cmd);

  unless($result / 256 == 0) {
    die "ERROR running the following command:  $cmd";
  }
}

##method to return the directory where the actual executable is found ... ie follows sym_links
sub getProgramDir {
  my $program = shift;
  return dirname($program) if $program =~ /^\//;     # return full path if provided 
  my $pgm = basename($program);
  my $wh = `which $pgm`;
  chomp $wh;
  my $link = readlink($wh);
  return dirname($link);
}


1;
