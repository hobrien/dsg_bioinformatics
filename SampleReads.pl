#!/usr/bin/perl -w

=head1 NAME

SampleReads.pl version 1, 23 June 2011

=head1 SYNOPSIS

SampleReads.pl -i read file -f format [-n number of reads | -c coverage -s genome size] [-o out file] 

=head1 DESCRIPTION


=head2 NOTES

This script will randomly subsample the indicated number of reads fro a fasta/fastq file. If the
-c option is used, it will calculate the number of reads needed for the indecated coverage.

At present, the script assumes paired reads. If the filename ends with _1.fa(st[aq]), it will also
use XX_2.fa(st[aq]). Otherwise, it will add _1 and _2 to the files

SEQUENCE MUST BE ON A SINGLE LINE

=head1 AUTHOR

 Heath E. O'Brien E<lt>heath.obrien-at-utoronto-dot-caE<gt>

=cut
####################################################################################################


use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $outfilename;
my $format = "fasta";
my $infilename;
my $coverage;
my $read_num;
my $read_len;
my $genome_size = 6112448;
my $help = 0;
my $paired = 0;                       #not implemented
my $count;                            #number of reads in each file

GetOptions(
'o|out|outfile:s' => \$outfilename,
'f|format:s' => \$format,
'i|in:s' => \$infilename,
'c|copy:s' => \$coverage,
'n|num|number:s' => \$read_num,
'l|len|length:s' => \$read_len,
't|tot|total:s' => \$count,
's|seq:s' => \$genome_size,
'p|paired|' => \$paired,  #paired reads in the same file
'h|help|?' => \$help,

);

my $usage = "type perldoc SampleReads.pl for help\n";
if( $help ) {
    die $usage;
}

$infilename = shift unless $infilename;
$read_num = shift unless $read_num;

my $infile2, my $outfile2;
(my $name, my $path, my $ext) = fileparse($infilename, qr/\.[^.]*/);
if ( $name =~ /(.*)_1$/ ) {
  $infile2 = $path . $1 . "_2" . $ext;
}
else {
  $infilename = $path . $name . "_1" . $ext; 
  $infile2 = $path . $name . "_2" . $ext;
}  
if ( $outfilename ) {
  (my $outname, my $outpath, my $outext) = fileparse($outfilename, qr/\.[^.]*/);
  $outname =~ /(.*)(_1)?$/;
  $outfilename = $outpath . $1 . "_1" . $outext; 
  $outfile2 = $outpath . $1 . "_2" . $outext; 
}
else { 
  if ( $name =~ /(.*)_1$/ ) {
    $outfilename = $path . $1 . "_sample_1" . $ext;
    $outfile2 = $path . $1 . "_sample_2" . $ext;
  }
  else {
    $outfilename = $path . $name . "_sample_1" . $ext;
    $outfile2 = $path . $name . "_sample_2" . $ext;
  }
}
####################################################################################################

#count number of reads in input_file
unless ( $count ) {
  $count = `wc -l < $infilename`; die "wc failed: $?" if $?; chomp($count);
  if ( $format =~ /fasta/i ) {$count = $count / 2;}
  elsif ( $format =~ /fastq/i){ $count = $count/4;}
  else {die "format $format not recognized\n$usage"; }
}

#determine read length to calculate read number from coverage
if ( $coverage) {
  if ( $read_num ) {die "both coverage and read number specified. pick one\n$usage\n";}
  else {
    unless ( $read_len) {
      my $wc = `tail -1 $infilename | wc`; die "wc failed: $?" if $?; chomp($wc);
      $wc =~ /(\d+)$/;
      $read_len = ($1 -1) * 2;
    }
    $read_num = int($genome_size * $coverage / $read_len );
  }
}

#create sorted array of random read selections
if ( $read_num > $count ) { die "requested read number is greater than the number of reads in the file\n$usage";}
my $x = 1;
my %num_hash;
while ( $x <= $read_num ) {
  my $rand_num = int(rand($count));
  unless ( $num_hash{$rand_num} ) {
    $num_hash{$rand_num} = 1;
    $x ++;
  }
}
my @lines = sort { $a <=> $b } keys %num_hash;

if ($format =~ /fasta/i ) {
  PrintFasta($infilename, $outfilename, @lines);
  PrintFasta($infile2, $outfile2, @lines);
}
else {
  PrintFastq($infilename, $outfilename, @lines);
  PrintFastq($infile2, $outfile2, @lines);
}

sub PrintFasta {
  my $infile = shift;
  open(IN, "<$infile");
  my $outfile = shift;
  open(OUT, ">$outfile");
  $x = 0;
  while (<IN>) {
    unless ( @_ ) { next; }
    if ( $x == $_[0] ) {
      print OUT "$_";
    }
    elsif ( $x == $_[0] + 0.5 ) {
      print OUT "$_";
      shift(@_);
    }
    $x += 0.5;
  }
}

sub PrintFastq {  
  my $infile = shift;
  open(IN, "<$infile");
  my $outfile = shift;
  open(OUT, ">$outfile");  
  $x = 0;
  while (<IN>) {
    unless ( @_ ) { next; }
    if ( $x >= $_[0] ) {
      print OUT "$_";
      if ( $x == $_[0] + 0.75 ) { shift(@_); }
    }
    $x += 0.25;
  }
}
