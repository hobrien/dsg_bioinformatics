#!/usr/bin/perl -w

=head1 NAME

GetT3SE.pl version 1, 26 June 2012

=head1 SYNOPSIS

GetT3SE.pl -u user_name -i blast_file (tabular format) -s strain [-p password]

=head1 DESCRIPTION

-Parses results of a blast of all T3SE amino acid sequences against a genome sequence, identifies
non-overlapping hits, then updates a MySQL database for the strain in question to add T3SEs annotations

-It will then add the nucleotide sequences of the T3SEs to a file of sequences for that family,
align them based on the amino acid translations, and examine the resulting alignment to determine
how many alignment gaps the sequence is adding to the alignment

Options:
    -u user_name for database
    -i filename (tabular blast output)
    -s strain (name must match the namespace stored in the db)
    -p password for db (will prompt for one if it is not specified).

=head2 NOTES

It will be easy to add functionality to do the blasting as part of the script

In cases where the start codon is between contigs within a scaffold, the script will pick the
last start before the contig break. This will be detected during the post processing

This is working pretty well to add annotations to the DB. I can easily add the code from AlignT3SE.pl
to write the seqs to file and align them. One I do this I can use CountGaps to figure out how well
they align.

=head1 AUTHOR

Heath E. O'Brien E<lt>heath.obrien-at-utoroanto-dot-caE<gt>

=cut

####################################################################################################


use strict;
use warnings;
use Getopt::Long;
use Bio::SearchIO;
use Bio::DB::SeqFeature::Store;
use Term::ReadKey;
use Bio::Range;
use Bio::SeqIO;
use Bio::AlignIO;
use File::Basename;
use List::Util qw[min max];

my $adaptor  = 'DBI::mysql';
#my $database = 'test';
my $database = 'bioperl_features';
my $host     = 'localhost';
my $user;
my $pass;
my $infilename;
my $strain;


GetOptions( 'adaptor=s'  => \$adaptor,
	    'database=s' => \$database,
	    'host=s'     => \$host,
	    'user=s'     => \$user,
	    'strain=s' => \$strain,
	    'password=s' => \$pass,
	    'infile=s' => \$infilename
	    )
  or die "failure to communicate\n";

unless ( $pass ) { $pass = GetPassword(); }

unless ( $strain ) { die "Strain name not provided\n"; }
my $NAMESPACE = $strain;

my $dsn = "$database:$host";

warn "using '$dsn'\n";

my $db = Bio::DB::SeqFeature::Store->
  new( -adaptor => $adaptor,
       -dsn     => $dsn,
       -user    => $user,
       -pass    => $pass,
    -namespace => $NAMESPACE,
     );

$infilename = shift unless $infilename;
my $blastIO = new Bio::SearchIO(-format => 'blasttable',
                           -file   => $infilename);

(my $name, my $path, my $in_ext) = fileparse($infilename, qr/\.[^.]*/);
my $logname = $path . $name ."_log.txt";

#my $strain_name = GetPhy(GetOriginal($NAMESPACE));
print "Parsing blast hits...\n";
my %t3se_homo = %{BinHits($blastIO)};
#foreach ( keys ( %t3se_homo ) ) { print "$_\n"; }
print "Adding features to database...\n";
UpdateDB($db, \%t3se_homo);
print "Writing sequences to files...\n";
#WriteSeqs($strain, $db, \%t3se_homo);
print  "Printing summary...\n";
SeqSummary($db, $logname, \%t3se_homo);

sub SeqSummary {
  my $db = shift;
  my $logname = shift;
  my $hash_ref = shift;
  my %t3se_homo = %{$hash_ref};
  open(LOG, ">$logname");
  foreach my $t3se ( keys ( %t3se_homo ) ) {
    my @features = $db->get_features_by_alias($t3se);
    my $gene;
    foreach my $feature ( @features ) {
      if ($feature->type =~ /gene/i ) {
        if ( $gene ) { print STDERR "multiple copies of $t3se in database\n"; }
        else { $gene = $feature; }
      }
    }
    unless ( $gene ) { print STDERR "no gene feature for $t3se\n"; next; }
    if ( $gene->has_tag('note') and ($gene->get_tag_values('note'))[0] =~ /frameshift/i ) {
      print LOG "$t3se\tpseudo\n";
    }
    else { print LOG "$t3se\n"; }
  }
}

#bin blast hits into different contigs
sub SortHits {
  my $blastIO = shift;
  my %hits;
  while ( my $result = $blastIO->next_result ) {
    my $query = $result->query_name;
    my $multi_copy = 0;
    my @hsps;
    my $num_hits;
    while ( my $hit = $result->next_hit ) {
      while ( my $hsp = $hit->next_hsp ) {
        $hsp->hit->location->seq_id($hit->name);
        $hsp->location->seq_id($result->query_name);
        if ( $hits{$hit->name} ) { push(@{$hits{$hit->name}}, $hsp); }
        else { $hits{$hit->name} = [ $hsp ]; }
      }
    }
  }
  return \%hits;
}

sub GetDif {
  my $feature = shift;
  my $hsp = shift;
  my @diff = ($feature->start - $hsp->hit->start, $feature->end - $hsp->hit->end);
  unless ( $feature->strand == 1 ) { @diff = reverse(@diff); $diff[0] = -$diff[0]; $diff[1] = -$diff[1];}
  return \@diff;
}

sub CountGaps {
  my $t3se = shift;
  my $seqname = shift;
  my $start_gaps = 0;
  my $start_insertions = 0;
  my $end_gaps = 0;
  my $end_insertions = 0;
  my $all_gaps = 0;
  my $all_insertions = 0;
  my $conserved = 0;
  my $alnfilename = "/Users/heo3/Bioinformatics/T3SE/DNA/Alignment/Fasta/$t3se.fa";
  my $alnfile = new Bio::AlignIO(-format => 'fasta',
                           -file   => $alnfilename);
  my $aln = $alnfile->next_aln;
  my $length = $aln->length;
  my @seqs;
  foreach ($aln->each_seq){
    unless ( $_->display_id eq $seqname ) {push(@seqs, $_->seq); }
  }
  my $query = ($aln->each_seq_with_id($seqname))[0];
  for ( my $x = 0; $x < $length; $x ++ ) {
    my $query_gap = 0;
    my $gaps = 0;
    if ( substr($query->seq,$x,1) eq '-' ) { $query_gap = 1; }
    for ( my $y = 0; $y < @seqs; $y ++ ) {
      if ( substr($seqs[$y],$x,1) eq '-' ) { $gaps ++; }
    }
    if ( $query_gap and $gaps == 0 ) {  #unique gap in query
      $all_gaps ++;
      unless ( $conserved ) { $start_gaps ++; }
    }
    elsif ( $gaps == @seqs ) {  #gap in all other seqs
      $all_insertions ++;
      unless ( $conserved ) { $start_insertions ++; }
    }
    else { $conserved ++; }
  }
  for ( my $x = $length; $x >= 0; $x -- ) {
    my $query_gap = 0;
    my $gaps = 0;
    if ( substr($query->seq,$x,1) eq '-' ) { $query_gap = 1; }
    for ( my $y = 0; $y < @seqs; $y ++ ) {
      if ( substr($seqs[$y],$x,1) eq '-' ) { $gaps ++; }
    }
    if ( $query_gap and $gaps == 0 ) {  #unique gap in query
      $end_gaps ++;
    }
    elsif ( $gaps == @seqs ) {  #gap in all other seqs
      $end_insertions ++;
    }
    else { last; }
  }
  $all_gaps = $all_gaps - $start_gaps - $end_gaps;
  $all_insertions = $all_insertions - $start_insertions - $end_insertions;
  my @results;
  if ( $start_gaps ) { push(@results, "first " . $start_gaps . "bp missing"); }
  if ( $start_insertions ) { push(@results, $start_insertions . " bp extra at start"); }
  if ( $end_gaps ) { push(@results, "last " . $end_gaps . "bp missing"); }
  if ( $end_insertions ) { push(@results, $end_insertions . " bp extra at end"); }
  if ( $all_gaps ) { push(@results, $all_gaps . " internal bp missing"); }
  if ( $all_insertions ) { push(@results, $all_insertions . " internal bp extra"); }
  return @results;
}

sub GetPassword {
  my $key = 0;
  my $pass = "";
  print "\nPlease input your password: ";
  # Start reading the keys
  ReadMode(4); #Disable the control keys
  while(ord($key = ReadKey(0)) != 10) {
  # This will continue until the Enter key is pressed (decimal value of 10)
  # For all value of ord($key) see http://www.asciitable.com/
    if(ord($key) == 127 || ord($key) == 8) {
        # DEL/Backspace was pressed
        #1. Remove the last char from the password
        chop($pass);
        #2 move the cursor back by one, print a blank character, move the cursor back by one
        print "\b \b";
    } elsif(ord($key) < 32) {
        # Do nothing with these control characters
    } else {
        $pass = $pass.$key;
        print "*";
    }
  }
  print "\n\n";
  ReadMode(0); #Reset the terminal once we are done
  return ($pass);
}

sub FindEnd {
  my $contig = shift;
  my $hsp = shift;
  my $stop;
  if ( $hsp->hit->strand == 1 ) {
    my $x = $hsp->hit->end;
    while ($x < length($contig) ) {
      if ( substr($contig, $x - 3, 3) =~ /(TAA)|(TAG)|(TGA)/ ) { $stop = $x; last; }
      if ( substr($contig, $x - 3, 3) =~ /NNN/ ) { last; }
      $x +=3;
    }
  }
  else {
    my $x = $hsp->hit->start;
    while ($x > 2 ) {
      if ( substr($contig, $x-1, 3) =~ /(TTA)|(CTA)|(TCA)/ ) { $stop = $x; last; }
      if ( substr($contig, $x-1, 3) =~ /NNN/ ) { last; }
      $x -= 3;
    }
  }
  print "End: $stop\n";
  return $stop;
}

sub GetTag {
  my $db = shift;
  my %names;
  foreach ( $db->get_features_by_type('gene') ){
    my @load_ids = $_->each_tag_value('locus_tag');
    if ( $load_ids[0] ) {
      $load_ids[0] =~ /([^.]+)/;
      $names{$1} = 1;
    }
  }
  (keys ( %names))[0] =~ /(.*_)/;
  my $prefix = $1;
  my $x = 1;
  while ( $names{$prefix . sprintf  "%04d", $x} ) { $x ++; }
  return $prefix . sprintf  "%04d", $x;
}

sub FindStart {
  my $contig = shift;
  my $hsp = shift;
  my $start;
  if ( $hsp->hit->strand == 1 ) {
    my $x = $hsp->hit->start;
    my $target_start = $x - ($hsp->query->start * 3) + 1;
    print "Target Start: $x - (", $hsp->query->start, " * 3) + 1, = $target_start\n";
    while ($x > 1 ) {
      if ( substr($contig, $x-1, 3) =~ /(TAA)|(TAG)|(TGA)|(NNN)/ ) { last; }
      if ( substr($contig, $x-1, 3) =~ /(AT[ACGT])|([CGT]TG)/ ) {
        unless ($start and abs($x - $target_start) > abs($start - $target_start ) ) { $start = $x; }
      }
      $x -=3;
    }
  }
  else {
    my $x = $hsp->hit->end;
    my $target_start = $x + ($hsp->query->start * 3) - 1;
    print "Target Start: $x + (", $hsp->query->start, " * 3) - 1, = $target_start\n";
    while ($x < length($contig) ) {
      if ( substr($contig, $x - 4, 3) =~ /(ATT)|(CTA)|(TCA)|(NNN)/ ) { last; }
      if ( substr($contig, $x - 4, 3) =~ /([ACGT]AT)|(CA[ACGT])/ ) {
        unless ($start and abs($x - $target_start) > abs($start - $target_start ) ) { $start = $x; }
      }
      $x -=3;
    }
  }
  print "Start: $start\n";
  return $start;
}

#update first range of start and end with union of ranges of second hsp
sub AdjustRange {
  my $hsp1 = shift;
  my $hsp2 = shift;
  (my $start, my $stop, my $strand) = $hsp1->hit->location->union($hsp2->hit->location);
  $hsp1->hit->location->start($start);
  $hsp1->hit->location->end($stop);
  ($start, $stop, $strand) = $hsp1->query->location->union($hsp2->query->location);
  $hsp1->query->location->start($start);
  $hsp1->query->location->end($stop);
  return $hsp1;
}

sub CheckEnds {
  my $t3se = shift;
  my $seqname = shift;
  my $start;
  my $end;
  my $novel_gaps;
  my $novel_bases;
  my $alnfilename = "/Users/heo3/Bioinformatics/T3SE/DNA/Alignment/Fasta/$t3se.fa";
  my $alnfile = new Bio::AlignIO(-format => 'fasta',
                           -file   => $alnfilename);
  my $aln = $alnfile->next_aln;
  foreach ($aln->each_seq){
    if ( $_->display_id eq $seqname ) {
      if ( $_->seq =~ /^---/ ) { $start = 'missing'; }
      if ( $_->seq =~ /---$/ ) { $end = 'missing' }
    }
    else {
      unless (  $_->seq =~ /^---/ or $start ) { $start = 'match'; }
      unless (  $_->seq =~ /---$/ or $end ) { $end = 'match'; }
    }
  }
  unless ( $start ) { $start = 'long'; }
  unless ( $end ) { $end = 'long'; }
  return "Start: $start, End: $end";
}

sub GetRange {
  my $feature = shift;
  my @range = ($feature->start, $feature->end);
  unless ( $feature->strand == 1 ) { @range = reverse(@range); }
  return @range;
}

sub GetOverlap {
  my $feature = shift;
  my $hsp = shift;
#  print "HSP Loaction: ", $hsp->hit->location, "\n";
#  print "Feature Loaction: ", $feature->location, "\n";
  (my $start, my $stop, my $strand) = $hsp->hit->location->intersection($feature->location);
  return abs($start -$stop);
#  if ( $feature->strand != $hsp->hit->strand ) { $overlap = 'wrong strand'; }
#  elsif ( ($feature->start - $hsp->hit->start) % 3 ) { $overlap = 'wrong frame'; }
#  elsif ( ${GetDif($feature, $hsp)}[1] < 0 ) { $overlap = "truncation"; }
#  elsif ( ${GetDif($feature, $hsp)}[0] > 0 ) { $overlap = "wrong start"; }
  #else { $overlap = "contians blast hit"; }
#  return $overlap;
}

#get hash of individual genomic intervals that are homologous to T3SEs (overlapping ranges are combined)
sub BinHits {
  my $blastIO = shift;
  my %hits = %{SortHits($blastIO)};
  my %t3se_homo;
  foreach ( keys %hits ) {
    my @contig = @{$hits{$_}};
    while (my $hsp = pop(@contig) ) {
      my $print = 1;
      my $x = 0;
      while ($x < @contig) {
        if ( $hsp->hit->location->overlaps($contig[$x]->hit->location)) {  #if hsps overlap
          if ( $hsp->bits < $contig[$x]->bits ) {	                     #if new hsp has higher score than current one (discard current result)
            $contig[$x] = AdjustRange($contig[$x], $hsp);
            $print = 0;    #do not print result
            $x = @contig;  #skip to end of hsp array
          } #if ( $hsp->bits < $contig[$x]->bits ) {
          else {                   #if current hsp has higher score than new one (discard new one from array)
            $hsp = AdjustRange($hsp, $contig[$x]);
            splice(@contig, $x, 1);
          } #else {
        } #if ( $hsp->hit->location->overlaps($contig[$x]->hit->location)) {
        else { $x ++; }
      } #while ($x < @contig) {
      if ( $print ) {
        if ( $hsp->location->seq_id =~ /(mltB)|(LysM)|(PSPTO)/ ) { next; }
        unless ( $hsp->location->seq_id =~ /_([^'_]+)'?$/ ) { die $hsp->location->seq_id, " not parsed correctly\n" ; }
        my $name = $1;
        $name =~ s/-\d//;
        if ( $t3se_homo{$name} ) {
          $t3se_homo{$name . '-1'} = $t3se_homo{$name};
          delete $t3se_homo{$name};
        } #if ( $t3se_homo{$name} ) {
        if ( $t3se_homo{$name . '-1'} ) {
          my $x = 2;
          while ( $t3se_homo{$name . '-' . $x} ) { $x ++; }
          $t3se_homo{$name . '-' . $x} = $hsp;
        } #if ( $t3se_homo{$name . '-1'} ) {
        else { $t3se_homo{$name} = $hsp; }
      } #if ( $print ) {
    } #while (my $hsp = pop(@contig) ) {
  } #foreach ( keys %hits ) {
  return \%t3se_homo;
}

sub FixGaps {
  my $seq = shift;
  my $id = shift;
  if ($seq =~ /[^ACGTN]/i ) { die "$id: non-canonical bases in seq\n"; }
  while ( length($seq) % 3 ) { $seq =~ s/NNN(?!.*NNN)/NN/; }
  my @segments;
  while ( $seq =~ s/([ACGT]+)(N+)//i ) { push(@segments, ($1, $2)); }
  $seq =~ s/([ACGT]{3})$//;
  push(@segments, $seq, $1);
  my $length = length($segments[0]);
  for ( my $x = 2; $x < @segments; $x += 2 ) {
    my $subseq = $segments[$x];
    $length += length($segments[$x-1]);
    my $frame1stop;
    my $frame2stop;
    my $frame3stop;
    $subseq =~ s/(TGA)|(TAA)|(TAG)/XXX/g;
    while  ( $subseq =~ s/([^X]+XXX)// ) {
      $length += length($1);
      if ( $length % 3 == 0 ) { $frame1stop ++; }
      elsif ( ($length + 2 ) % 3 == 0 ) { $frame2stop ++; }
      else { $frame3stop ++ ; }
    }
    $length += length($subseq);
    if ( $frame1stop ) {
      if ( $frame2stop ) {
        if ( $frame3stop ) { die "$id: Stop codons in all frames\n"; }
        else { $segments[$x-1] = substr($segments[$x-1], 2); $length -= 2; }    #frame 3 is correct
      }
      elsif ( $frame3stop ) {                                                   #frame 2 is correct
        $segments[$x-1] = substr($segments[$x-1], 1);
        $length -= 1;
      }
      else { warn "$id: multiple frames without stop codons\n"; }
    }
  }
  return join("", @segments);
}

sub WriteSeqs {
  my $strain = shift;
  my $db = shift;
  my $hash_ref = shift;
  my %t3se_homo = %{$hash_ref};
  foreach my $t3se ( keys ( %t3se_homo ) ) {
    $t3se =~ /([a-zA-Z]+)/;
    my $outfilename = "/ntfs/Genomics_resources/T3SE/$1.fa";
    my @features = $db->get_features_by_alias($t3se);
    my $gene;
    foreach my $feature ( @features ) {
      if ($feature->type =~ /gene/i ) {
        if ( $gene ) { print STDERR "multiple copies of $t3se in database\n"; }
        else { $gene = $feature; }
      }
    }
    unless ( $gene ) { print STDERR "no gene feature for $t3se\n"; next; }
    my $seq = $gene->seq;
    if ( $gene->has_tag('note') and ($gene->get_tag_values('note'))[0] =~ /frameshift/i ) { $outfilename =~ s/.fa/_FS.fa/; }
    elsif ( $seq->seq =~ /NNN/ ) { $seq->seq(FixGaps($seq->seq, $t3se)); }
    $seq->display_id($strain . "_" . $t3se);
    my $outfile = Bio::SeqIO->new('-file' => ">>$outfilename",
         '-format' => 'fasta' ) or die "could not open seq file $outfilename\n";
    $outfile->write_seq($seq);
  }
}

sub UpdateDB { #identify gene and CDS features in Database corresponding to blast hits and add aliases
  my $db = shift;
  my $hash_ref = shift;
  my %t3se_homo = %{$hash_ref};
  foreach my $t3se ( sort ( keys %t3se_homo ) ) {
    my $hsp = $t3se_homo{$t3se};
    my $cds;
    my $gene;
    my @locus_tags;
    my $location = $hsp->hit->location->seq_id . ", " . join('-',$hsp->range('hit'));
    #print "looking for features\n";
    foreach ( $db->get_features_by_location($hsp->hit->location->seq_id, $hsp->range('hit'))) {
      print GetOverlap($_, $hsp), " > 50\n";
      if ( GetOverlap($_, $hsp) > 50 ) {
        @locus_tags = (@locus_tags, $_->get_tagset_values('locus_tag'));
        if ( $_->type =~ /gene/i or $_->type =~ /CDS/i ) {$db->delete(($_)); }
      }
    }
    my $locus_tag;
    if ( @locus_tags ) { $locus_tag = $locus_tags[0]; print "$locus_tag\n"; }
    else { $locus_tag = GetTag($db); }
    AddGene($db, $hsp, $locus_tag, $t3se);
  }
}
#this will determine the codon_start for features that start before/after the contig
sub GetCodonStart {
  my $hsp = shift;
  my $contig = shift;
  my $codon_start = 1;
  if ( $hsp->hit->strand == 1 ) {
    while ( ($hsp->hit->start - 1 + $codon_start ) % 3 ) { $codon_start ++; }
  }
  else {
    while ( (length($contig) - 1 - $hsp->hit->end + $codon_start) % 3 ) { $codon_start ++; }
  }
  return $codon_start;
}

#This will correct sequences with gaps so that they will translate correctly.
#Some of the same logic an be used to figure out if I have a pseudogene or not
sub CheckORF {
  #print "checking ORF\n";
  my $seq = shift;
  my $id = shift;
  $seq = substr($seq, 0, length($seq) - 3);                           #remove stop codon
  my $no_pseudo = 1;
  if ($seq =~ /[^ACGTN]/i ) { die "non-canonical bases in seq\n"; }
  #while ( length($seq) % 3 ) { $seq =~ s/NNN(?!.*NNN)/NN/; }
  my @segments;
  #print "removing Ns\n";
  while ( $seq =~ s/([ACGT]+)N+//i ) { push(@segments, $1); }
  push(@segments, $seq);
  foreach my $subseq ( @segments ) {
    my $length;
    my $frame1stop;
    my $frame2stop;
    my $frame3stop;
    $subseq =~ s/(TGA)|(TAA)|(TAG)/XXX/g;
    while  ( $subseq =~ s/([^X]+XXX)// ) {
      $length += length($1);
      if ( $length % 3 == 0 ) { $frame1stop ++; }
      elsif ( ($length + 2 ) % 3 == 0 ) { $frame2stop ++; }
      else { $frame3stop ++ ; }
    }
    if ( $frame1stop and $frame2stop and $frame3stop ) {
      $no_pseudo = 0;
    }
  }
  return $no_pseudo;
}

sub RevCom {
  my $seq = shift;
  $seq =~ tr/ACGTacgt/TGCAtgca/;
  $seq = reverse $seq;
  return $seq;
}

sub AddGene {
  my $db = shift;
  my $hsp = shift;
  my $tag = shift;
  my$t3se = shift;
  my $t3_type;
  if ( $t3se =~ /hrp/ or $t3se =~ /hopP'/ or $t3se =~ /hopA[JK]/ ) { $t3_type = 'Type III helper'; }
  elsif ( $t3se =~ /hop[JL]/ or $t3se =~ /hopA[CJNP]/ or $t3se =~ /hopAH2/ ) { $t3_type = 'Discontinued type III effector'; }
  else { $t3_type = 'Type III effector'; }
  #print "Adding $t3se\n";
  print "getting ", $hsp->hit->location->seq_id, "\n";
  my $contig = $db->fetch_sequence($hsp->hit->location->seq_id);
  #my @contigs = $db->segment($hsp->hit->location->seq_id);
  #if ( @contigs > 1 ) {  $contig = $contigs[0]->seq; }
  #else { $contig = $db->segment($hsp->hit->location->seq_id)->seq; }
  my $contig_break = 0;
  my $codon_start = 1;
  my $end = FindEnd($contig, $hsp);
  my $start = FindStart($contig, $hsp);
  unless ( $end ) {
    if ( $hsp->hit->strand == 1 ) { $end = $hsp->hit->end; }
    else { $end = $hsp->hit->start; }
    $contig_break = 1;
  }
  unless ( $start ) {
    if ( $hsp->hit->strand == 1 ) { $start = $hsp->hit->start; }
    else { $start = $hsp->hit->end; }
    $contig_break = 1;
    $codon_start = GetCodonStart($hsp, $contig);
  }
  #print "start = $start, end = $end\n";
  if ( $start > $end ) {
    ($start, $end) = ($end,$start);
  }
  my $cds = substr($contig, $start, $end - $start + 1);
  unless ( $hsp->hit->strand == 1 ) {  $cds = RevCom($cds); }
  if ( $cds =~ /NN/) { $contig_break = 1; }
  my $note;
  unless ( CheckORF($cds, $t3se) ) { $note = "$t3_type $t3se disrupted by one or more putative frameshift or nonsense mutations"; }
  if ( $contig_break ) {
  if ( $note ) { $note .= 'and contig break'; }
    else { $note = '$t3_type $t3se disrupted by contig break'; }
  }
  my %attributes = ( 'Alias' => $t3se,
                     'locus_tag' => $tag
                   );
  if ( $note ) { $attributes{'note'} = $note; }
  #print "Adding feature\n";
  $db->new_feature(-seq_id => $hsp->hit->location->seq_id,
                   -start => $start,
                   -end => $end,
                   -primary_tag => 'gene',
                   -display_name => $tag,
                   -strand => $hsp->hit->strand,
                   -source => 'tBLASTn',
                   -score => $hsp->score,
                   -load_id => $tag . '.gene',
                   -phase => 1,
                   -index => 1,
                   -attributes => \%attributes,
                  );
  unless ( $note ) {
    $attributes{'product'} = "$t3_type " . ucfirst($t3se);
    $attributes{'codon_start'} = $codon_start;
    $db->new_feature(-seq_id => $hsp->hit->location->seq_id,
                     -start => $start,
                     -end => $end,
                     -primary_tag => 'CDS',
                     -display_name => $tag,
                     -strand => $hsp->hit->strand,
                     -source => 'tBLASTn',
                     -score => $hsp->score,
                     -load_id => $tag,
                     -phase => 1,
                     -index => 1,
                     -attributes => \%attributes,
                    );
  }
}
