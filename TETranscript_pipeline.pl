#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) TETranscript_pipeline.pl
##  Author:
##      onishi@ie-freiburg.mpg.de
##  Description:
##	Pipeline to take in fastq files, trim, align to designated genome,
##	get counts over genes and repeats using TETranscript, as well as
##	create bigwig files for visualization.
##
##
#******************************************************************************
#

=head1 NAME

TETranscript_pipeline.pl

=head1 SYNOPSIS

  TETranscript_pipeline.pl [-options] <fastq dir>

=head1 DESCRIPTION

The options are:

=over 4

=item -h(elp)

Detailed help

=back

Default settings are for mouse.

=over 4

=item -dir [directory name]

Writes output to this directory (default is query file directory,
"-dir ." will write to current directory).

=back

=over 4

=item -g(enome) [mm10|mm10_L1|hg38]

Use an alternate search engine to the default.

=back

=over 4

=item -trim

Trim fastq files

=back

=over 4

=item -hisat2

Align fastq files

=back

=over 4

=item -s(amples) [sample info file]

File containing list of "treatment" and "control" samples.
Necessary for DESeq.

Row1: Treatment1.bam Treatment2.bam (space in between)

Row2: Control1.bam Control2.bam

For multiple comparisons, list sample files, separated by commas.

=back

=over 4

=item -countsOnly

Produce only counts table with TEtranscript. Do not run DESeq2.

=back

=head1 AUTHORS

<onishi@ie-freiburg.mpg.de>

=cut

#
# Module Dependence
#
use strict;
use warnings;
use FindBin;
use lib $FindBin::RealBin;
use Carp;
use Getopt::Long;
use POSIX qw(:sys_wait_h);
use Storable qw(nstore retrieve);
use File::Copy;
use File::Spec;
use File::Path;
use Data::Dumper;
use Cwd 'abs_path';
use File::Basename;

# Libraries
use TEpipelineConfig;


my $version = "1.0";
print "TEtranscript_pipeline version $version\n";

if ( $ARGV[ 0 ] && $ARGV[ 0 ] eq '-v' ) {
  exit( 0 );
}

my $cmdLine = $0 . join( ' ', @ARGV );

my $load_samtools = "echo \"module load samtools1.3.1\"";
#print "\t\t$load_samtools\n";
system($load_samtools);

my $load_R = "echo \"module load R/3.2.3\"";
system($load_R);

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @opts = qw( genome|g=s dir=s trim hisat2 debug=s samples|s=s countsOnly );

#
# Get the supplied command line options, and set flags
#
my %options = ();
unless ( &GetOptions( \%options, @opts ) ) {
  exec "pod2text $0";
  exit( 0 );
}

# Print the internal POD documentation if something is missing
if ( $#ARGV == -1 && !$options{'help'} ) {
  print "No fastq dir indicated\n\n";

  # This is a nifty trick so we don't have to have
  # a duplicate "USAGE()" subroutine.  Instead we
  # just recycle our POD docs.  See PERL POD for more
  # details.
  exec "pod2text $0";
  die;
}

#
# Get the date
#
my $date = localtime( time() );

# Debugging flag
my $DEBUG = 0;
$DEBUG = 1 if ( exists $options{'debug'}
                && $options{'debug'} & 1 );

# Windows does not support the use of ":" in a filename.
$date =~ s/[ ,\t,\n:]//g;

#load samtools
my $samtools = "$TEpipelineConfig::SAMTOOLS";

# get output dir
my $outdir;
if (defined $options{'dir'}) {
  $outdir = abs_path($options{'dir'});
  mkdir $outdir if ( ! -d $outdir );
} else {
  $outdir = dirname(abs_path($0));
}

print "\t\toutdir is $outdir\n";

# Setup the genome
my $genomeToUse;

if ( defined $options{'genome'} && $options{'genome'} =~ m/mm10|hg38|mm10_L1/ ){
    $genomeToUse = $options{'genome'};
    print "\t\t$genomeToUse\n"
 } else {
   print "No genome selected (mm10, hg38, or mm10_L1)\n\n";
   exec "pod2text $0";
   die;
}

my $hisat2_index; my $gene_gtf; my $rep_gtf;

if ( $genomeToUse =~ m/hg38/ ) {
  $hisat2_index = "$TEpipelineConfig::hisat2_hg38";
  $gene_gtf = "$TEpipelineConfig::tetrans_hg38_gene";
  $rep_gtf = "$TEpipelineConfig::tetrans_hg38_repeat";
} elsif ( $genomeToUse =~ m/mm10_L1/ ){
    $hisat2_index = "$TEpipelineConfig::hisat2_mm10L1";
    $gene_gtf = "$TEpipelineConfig::tetrans_mm10L1_gene";
    $rep_gtf = "$TEpipelineConfig::tetrans_mm10L1_repeat";
  } elsif ( $genomeToUse =~ m/mm10/ ){
    $hisat2_index = "$TEpipelineConfig::hisat2_mm10";
    $gene_gtf = "$TEpipelineConfig::tetrans_mm10_gene";
    $rep_gtf = "$TEpipelineConfig::tetrans_mm10_repeat";
  } else {
    print "\t\tNo Genome Specified!!\n";
    exit;
  }

  print "\t\thisat2 index: $hisat2_index\n";

##get list of fastq files
my @samplesR1;
my $paired = 1; # paired-end = 1, single-end =0

my $inputdir = $ARGV[0];
opendir(DIR, "$inputdir") || die "Cannot open DIR $inputdir";
my @fqfiles = readdir(DIR);
closedir DIR;

for my $file (@fqfiles) {
  if ( $file =~ /fastq|fq|fastq.gz|fq.gz$/ && $file =~ /_R1/ ) {
    push (@samplesR1, $file);
    #check if R2 file exists
    my $file2 = $file;
    $file2 =~ s/_R1/_R2/;
    if ( ! -e $inputdir."/".$file2 ) {
      $paired = 0;
      print "$file2\n";
    }
  }
}

if ( $paired == 1) {
  print "\t\tpaired-end detected\n";
} else {
  print "\t\tsingle-end detected\n";
}

##################################################################
## trim if necessary
##################################################################
my $trim_dir = $outdir."/00_trimmed_fq";
if ( defined $options{'trim'} ){
    mkdir $trim_dir if ( ! -d $trim_dir);
    print "\t\t---------------------------------\n\t\tFQ_TRIMMING\n";
    my $prog = "$TEpipelineConfig::TRIM_PROG";
    for my $fileR1 (@samplesR1 ) {
      if ( $fileR1 =~ /fastq.gz$/ ){
        my $sample = $fileR1;
        $sample =~ s/(.+)_R1.+/$1/;
        my $fileR2 = $fileR1;
        $fileR2 =~ s/R1/R2/;
        my $outfile = $trim_dir."/".$sample."_R1.fastq.gz";
        if ( ! -e $outfile ){
          my $trim_cmd;
          if ( $paired == 0 ){
            $trim_cmd = "$prog --paired --stringency 2 -o $trim_dir $inputdir/$fileR1 >> $trim_dir/log.txt 2>&1";
          } else {
            $trim_cmd = "$prog --paired --stringency 2 -o $trim_dir $inputdir/$fileR1 $inputdir/$fileR2 >> $trim_dir/log.txt 2>&1";
          }
          print "\t\t$trim_cmd\n";
          system($trim_cmd);

          my $tempfile = $trim_dir."/".$sample."_R1_val_1.fq.gz";
          my $mv1_cmd = "mv $tempfile $outfile";
          print "\t\t$mv1_cmd\n";
          system($mv1_cmd);

          if ( $paired == 1 ) {
            my $tempfile2 = $trim_dir."/".$sample."_R2_val_2.fq.gz";
            my $outfile2 = $trim_dir."/".$sample."_R2.fastq.gz";
            my $mv2_cmd = "mv $tempfile2 $outfile2";
            print "\t\t$mv2_cmd\n";
            system($mv2_cmd);
          }
        }
      }
    }
    $inputdir = $trim_dir;
}

##################################################################
## align with hisat2 (use hisat2 as opposed to star because soft-
## clipping can be turned off)
## use --no-softclip to turn off softclipping
## report up to 100 positions with -k 100
## output will be unsorted!
##################################################################

my $libtype = "RF"; #same as fr-firststrand for Illumina
$libtype = "R" if ( $paired == 0 );
my $numk = 100;

if ( defined $options{'hisat2'} ){
  my $hisat2_dir = $outdir."/01_hisat2";
  mkdir $hisat2_dir if ( ! -d $hisat2_dir );

  print "\t\t---------------------------------\n\t\tHISAT2\n\n";
  my $hisat2_prog = "$TEpipelineConfig::HISAT_PROG";
  for my $fileR1 ( @samplesR1 ){
    if ( $fileR1 =~ /fastq.gz$|fq.gz$|fastq$|fq$/ ) {
      my $sample = $fileR1;
      $sample =~ s/(.+)_R1.+/$1/;
      my $fileR2 = $fileR1;
      $fileR2 =~ s/R1/R2/;
      my $hisat_outdir = $hisat2_dir."/".$sample;
      mkdir $hisat_outdir;
      my $summary_file = $hisat_outdir."/align_summary.txt";
      my $outbam = $hisat_outdir."/".$sample.".bam";

      my $linkbam = $hisat2_dir."/".$sample.".bam";
      if ( ! -e $linkbam ) {
        my $hisat2_cmd;
        if ( $paired == 1) {
          $hisat2_cmd = "$hisat2_prog -p 8 -x $hisat2_index --rna-strandness $libtype -1 $inputdir/$fileR1 -2 $inputdir/$fileR2 --novel-splicesite-outfile $hisat_outdir/splice_sites.txt -k $numk --sp 1000,1000 2> $summary_file | $samtools view -uh - | $samtools sort -n -@ 2 -m 2G - -o $outbam";
        } else {
          $hisat2_cmd = "$hisat2_prog -p 8 -x $hisat2_index --rna-strandness $libtype -1 $inputdir/$fileR1 --novel-splicesite-outfile $hisat_outdir/splice_sites.txt -k $numk --sp 1000,1000 2> $summary_file | $samtools view -uh - | $samtools sort -n -@ 2 -m 2G - -o $outbam";
        }

        print "\t\t$hisat2_cmd\n";
        system($hisat2_cmd);

      #add to header
#      my $header = $hisat_outdir."/header.sam";
#      my $header_cmd = "cat <($samtools view -H $outbam) <(echo \'\@PG\tCL:\"$hisat2_cmd\"\') > $header";
#      print "\t\t$header_cmd\n";
#      system($header_cmd);
#      next;
#      my $tempbam = $hisat_outdir."/".$sample.".reheader.bam";
#      my $add_header_cmd = "$samtools reheader $header $outbam | $samtools view - -bo $tempbam";
#      print "\t\t$add_header_cmd\n";
#      system($add_header_cmd);

#      my $move_cmd = "mv $tempbam $outbam";
#      print "\t\t$move_cmd\n";
#      system($move_cmd);

#      my $cleanup_cmd = "rm $header";
#      print "\t\t$cleanup_cmd\n";
#      system($cleanup_cmd);

#      my $index_cmd = "$samtools index $outbam";
#      print "\t\t$index_cmd\n";
#      system($index_cmd);

        my $link_cmd = "ln -s $outbam $linkbam";
        print "\t\t$link_cmd\n";
        system($link_cmd);
      }
	   }
  }
}

##################################################################
## TETranscript using treatment/control from samples list
##################################################################
if ( defined $options{'samples'} ){
  print "\t\t---------------------------------\n\t\tTETRANSCRIPT\n\n";
  my $te_dir = $outdir."/02_TEtranscript";
  mkdir $te_dir if (! -d $te_dir);
  my $strand = "reverse";
  my $minRead = 100;

  my $teprog = "$TEpipelineConfig::TETRANS_PROG";
  
  my @sampleFiles = split(/,/, $options{'samples'});

  foreach my $sampleFile (@sampleFiles) {
	print "$sampleFile\n";
	  my $sampleName = basename($sampleFile);
	  print "\t\t$sampleName\n";
	  $sampleName =~ s/(.+)\.txt/$1/;

	  my $sample_dir = $te_dir."/".$sampleName;
	  mkdir $sample_dir if ( ! -d $sample_dir);

	  my $projectName = $sample_dir."/".$sampleName;

	  my $tsamples; my $csamples;

	  #get sample info
	  #first line is list of treatment bam files
	  #second line is list of control bam files
	  my $count=0;
	  open FH, "$sampleFile" or die "cannot open samples file $sampleFile $!\n";
	  while (my $line = <FH> ){
	    chomp $line;
	    if ( $count == 0 ) {
	      $tsamples=$line;
	      $count++;
	    } elsif ( length($line) > 0 ){
	      $csamples=$line;
	    }
	  }

	  #run TEtranscripts
	  my $summary_out = $sample_dir."/".$sampleName."_summary.txt";
	  my $tet_cmd;

	  if ( defined $options{'countsOnly'} ){
		  $tet_cmd = "$teprog -t $tsamples -c $csamples --GTF $gene_gtf --TE $rep_gtf --format BAM --mode multi --minread $minRead --stranded $strand --project $projectName --countsOnly 2> $summary_out";
	  } else {
		  $tet_cmd = "$teprog -t $tsamples -c $csamples --GTF $gene_gtf --TE $rep_gtf --format BAM --mode multi --minread $minRead --stranded $strand --project $projectName 2> $summary_out";
	  }
	  print "\t\t$tet_cmd\n";
	  system($tet_cmd);
	}
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
