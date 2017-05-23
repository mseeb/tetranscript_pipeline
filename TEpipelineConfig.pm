#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) TEpipelineConfig.pm
##  Author:
##      <onishi@ie-freiburg.mpg.de>
#
##  Description:
##      This is the main configuration file for the TEtranscript
##      pipeline.  Before you can run the programs included
##      in this package you will need to edit this file and
##      configure for your site.
##
#******************************************************************************
#*
###############################################################################
package TEpipelineConfig;

use FindBin;
require Exporter;

@EXPORT_OK = qw( $SCRIPT_DIR $SCRIPT_LIB_DIR $SAMTOOLS
    $TRIM_PROG $STAR_PROG $star_mm10 $star_hg38 $star_mm10L1 $HISAT_PROG $hisat2_mm10 $hisat2_hg38 $hisat2_mm10L1 $TETRANS_PROG
    $tetrans_hg38_gene $tetrans_hg38_repeat $tetrans_mm10L1_gene $tetrans_mm10L1_repeat
    $tetrans_mm10_gene $tetrans_mm10_repeat $DEEP_PROG $genome_length_hg38 $genome_length_mm10);

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );
@ISA         = qw(Exporter);

BEGIN {
##----------------------------------------------------------------------##
##     CONFIGURE THE FOLLOWING PARAMETERS FOR YOUR INSTALLATION         ##
##----------------------------------------------------------------------##
## TEpipeline Location
## ======================
## The path to the TETranscript programs and support files
##
##    i.e. Typical UNIX installation
##     $REPEATMASKER_DIR = "/usr/local/RepeatMasker";
##    Windows w/Cygwin example:
##     $REPEATMASKER_DIR = "/cygdrive/c/RepeatMasker";
##
  $SCRIPT_DIR = "$FindBin::RealBin";
  $SCRIPT_LIB_DIR = "$SCRIPT_DIR/Libraries";

## Fastq Trimming
## ======================
## Determines whether reads are paired or single-end
## Trims fastq files
## Renames files
## Requires <fastq dir> <out dir> as input
  $TRIM_PROG  = "source /package/cutadapt-1.8.1/bin/activate && /package/trim_galore_v0.4.0/trim_galore";

## samtools
  $SAMTOOLS = "/package/samtools-1.3.1/bin/samtools";

## STAR
## ======================
  $STAR_PROG = "/package/STAR-2.5.2b/bin/STAR";
  $star_hg38 = "/data/repository/organisms/GRCh38_ensembl/STARIndex";
  $star_mm10 = "/data/repository/organisms/GRCm38_ensembl/STARIndex";
  $star_mm10L1 = "/data/repository/organisms/GRCm38_ensembl/STARIndex";

## HISAT2
## ======================
##
  $HISAT_PROG = "/package/hisat2-2.0.4/hisat2";
  $hisat2_mm10 = "/data/repository/organisms/GRCm38_ensembl/HISAT2Index/genome";
  $hisat2_hg38 = "/data/repository/organisms/GRCh38_ensembl/HISAT2Index/genome";
  $hisat2_mm10L1 = "/data/jenuwein/group/onishi/ref_seqs/indexes/hisat2/mm10L1";

  ## TETranscript
  ## ======================
  ##
  $TETRANS_PROG = "/data/jenuwein/group/onishi/software/tetoolkit-master/bin/TEtranscripts";
  $tetrans_hg38_gene = "/data/repository/organisms/GRCh38_ensembl/gencode/release_24/genes.gtf";
  $tetrans_hg38_repeat = "/data/jenuwein/group/onishi/ref_seqs/hg38/tetranscript/hg38_rmsk_TE.nochr.gtf";
  $tetrans_mm10_gene = "/data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.gtf";
  $tetrans_mm10_repeat = "/data/jenuwein/group/onishi/ref_seqs/mm10/tetranscript/mm10_rmsk_TE.nochr.gtf.gz";
  $tetrans_mm10L1_gene = "/data/repository/organisms/GRCm38_ensembl/gencode/m9/genes.gtf";
  $tetrans_mm10L1_repeat = "/data/jenuwein/group/onishi/ref_seqs/mm10/tetranscript/mm10L1_rmsk_TE.nochr.gtf";
}

## deeptools
## =======================
##
$DEEP_PROG = "/data/jenuwein/group/onishi/software/miniconda/bin/bamCoverage";  #version2.4 installed here
$genome_length_hg38 = "2569321340";  #for hg38 unique regions via GEM mappability
$genome_length_mm10 = "2207072725";  #for mm10 unique regions via GEM mappability

1;
