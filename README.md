TEtranscript_pipeline version 1.0
module load samtools1.3.1
module load R/3.2.3
No fastq dir indicated

NAME
    TETranscript_pipeline.pl

SYNOPSIS
      TETranscript_pipeline.pl [-options] <fastq dir>

DESCRIPTION
    The options are:

    -h(elp)
        Detailed help

    Default settings are for mouse.

    -dir [directory name]
        Writes output to this directory (default is query file directory,
        "-dir ." will write to current directory).

    -g(enome) [mm10|mm10_L1|hg38]
        Use an alternate search engine to the default.

    -trim
        Trim fastq files

    -hisat2
        Align fastq files

    -s(amples) [sample info file]
        File containing list of "treatment" and "control" samples. Necessary
        for DESeq.

        Row1: Treatment1.bam Treatment2.bam (space in between)

        Row2: Control1.bam Control2.bam

        For multiple comparisons, list sample files, separated by commas.

    -countsOnly
        Produce only counts table with TEtranscript. Do not run DESeq2.

AUTHORS
    <onishi@ie-freiburg.mpg.de>

