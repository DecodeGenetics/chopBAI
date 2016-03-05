chopBAI
=======

Chopping a BAM index file into pieces.

chopBAI implements a reduction of a BAM index file to a specified genomic interval. The resulting index is much smaller in size and is semantically equivalent to the complete index, in the sense that it will give the same answers for all queries to reads within the interval of interest. Queries outside of the interval do *not* result in an error; instead, the queried regions may appear as empty even if the BAM file contains aligned reads.

chopBAI supports reduction of BAM file indices in the BAI format only. CSI is currently not supported.

Installation
------------

1. Download the SeqAn library version 2.0.0 or later from https://github.com/seqan/seqan.
   You do *not* need to follow the SeqAn install instructions.
   You only need the directory `.../include/seqan` of the SeqAn core library with all its content.
2. Set the path to the SeqAn core library in the `Makefile` if you chose not to save it in the chopBAI source directory.
   For example, set it to `SEQAN_LIB=/home/<user>/libraries/` if your directory containing all seqan header files is
   `/home/<user>/libraries/seqan`.
3. Run `make` in the chopBAI directory. If everything is setup correctly, this will create the binary `chopBAI`.


Usage
-----

chopBAI expects a BAM file and a list of regions (genomic intervals) as parameters.
You can either specify one or more regions on the command line or pass a file listing one region per line:

    ./chopBAI [OPTIONS] BAM-FILE REGION1 [... REGIONn]
    ./chopBAI [OPTIONS] BAM-FILE REGION-FILE


Regions should be in the format `CHR:BEGIN-END`, e.g. `chr4:15000000-16000000`.

The program looks for a BAI file at `BAM-FILE.bai`.


Example use case
----------------

Assume we have given a bamfile and its index, e.g. `NA12878.bam` and `NA12878.bam.bai` downloaded from the 1000 genomes project and indexed using samtools, in the current directory:

    wget -O NA12878.bam ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
    samtools index NA12878.bam

We can create a reduced BAI file using

    ./chopBAI NA12878.bam 4:15000000-16000000

This will create the folder `4:15000000-16000000` and inside it the BAI file `NA12878.bam.bai`, which is reduced in size compared to the original index file:

    du -h NA12878.bam.bai 4\:15000000-16000000/NA12878.bam.bai
    
    8.3M    NA12878.bam.bai
    4.0K    4:15000000-16000000/NA12878.bam.bai

If you would like to use a tool that assumes the BAM and the BAI file to share a common prefix (such as `samtools`), you can use chopBAI's `-s` option to create a symbolic link in the output folder:

    ./chopBAI -s NA12878.bam 4:15000000-16000000

or create a symbolic link yourself, e.g.

    ln -s NA12878.bam 4:15000000-16000000/NA12878.bam
    
A samtools command then uses the original bam file but the reduced index, e.g.

    samtools view -c 4:15000000-16000000/NA12878.bam 4:15501000-15506000
    
    263

If we compare the difference in I/O, e.g. the read system calls, when using the full index and the reduced index, we see a reduction by 96.21 % for this particular example:

    strace -e trace=read -o full.log samtools view -c NA12878.bam 4:15501000-15506000 > /dev/null
    awk 'BEGIN {FS="="}{ sum += $2} END {print sum}' full.log
    
    8033492
    
    strace -e trace=read -o reduced.log samtools view -c 4:15000000-16000000/NA12878.bam 4:15501000-15506000 > /dev/null
    awk 'BEGIN {FS="="}{ sum += $2} END {print sum}' reduced.log
    
    304300


References
----------

Birte Kehr and Páll Melsted (2015).
chopBAI: BAM index reduction solves I/O bottlenecks in the joint analysis of large sequencing cohorts.
*bioRxiv*, 030825.

http://biorxiv.org/content/early/2015/11/06/030825


Contact
-------

For questions and comments contact birte.kehr [ at ] decode.is
