chopBAI
=======

Chopping a BAM index file into pieces.


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

Assume we have given the file `sequence.bam` and `sequence.bam.bai` in the current directory.
We can create a reduced BAI file using

    ./chopBAI sequence.bam chr4:15000000-16000000

This will create the folder `chr4:15000000-16000000` with the reduced BAI file `sequence.bam.bai` inside it.

If you would like to use a tool that assumes the BAM and the BAI file to share a common prefix (such as `samtools`),  you can create a symbolic link, e.g

    ln -s sequence.bam chr4:15000000-16000000/sequence.bam
    samtools view -c chr4:15000000-16000000/sequence.bam chr4:15500000-15600000

which then uses the original bam file, but the reduced index.

References
----------

Birte Kehr and Páll Melsted (2015).
chopBAI: BAM index reduction solves I/O bottlenecks in the joint analysis of large sequencing cohorts.
*bioRxiv*, 030825.

http://biorxiv.org/content/early/2015/11/06/030825


Contact
-------

For questions and comments contact birte.kehr [ at ] decode.is
