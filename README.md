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
You can either specify one or more regions on the command line or pass a file listing one region per file:

    ./chopBAI [OPTIONS] BAM-FILE REGION1 [... REGIONn]
    ./chopBAI [OPTIONS] BAM-FILE REGION-FILE


Regions should be in the format `CHR:BEGIN-END`, e.g. `chr4:15000000-16000000`, one region per line.

The program looks for a BAI file at `BAM-FILE.bai`.

Running with chopped indices
----------------------------

Assuming we have the file `sequence.bam` and `sequence.bam.bai` in the current directory we can chop up the BAI file using

    ./chopBAI sequence.bam chr4:15000000-16000000

this will create the folder `chr4:15000000-16000000` with the smaller file `sequence.bam.bai` inside it. Many tools assume that the BAM and the BAI file share a common prefix. In that case you can get around this with a symbolic link, e.g

    cd chr4:15000000-16000000
    ln -s ../sequence.bam
    samtools view -c chr4:15500000-15600000 sequence.bam
    cd ../

which uses the original bam file, but the reduced index.

References
----------

Birte Kehr and Páll Melsted (2015).
chopBAI: BAM index reduction solves I/O bottlenecks in the joint analysis of large sequencing cohorts.
Submitted.


Contact
-------

For questions and comments contact birte.kehr [ at ] decode.is
