chopBAI
=======

Chopping a BAM index file into pieces.


Installation
------------

1. Download the SeqAn library version 2.0.0 or later from https://github.com/seqan/seqan.
   You do *not* need to follow the SeqAn install instructions.
   You only need the directory .../include/seqan of the SeqAn core library with all its content.
2. Set the path to the SeqAn core library in the 'Makefile' if it does not reside into the chopBAI source directory,
   e.g. to 'SEQAN_LIB=/home/<user>/libraries/' if it you have the seqan directory containing all header files stored in
   '/home/<user>/libraries/seqan/'.
3. Run 'make' in the chopBAI directory. If everything is setup correctly, this will create the binary 'chopBAI'.


Usage
-----

chopBAI expects a BAM file and a list of regions (genomic intervals) as paramters.
You can either specify one or more regions on the command line or pass a file listing one region per file:

    ./chopBAI [OPTIONS] BAM-FILE REGION1 [... REGIONn]
    ./chopBAI [OPTIONS] BAM-FILE REGION-FILE


Regions should be in the format 'CHR:BEGIN-END', e.g. 'chr4:15000000-16000000'.

The program looks for a BAI file at 'BAM-FILE.bai'.


References
----------

Birte Kehr and Páll Melsted (2015).
chopBAI: BAM index reduction solves I/O bottlenecks in the joint analysis of large sequencing cohorts.
Submitted.


Contact
-------

For questions and comments contact birte.kehr [ at ] decode.is
