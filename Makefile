# Use version 4.8.2 of g++
CXX=g++
CC=$(CXX)

# Set this to include SeqAn libraries, either system wide
# or download into current folder and set to .
SEQAN_LIB=/nfs/odinn/users/birtek/code/libraries
#SEQAN_LIB=.


CXXFLAGS+=-I$(SEQAN_LIB) -DSEQAN_HAS_ZLIB=1
LDLIBS=-lz -lpthread
DATE=on `date +'%Y-%m-%d %H:%M'`
CXXFLAGS+=-DDATE=\""$(DATE)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros

# DEBUG build
#CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0


all: chopBAI

chopBAI: chopBAI.o 

chopBAI.o: chopBAI.cpp bam_index_csi.h

clean:
	rm -f *.o chopBAI
