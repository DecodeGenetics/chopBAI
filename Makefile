# Use version 4.8.2 of g++
CXX=/opt/rh/devtoolset-2/root/usr/bin/g++

# include SeqAn libraries
CXXFLAGS+=-I/nfs/odinn/users/birtek/code/libraries -DSEQAN_HAS_ZLIB=1 -lz -pthread

DATE=on `date +'%Y-%m-%d %H:%M'`
CXXFLAGS+=-DDATE=\""$(DATE)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros

# DEBUG build
#CXXFLAGS+=-g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0

# set std to c++0x to allow using 'auto' etc.
CXXFLAGS+=-std=c++0x
