#!/bin/bash
set -eo pipefail

REGION=$1
#clean up directory if it exists
if [[ -d "${REGION}" ]]; then
  rm -rf ./${REGION}
fi

shift 1
OPTS=$@
echo "Running command: ../chopBAI ${OPTS} test.sorted.bam ${REGION}"
../chopBAI ${OPTS} test.sorted.bam ${REGION}

if [[ -d ${REGION} ]]; then
  #ok, directory exists
  cd ${REGION}
  if [[ ! -f test.sorted.bam ]]; then
    if [[ $OPTS =~ "--symlink" || $OPTS =~ "-s" ]]; then
      echo "Symlink not created!"
      exit 1
    else
      #we'll just create it ourselves
      ln -s ../test.sorted.bam
    fi

    samtools view test.sorted.bam ${REGION} > out.chopBAI.sam
    samtools view ../test.sorted.bam ${REGION} > out.samtools.sam
    diff -q  out.chopBAI.sam out.samtools.sam
  fi
else
  echo "Directory ${REGION} was not created."
  exit 1
fi
