#!/bin/bash

set -eo pipefail

# Test all options
for opts in "" "--linear" "--symlink" "--linear --symlink"; do
  echo "Testing chopBAI with option $opts"


  for chr in chrA:B:C:D chrA:B chrB; do
    echo "  Testing chromosome $chr"
    ./test.sh $chr $opts
    for reg in "100" "2,000" "1-100" "1,000-10,000"; do
      chrreg="${chr}:${reg}"
      ./test.sh $chrreg $opts
    done
  done
done
