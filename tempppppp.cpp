#!/bin/bash

if [ $# -lt 1 ]; then
  echo "Usage: $0 <executable_name_without_dot_out>"
  exit 1
fi

EXENAME=$1
INPDIR=../../inputs/selected-inputs
OUTDIR=out${EXENAME}-$(date +%d-%b-%Y-%H%M%S)

mkdir -p "$OUTDIR"

for file in $(ls -Sr $INPDIR/*.vrp); do
  fileName=$(basename "$file" .vrp)
  ./${EXENAME}.out "$file" > "$OUTDIR/$fileName.sol" 2>> "$OUTDIR/time.txt"
  echo "$file - Done"
done

sort "$OUTDIR/time.txt"
