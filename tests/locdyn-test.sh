#!/bin/bash

CPUS="2 3 4"

echo -n "(processors: 1, "
solfec -s 1 ./tests/locdyn.py
sort ./out/tests/locdyn/locdyn-1 > ./out/tests/locdyn/locdyn-1-sorted

DIFF=""
for NCPU in $CPUS; do
  if [ $NCPU -lt "4" ]; then echo -n $NCPU", "
  else echo -n $NCPU") " 
  fi
  mpirun -np $NCPU solfec-mpi -s $NCPU ./tests/locdyn.py
  sort ./out/tests/locdyn/locdyn-$NCPU > ./out/tests/locdyn/locdyn-$NCPU-sorted
  D=`diff ./out/tests/locdyn/locdyn-1-sorted ./out/tests/locdyn/locdyn-$NCPU-sorted`
  DIFF="$DIFF$D"
done

if [ -n "$DIFF" ]; then
  echo "FAILED"
else
  echo "PASSED"
fi
