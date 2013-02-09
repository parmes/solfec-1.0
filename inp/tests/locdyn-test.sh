#!/bin/bash

CPUS="2 3 4 5 6 7 8"

echo -n "(processors: 1, "
solfec -s 1 ./inp/tests/locdyn.py
sort ./out/tests/locdyn/locdyn-1 > ./out/tests/locdyn/locdyn-1-sorted

DIFF=""
for NCPU in $CPUS; do
  if [ $NCPU -lt "8" ]; then echo -n $NCPU", "
  else echo -n $NCPU") " 
  fi
  mpirun -np $NCPU solfec-mpi -s $NCPU ./inp/tests/locdyn.py
  sort ./out/tests/locdyn/locdyn-$NCPU > ./out/tests/locdyn/locdyn-$NCPU-sorted
  D=`diff ./out/tests/locdyn/locdyn-1-sorted ./out/tests/locdyn/locdyn-$NCPU-sorted`
  DIFF="$DIFF$D"
done

if [ -n "$DIFF" ]; then
  echo "FAILED"
else
  echo "PASSED"
fi
