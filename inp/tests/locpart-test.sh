#!/bin/bash

CPUS="2 3 4 5 6 7 8"

echo -n "(processors: 1, "
solfec -s 1 ./inp/tests/locpart.py
sort ./out/tests/locpart/locpart-1 > ./out/tests/locpart/locpart-1-sorted

DIFF=""
for NCPU in $CPUS; do
  if [ $NCPU -lt "8" ]; then echo -n $NCPU", "
  else echo -n $NCPU") " 
  fi
  mpirun -np $NCPU solfec-mpi -s $NCPU ./inp/tests/locpart.py
  sort ./out/tests/locpart/locpart-$NCPU > ./out/tests/locpart/locpart-$NCPU-sorted
  D=`diff ./out/tests/locpart/locpart-1-sorted ./out/tests/locpart/locpart-$NCPU-sorted`
  DIFF="$DIFF$D"
done

if [ -n "$DIFF" ]; then
  echo "FAILED"
else
  echo "PASSED"
fi
