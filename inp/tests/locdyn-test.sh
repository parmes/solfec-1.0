#!/bin/bash

solfec ./inp/tests/locdyn.py

mpirun -np 2 solfec-mpi ./inp/tests/locdyn.py

sort ./out/tests/locdyn-1 > ./out/tests/locdyn-1-sorted

sort ./out/tests/locdyn-2 > ./out/tests/locdyn-2-sorted

diff ./out/tests/locdyn-1-sorted ./out/tests/locdyn-2-sorted
