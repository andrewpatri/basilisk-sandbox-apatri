#!/bin/bash
#$ -q chimica2.q
#$ -N HF0
#$ -cwd
#$ -V
#$ -pe mpi 16
#$ -o gasification.out
#$ -e gasification.err

set-gcc-13.2
export CC=/software/chimica2/tools/gcc/gcc-13.2.0/bin/gcc
export BASILISK=/home/chimica2/rcaraccio/basilisk/basilisk/src
export PATH=$BASILISK:$PATH

CC="mpicc -D_MPI=16" CFLAGS+="-DTRACE=2 -DLB_AUTO=1"  make slab-q.tst
