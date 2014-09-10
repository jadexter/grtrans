#!/bin/bash
#th=(1 2 4 8 16 32)
#for i in `seq 1 6`;
#do
    export OMP_NUM_THREADS=$1
#echo $1
    ./grtrans
#done
