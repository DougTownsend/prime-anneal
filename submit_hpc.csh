#!/bin/tcsh
#BSUB -n 32
#BSUB -W 10000
#BSUB -R select[avx]
#BSUB -q shared_memory
#BSUB -J test_eq
#BSUB -o stdout.%J
#BSUB -e stderr.%J