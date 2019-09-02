#!/usr/bin/env bash

# IDBA
srun --pty --nodes=1 --cpus-per-task=16 --mem-per-cpu=2 --tasks-per-node=1 shifter --image=cimendes/idba:31.12.2016-1
# merge (uncompressed) reads and convert to fasta
fq2fa --merge ERR2984773_1.fq ERR2984773_2.fq reads.fa

