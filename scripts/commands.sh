#!/usr/bin/env bash

# This file has all the commands used to assebmle all datasets. It's specific to the cluster available to be (with
# SLURM controller). It's not meant to be used, just for consultation purposes

# IDBA
srun --pty --nodes=1 --cpus-per-task=16 --mem-per-cpu=2 --tasks-per-node=1 shifter --image=cimendes/idba:31.12.2016-1
# merge (uncompressed) reads and convert to fasta
fq2fa --merge ERR2984773_1.fq ERR2984773_2.fq reads.fa


# MEGAHIT
srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=cimendes/megahit-assembler:12.08.19-1
# even
/NGStools/megahit/bin/megahit -1 ERR2935805_1.fq.gz -2 ERR2935805_2.fq.gz -o out -t 16
# log
/NGStools/megahit/bin/megahit -1 ERR2984773_1.fq.gz -2 ERR2984773_2.fq.gz -o out_ERR298477 -t 16



# MAPPING
srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=mcfonsecalab/minimap2:latest
for file in $(ls /home/cimendes/Binning_assessment/metagenomic-assembler-comparison/data/references/ZymoBIOMICS.STD.refseq.v2/Genomes/*triple_chromosome.fasta); do minimap2 -c -t 16 -r 10000 -g 10000 -x asm20 --eqx $file out_ERR2935805/final.contigs.fa > ERR2935805_$(basename $file).paf; done
for file in $(ls /home/cimendes/Binning_assessment/metagenomic-assembler-comparison/data/references/ZymoBIOMICS.STD.refseq.v2/Genomes/*triple_chromosome.fasta); do minimap2 -c -t 16 -r 10000 -g 10000 -x asm20 --eqx $file out_ERR298477/final.contigs.fa  > ERR298477_$(basename $file).paf; done

