#!/usr/bin/env bash

# This file has all the commands used to assebmle all datasets. It's specific to the cluster available to be (with
# SLURM controller). It's not meant to be used, just for consultation purposes

# IDBA
srun --pty --nodes=1 --cpus-per-task=40 --mem-per-cpu=6 --tasks-per-node=1 shifter --image=cimendes/idba:31.12.2016-1
# merge (uncompressed) reads and convert to fasta
fq2fa --merge ERR2984773_1.fq ERR2984773_2.fq reads.fa
idba_ud -l reads.fa --num_threads 16 -o out # -l reads longer than 128 nucleotides. sample had 150nt

srun --pty --nodes=1 --cpus-per-task=40 --mem-per-cpu=16 --tasks-per-node=1 shifter --image=loneknightpy/idba
idba_ud -l reads.fa --num_threads 16 -o out2

#bbtools
reformat.sh in1=mockSample_fwd_shuffled.fastq.gz in2=mockSample_rev_shuffled.fastq.gz out=mockSample_reads.fasta
idba_ud -l mockSample_reads.fasta --num_threads 40 -o out


# MEGAHIT
srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=cimendes/megahit-assembler:12.08.19-1
# even
/NGStools/megahit/bin/megahit -1 ERR2935805_1.fq.gz -2 ERR2935805_2.fq.gz -o out -t 16
# log
/NGStools/megahit/bin/megahit -1 ERR2984773_1.fq.gz -2 ERR2984773_2.fq.gz -o out_ERR298477 -t 16



# MAPPING
srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=cimendes/minimap2:2.17-1
minimap2 -c -t 4 -r 10000 -g 10000 -x asm20 --eqx --secondary=no $file out_ERR2935805/final.contigs.fa > ERR2935805_$(basename $file).paf


srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=cimendes/metaspades:11.10.2018-1
metaspades.py -o out2 -1 mockSample_fwd_shuffled.fastq.gz -2 mockSample_rev_shuffled.fastq.gz --only-assembler -t 16 -m 32


srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=cimendes/snowball:24.09.2017-1
export PYTHONPATH=/NGStools/snowball:$PYTHONPATH
python2 /NGStools/snowball/algbioi/ga/run.py -f mockSample_fwd_shuffled.fastq.gz -s mockSample_rev_shuffled.fastq.gz -m /NGStools/Pfam-A.hmm -o out.fna -a -p 16 -i 225


# SKESA
srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=flowcraft/skesa:2.3.0-1
skesa --cores 16 --fastq <1> <2> --use_paired_ends > out.fasta


#SPAdes
srun --pty --nodes=1 --tasks-per-node=1 --cpus-per-task=16 --mem-per-cpu=2GB shifter --image=cimendes/metaspades:11.10.2018-1
spades.py -o out -t 16 --only-assembler --careful -1 mockSample_fwd_shuffled.fastq.gz -2 mockSample_rev_shuffled.fastq.gz


minimap2 -x sr --secondary=no ../final.contigs.fa /home/cimendes/Binning_assessment/Mock_in_silico/mockSample_fwd_shuffled.fastq /home/cimendes/Binning_assessment/Mock_in_silico/mockSample_rev_shuffled.fastq


for file in $(ls *.fasta); do filename="${file%%.*}"; minimap2 -c -t 4 -r 10000 -g 10000 -x asm20 --eqx --secondary=no /home/ines/git/metagenomic-assembler-comparison/data/references/Zymos_Genomes_triple_chromosomes.fasta $file > ${filename}.paf; done

bbsplit.sh in=../fastq/ERR2935805_1.fq.gz in2=../fastq/ERR2935805_2.fq.gz ref=Escherichia_coli_plasmid.fa,Staphylococcus_aureus_plasmid1.fa,Staphylococcus_aureus_plasmid2.fa,Staphylococcus_aureus_plasmid3.fa,Cryptococcus_neoformans_draft_genome.fasta,Saccharomyces_cerevisiae_draft_genome.fa basename=bbsplit%.fq.gz out1=bbsplit_ERR2935805_1.fq.gz out2=bbsplit_ERR2935805_2.fq.gz