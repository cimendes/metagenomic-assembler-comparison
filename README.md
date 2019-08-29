#  Benchmarking of de novo metagenomic assembly software

While searching for a benchmark with relevant metagenomic assembly software, I gave up and decided to make my own. WORK IN PROGRESS 

## Table of Contents

TODO

## Introduction

TODO

## Methods

### Metagenomic datasets

### Assembly softwares and commands

I've followed the following [metagenomic assembly tools table](https://academic.oup.com/view-large/131667617), published by [Ayling et al. 2019](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz020/5363831), as base for selecting metagenomic software to bee tested. In total 17 tools are presented. [IVA]() and [SAVAGE]() were removed as they were aimed at viruses, as well as [Genovo]() and [MAP]() due to inavilability of software, [VICUNA]() due to requiring registeration, [Omega]() as it was an assembly pipeline, abd [BBAP](), [MetaVelvet](), [MEtaVelvet-SL](), [PRICE]() and [Ray Meta]() due to no update sinde 2016. The following tools will be tested:

#### IDBA-UD
Published by [Peng et al. 2012](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts174), it's a De Brujin graph assembler for assembling reads from single-cell sequencing or metagenomic sequencing technologies with uneven sequencing depths. It employs multiple depthrelative thresholds to remove erroneous k-mers in both low-depth and high-depth regions. The technique of local assembly with paired-end information is used to solve the branch problem of low-depth short repeat regions. To speed up the process, an error correction step is conducted to correct reads of high-depth regions that can be aligned to highconfident contigs. The latest version is available at https://github.com/loneknightpy/idba, and an official docker image at https://hub.docker.com/r/loneknightpy/idba
Last update: 31/12/2016 (GitHub)

#### MegaGTA
Published by [Le et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5657035/)


#### MEGAHIT

#### 

 
