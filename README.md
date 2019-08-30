#  Benchmarking of de novo metagenomic assembly software

While searching for a benchmark with relevant metagenomic assembly software, I gave up and decided to make my own. WORK IN PROGRESS 

## Table of Contents

TODO

## Introduction

TODO

## Methods

### Metagenomic datasets

### Assembly softwares and commands

I've followed the following [metagenomic assembly tools table](https://academic.oup.com/view-large/131667617), published by [Ayling et al. 2019](https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bbz020/5363831), as base for selecting metagenomic software to bee tested. In total 17 tools are presented. [IVA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4495290/) and [SAVAGE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411778/) were removed as they were aimed at viruses, as well as [Genovo](https://www.liebertpub.com/doi/abs/10.1089/cmb.2010.0244?rfr_dat=cr_pub%3Dpubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&journalCode=cmb) and [MAP](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts162) due to inavilability of software, [VICUNA](https://www.broadinstitute.org/viral-genomics/viral-genomics-analysis-software-registration) due to requiring registeration, [Omega](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu395) as it was an assembly pipeline, abd [BBAP](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5406902/), [MetaVelvet](http://metavelvet.dna.bio.keio.ac.jp/), [MEtaVelvet-SL](http://metavelvet.dna.bio.keio.ac.jp/MSL.html), [PRICE](http://derisilab.ucsf.edu/software/price/) and [Ray Meta](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-12-r122) due to no update sinde 2016. The following tools will be tested:

#### IDBA-UD
Published by [Peng et al. 2012](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts174), it's a De Brujin graph assembler for assembling reads from single-cell sequencing or metagenomic sequencing technologies with uneven sequencing depths. It employs multiple depthrelative thresholds to remove erroneous k-mers in both low-depth and high-depth regions. The technique of local assembly with paired-end information is used to solve the branch problem of low-depth short repeat regions. To speed up the process, an error correction step is conducted to correct reads of high-depth regions that can be aligned to highconfident contigs. The latest version is available at https://github.com/loneknightpy/idba, and an official docker image at https://hub.docker.com/r/loneknightpy/idba
Last update: 31/12/2016 (GitHub)

#### MegaGTA
Published by [Le et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5657035/), MegaGTA is a gene-targeted assembler that utilizes iterative de Bruijn graphs. It tries to improve on Xander assembler, implementing the same method of using the trained Hidden Markov Model (HMM) to guide the traversal of de Bruijn graph, but using mutiple k-mer sizes to take full advantage of multiple k-mer sizes to make the best of both sensitivity and accuracy. 
The latest version is available at https://github.com/HKU-BAL/MegaGTA. 
Last update: 16/06/2016 (GitHub)


#### MEGAHIT
MEGAHIT, published by [Li et al. 2015](https://academic.oup.com/bioinformatics/article/31/10/1674/177884), *de novo* assembler for assembling large and complex metagenomics data in a time- and cost-efficient manner. It makes use of succinct de Bruijn graph, with a a multiple k-mer size strategy. In each iteration, MEGAHIT cleans potentially erroneous edges by removing tips, merging bubbles and removing low local coverage edges,specially useful for metagenomics which suffers from non-uniform sequencing depths. 
The latest version is available at https://github.com/voutcn/megahit
Last update: 12/08/2019 (GitHub)

#### Snowball
Published by [Gregot et al. 2016](https://academic.oup.com/bioinformatics/article/32/17/i649/2450756), Snowball is a novel strain aware gene assembler for shotgun metagenomic data that does not require closely related reference genomes to be available. Like MegaGTA and Xander, it uses profile hidden Markov models (HMMs) of gene domains of interest to guide the assembly. 
The latest version is available at https://github.com/hzi-bifo/snowball
Last update: 24/09/2017 (GitHub)

### MetaSPAdes
SPAdes started out as a tool aiming to resolve uneven coverage in single cell genome data, but later metaSPAdes was released, building specific metagenomic pipeline on top of SPAdes. IT was published by [Nurk et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411777/), it uses multiple k-mer sizes of de Bruijn graph, starting with lowest kmer size and adding hypothetical kmers to connect graph. It's available at http://cab.spbu.ru/software/spades/ and https://github.com/ablab/spades
Last update: 11/10/2018 (release) and 24/04/2019 (GitHub)

### Xander
Like Snowball and MegaGTA, Xander employs HMM profile model to perform a guided assebly targeting specific genes. These are used to create a novel combined weighted assembly graph. Xander performs both assembly and annotation concomitantly using information incorporated in this graph. It was published by [Wang et al. 2015](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-015-0093-6) and it's available at https://github.com/rdpstaff/Xander_assembler.
Last update: 27/10/2017 (GitHub)
 
