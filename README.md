#  Benchmarking of *de novo* (meta)genomic assembly software

 :warning: WORK IN PROGRESS :warning:

## Table of Contents

* [Introduction](#introduction)
* [Methods](#methods)
    * [de novo Assembly tools](#de-novo-assembly-tools)
    * [Benchmark datasets](#benchmark-datasets)
    * [Assessing Assembly Success](#assessing-assembly-success)
* [Results](#resuls)
* [Authors](#authors)
* [References](#references)


## Introduction

Shotgun metagenomics can offer relatively unbiased pathogen detection and characterization, potentially able to provide 
genotyping, antimicrobial resistance and virulence profiling in a single methodological step. This comes with the cost 
of producing massive amounts of information that require expert handling and processing, as well as capable 
computational infrastructures. One of the biggest challenges when dealing with metagenomic data is the lack of gold 
standards, although major efforts are being made on the standardization and assessment of software, both commercial 
and open source (Angers-Loustau et al., 2018; Gruening et al., 2018; Sczyrba et al., 2017: Couto et al., 2018).

The most common strategy for analysing metagenomic data is through *de novo* assembly followed by a combination of 
different tools for characterization. The assembly methods provide longer sequences that are more informative than 
shorter sequencing data and can provide a more complete picture of the microbial community in a given sample. 

The contigs obtained in an assembly, ideally each collecting the sequences that belong to a single 
microorganism present in the sample, represents one of the greatest bottlenecks when trying to obtain fiable, 
reproducible results, not only in metagenomic samples but also in classical whole genome sequencing methodologies. 

Several **dedicated** metagenomic assembly tools for short-read data are available that, n comparison to traditional 
assemblers, are supposed to be better at dealing with the combination of intragenomic and intergenomic repeats and 
uneven sequencing coverage (Olson et al., 2017). The use of non dedicated assemblers for metagenomics may come with the 
cost of wrongly interpret variation as error, especially in samples that contained closely related species, and the 
construction of chimeric sequences (Teeling & Glockner, 2012) as traditional assemblers follow the basic principle 
that the coverage in a sample is constant. Despite these caveats, the distribution of the community is often unknown, 
and a reproducible comparison of methods has yet to be performed.

The *de novo* assembly is one of the key processesses when analysing (meta)genomic data as, in theory, it allows the 
re-build of complete genomes from a pool of raw sequence. To assess the performance and limitations of currently 
available de novo assembly algorithms, we propose a benchmarking of **currently available and recently mantained** 
*de novo* assembly tools, both traditional and dedicated for metagenomic data. 


## Methods

One of the main goals is to ensure that the results obtained in this assessment are reproducible through the many 
systems and environments available. Several steps that have been implemented to ensure the transparency and 
reproducibility of the results. Favouring open-source tools, with clear documentation describing the methodology 
implemented, and stating the version of the software used and which parameters were used enables the comparison of 
results. This is simplified by containerizing all the software tools with [Docker](https://www.docker.com/). The use of 
the Nextflow (Tommaso et al., 2017) workflow managers pushes reproducibility to the next level by taking advantage of 
containerization and scalability, enabling the workflow to be executed with the exact same parameters in the same 
conditions in a multitude of different environments.

The [scripts](scripts/) of the analysis of the assembly results are provided in this repository, alongside a [Jupyter Notebook]() 
describing the steps in which the scripts were used. 


### *de novo* Assembly tools

We've compiled a [collection](data/docs/de%20novo%20metagenomic%20assembly%20comparison%20-%20Tools%20available%20-%20short-read.pdf) 
of 24 *de novo* assembly tools, including 5 implementing Overlap, Layout and Consensus (OLC) assembly algorithm, 18 
De Bruijn graph assemblers, of which 10 are single k-mer value assemblers and 8 implement a multiple k-mer value 
approach, and a hybrid assembler using both OLC and single k-mer De Bruijn algorythms. Of these, 12 were developed 
explicitly to handle metagenomic datasets. 

The tools were ordered by **date of last update** and [Docker containers](docker/) for the top 11 assemblers were 
created with the **latest released version**. The version was used as tag. In the case of tools with no release, the 
container was compiled with latest version in the master branch, using the date of the last update as tag.     

An [assembly pipeline](nextflow/short) was developed to perform concurrent assemblies with all tools in the benchmark and obtain 
performance metrics for each one. This pipeline is implemented in [Nextflow](https://www.nextflow.io/). All resulting 
assemblies were filtered for a minimum contig length of 1000bp with BBTool's `reformat.sh`.

The command used to run the Nextflow pipeline was `nextflow run main.nf -profile slurm_shifter --fastq="fastq/*_{1,2}*`,
using nextflow version 19.04.1.5072. 

Assemblers benchmarked, ordered by date of last update:
* **Unicycler**
An assembly pipeline for bacterial genomes that can do long-read assembly, hybrid assembly and short-read assembly. 
When asseblying Illumina-only read sets where it functions as a SPAdes-optimiser, using a De Bruijn 
algorithm with multiple k-mer values. It was published by [Wick et al. 2017](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595).
    * Source code: [https://github.com/rrwick/Unicycler](https://github.com/rrwick/Unicycler)
    * Date of last update: 20/07/2020
    * Container: `cimendes/unicycler:0.4.8-1`

* **SPAdes**
A tool aiming to resolve uneven coverage in single cell genome data through multiple k-mer sizes of De Brujin graphs. It
starts with the smallest k-mer size and and adds hypotetical k-mers to connect graph. 
    * Source code: [http://cab.spbu.ru/software/spades/](http://cab.spbu.ru/software/spades/)
    * Date of last update: 20/06/2020
    * Container: `cimendes/spades:3.14.1-1`

* **MetaSPAdes**
[SPAdes](#spades) started out as a tool aiming to resolve uneven coverage in single cell genome data, but later metaSPAdes was 
released, building specific metagenomic pipeline on top of [SPAdes](#spades). It was published by [Nurk et al. 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411777/), 
and like SPAdes, it uses multiple k-mer sizes of de Bruijn graph, starting with lowest kmer size and adding hypothetical 
kmers to connect graph. It's available at http://cab.spbu.ru/software/spades/ and https://github.com/ablab/spades
    * Source code: [http://cab.spbu.ru/software/spades/](http://cab.spbu.ru/software/spades/)
    * Date of last update: 20/06/2020
    * Container: `cimendes/spades:3.14.1-1`

* **Minia**
This tool, published by [Chikhi & Rizk, 2013](https://almob.biomedcentral.com/articles/10.1186/1748-7188-8-22) in 
*Algorithms for Molecular Biology, performs the assembly on a data structure based on unitigs produced by the 
BCALM software and using graph simplifications that are heavily inspired by the [SPAdes](#spades) assembler.
    * Source code: [https://github.com/GATB/gatb-minia-pipeline](https://github.com/GATB/gatb-minia-pipeline)
    * Date of last update: 07/06/2020
    * Container: `cimendes/minia:3.3.4-1`

* **GATB-Minia Pipeline**
GATB-Minia is an assembly pipeline, still unpublished, that consists Bloocoo for error correction,  Minia 3 for contigs 
assembly, which is based on the BCALM2 assembler, and the BESST for scaffolding. It was developed to extend Minia 
assembler to use multiple k-mer values. It was developed to extend the [Minia](#minia) assembler to use De Bruijn 
algorithm with multiple k-mer values. It's explicit for metagenomic data. 
    * Source code: [https://github.com/GATB/gatb-minia-pipeline](https://github.com/GATB/gatb-minia-pipeline)
    * Date of last update: 31/05/2020
    * Container: `cimendes/gatb-minia-pipeline:31.05.2020-1`

* **BCALM2**
This assembler, publiched by [Chikhi et al, 2016](https://academic.oup.com/bioinformatics/article/32/12/i201/2289008) 
in *Bioinformatics*, is a fast and low memory algorithm for graph compaction, consisting of three stages: careful 
distribution of input k-mers into buckets, parallel compaction of the buckets, and a parallel reunification step to 
glue together the compacted strings into unitigs. It's a traditional single k-mer value De Bruijn assembler.
    * Source code: [https://github.com/GATB/bcalm](https://github.com/GATB/bcalm)
    * Date of last update: 22/05/2020
    * Container: `cimendes/bcalm:2.2.3-1`

* **SKESA**
This *de novo* sequence read assembler is based on De Bruijn graphs and uses conservative heuristics and is designed to
create breaks at repeat regions in the genome, creating shorter assemblies but with greater sequence quality. It tries
to obtain good contiguity by using k-mers longer than mate length and up to insert size. It was recently published by
[Souvorov et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1540-z) and it's available at
https://github.com/ncbi/SKESA.
    * Source code: [https://github.com/ncbi/SKESA/releases](https://github.com/ncbi/SKESA/releases)
    * Date of last update: 14/03/2020
    * Container: `flowcraft/skesa:2.4.0-1`

* **MEGAHIT**
MEGAHIT, published by [Li et al. 2015](https://academic.oup.com/bioinformatics/article/31/10/1674/177884), *de novo* 
assembler for assembling large and complex metagenomics data in a time- and cost-efficient manner. It makes use of 
succinct de Bruijn graph, with a a multiple k-mer size strategy. In each iteration, MEGAHIT cleans potentially erroneous 
edges by removing tips, merging bubbles and removing low local coverage edges,specially useful for metagenomics which 
suffers from non-uniform sequencing depths. 
    * Source code: [https://github.com/voutcn/megahit](https://github.com/voutcn/megahit)
    * Date of last update: 15/10/2019
    * Container: `cimendes/megahit-assembler:1.2.9-1`

* **PANDAseq**
This assembler, published by [Masella et al. 2012] in *BMC Bioinformatics*, implements an OLC algorithm to assemble 
genomic data. It align Illumina reads, optionally with PCR primers embedded in the sequence, and reconstruct an 
overlapping sequence. 
    * Source code: [https://github.com/neufeld/pandaseq](https://github.com/neufeld/pandaseq)
    * Date of last update: 30/08/2018
    * Container: `cimendes/pandaseq:2.11-1`

* **VelvetOptimier**
This optimizing pipeline, developed by Torsten Seeman, is still unpublished but extends the original [Velvet]() 
assembler by performing several assemblies with variable k-mer sizes. It searches a supplied hash value range for the 
optimum, estimates the expected coverage and then searches for the optimum coverage cutoff. It uses Velvet's internal 
mechanism for estimating insert lengths for paired end libraries. It can optimise the assemblies by either the default 
optimisation condition or by a user supplied one. It outputs the results to a subdirectory and records all its 
operations in a logfile.
    * Source code: [https://github.com/tseemann/VelvetOptimiser](https://github.com/tseemann/VelvetOptimiser)
    * Date of last update: 21/01/2017
    * Container: `cimendes/velvetoptimiser:2.2.6-1`

* **IDBA-UD**
Published by [Peng et al. 2012](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts174), 
it's a De Brujin graph assembler for assembling reads from single-cell sequencing or metagenomic sequencing technologies 
with uneven sequencing depths. It employs multiple depth relative thresholds to remove erroneous k-mers in both 
low-depth and high-depth regions. The technique of local assembly with paired-end information is used to solve the 
branch problem of low-depth short repeat regions. To speed up the process, an error correction step is conducted to 
correct reads of high-depth regions that can be aligned to high confidence contigs. 
    * Source code: [https://github.com/loneknightpy/idba](https://github.com/loneknightpy/idba)
    * Date of last update: 31/12/2016
    * Container: `cimendes/idba:31.12.2016-3`


### Benchmark datasets

* **Zymobiomics Community Standard**
Two commercially available mock communities containing 10 microbial species (ZymoBIOMICS Microbial Community Standards) 
were sequences by [Nicholls et al. 2019](https://academic.oup.com/gigascience/article/8/5/giz043/5486468). Shotgun 
sequencing of the **Even and Log communities** was performed with the same protocol, with the exception that the Log 
community was sequenced individually on 2 flowcell lanes and the Even community was instead sequenced on an Illumina 
MiSeq using 2×151 bp (paired-end) sequencing. They are available under accession numbers: 
* [ERR2984773](https://www.ebi.ac.uk/ena/data/view/ERR2984773) (even) 
* [ERR2935805](https://www.ebi.ac.uk/ena/data/view/ERR2935805) (log)

With the [M3S3 tool](http://medweb.bgu.ac.il/m3s3/), a simulation sample was obtained made up of the Zymobiomics 
microbial community standasrd species of **bacteria** (n=8).


### Assessing Assembly Success

The assembly performance is evaluated through mapping of the obtained assemblies to the [triple bacterial reference sequences](data/references/Zymos_Genomes_triple_chromosomes.fasta), 
with [minimap2](https://github.com/lh3/minimap2), using the [docker image](docker/minimap2/Dockerfile) 
`cimendes/minimap2:2.17-1`. 

To map the contigs to the triple references, the ` minimap2 -c -t 4 -r 10000 -g 10000 -x asm20 --eqx --secondary=no 
Zymos_Genomes_triple_chromosomes.fasta <in_file>.fasta <out_file>.paf` command was used. To map the read dataset to the 
assemblies we used the command `minimap2 -x sr --secondary=no <assembly> <forward_read> <reverse_read>`.
For each command, the resulting mapping files, in [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md) describes the mapping 
positions between each mapped contig and the reference, and between each read and the assemby, respectively.

The sets of python scripts used to obtain the evaluation metrics are available in the [scripts](scripts) folder, 
with the step by step analysis available as a [Jyputer Notebook](analysis/run_analysis.ipynb).

#### Main metrics implemented

* **Assembly Statistics**
This is mostly obtained with the [assembly_stats_global](scripts/assembly_stats_global.py) and the 
[assembly_mapping_stats_global.py](scripts/assembly_mapping_stats_global.py) scripts.

* **Accuracy of Assembly** 
This information is obtained through the [assembly_mapping_stats_per_ref](scripts/assembly_mapping_stats_per_ref.py) 
python script. This script produces a boxplot of the mapped contig size distribution for each assembler, with unmmaped
contigs showed as a scatter plot.

On the metrics used, [Rick et al. 2019](https://github.com/rrwick/Long-read-assembler-comparison) proposed the use of a 
triple reference to assess **chromosome contiguity** while benchmaking long-read genomic assemblers. This measure is the longest single 
alignment between the assembly and the reference. More information is available on 
[Rick's Assembly Benchmark GitHub page](https://github.com/rrwick/Long-read-assembler-comparison#assessing-chromosome-contiguity).

:warning:TODO:warning: New/other metrics implemented

* **Accuracy of Assembly** 
:warning:TODO:warning:

## Results
:warning:TODO:warning:
The results for each dataset are available in the [results](results) folder.

## Authors
* Inês Mendes <cimendes@medicina.ulisboa.pt> 
    * Instituto de Microbiologia, Instituto de Medicina Molecular, Faculdade de Medicina, Universidade de Lisboa, 
    Lisboa, Portugal; 
    * University of Groningen, University Medical Center Groningen, Department of Medical Microbiology and Infection 
    Prevention, Groningen, The Netherlands
* Rafael Mamede
    * Instituto de Microbiologia, Instituto de Medicina Molecular, Faculdade de Medicina, Universidade de Lisboa, 
    Lisboa, Portugal; 
* Yair Motro 
    * Department of Health System Management, School of Public Health, Faculty of Health Sciences, Ben-Gurion 
University of the Negev, Beer-Sheva, Israel.

## References
* Angers-Loustau, A., Petrillo, M., Bengtsson-Palme, J., Berendonk, T., Blais, B., Chan, K.-G., et al. (2018). The challenges of designing a benchmark strategy for bioinformatics pipelines in the identification of antimicrobial resistance determinants using next generation sequencing technologies. F1000Research 7, 459. doi:10.12688/f1000research.14509.1.
* N Couto, L Schuele, EC Raangs, MP Machado, ​CI Mendes​, TF Jesus, M Chlebowicz, S Rosema, M Ramirez, JA Carriço, IB Autenrieth, AW Friedrich, S Peter, JW Rossen. (2018) Critical steps in clinical shotgun metagenomics for the concomitant detection and typing of microbial pathogens. Scientific Reports 8: 13767. Doi: 10.1038/s41598-018-31873-w
* Gruening, B., Sallou, O., Moreno, P., Leprevost, F. D. V., Ménager, H., Søndergaard, D., et al. (2018). Recommendations for the packaging and containerizing of bioinformatics software. F1000Research 7, 742. doi:10.12688/f1000research.15140.1.
* Kurtzer, G. M., Sochat, V., and Bauer, M. W. (2017). Singularity: Scientific containers for mobility of compute. Plos One 12. doi:10.1371/journal.pone.0177459
* Li, D., Liu, C.-M., Luo, R., Sadakane, K., and Lam, T.-W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics 31, 1674–1676. doi:10.1093/bioinformatics/btv033.
* Samuel M Nicholls, Joshua C Quick, Shuiquan Tang, Nicholas J Loman, Ultra-deep, long-read nanopore sequencing of mock microbial community standards, GigaScience, Volume 8, Issue 5, May 2019, giz043, https://doi.org/10.1093/gigascience/giz043 
* Nurk, S., Meleshko, D., Korobeynikov, A., and Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome Research 27, 824–834. doi:10.1101/gr.213959.116.
* Olson, N. D., Treangen, T. J., Hill, C. M., Cepeda-Espinoza, V., Ghurye, J., Koren, S., et al. (2017). Metagenomic assembly through the lens of validation: recent advances in assessing and improving the quality of genomes assembled from metagenomes. Briefings in Bioinformatics. doi:10.1093/bib/bbx098.
* Sczyrba, A., Hofmann, P., Belmann, P., Koslicki, D., Janssen, S., Dröge, J., et al. (2017). Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software. Nature Methods 14, 1063–1071. doi:10.1038/nmeth.4458
* Teeling, H., and Glockner, F. O. (2012). Current opportunities and challenges in microbial metagenome analysis--a bioinformatic perspective. Briefings in Bioinformatics 13, 728–742. doi:10.1093/bib/bbs039.
* Tommaso, P. D., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., and Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319. doi:10.1038/nbt.3820
