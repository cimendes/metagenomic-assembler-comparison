#!/usr/bin/env nextflow

import Helper
import CollectInitialMetadata

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

def infoMap = [:]
if (params.containsKey("fastq")){
    infoMap.put("fastq", file(params.fastq).size())
}

//Help.start_info(infoMap, "$workflow.start", "$workflow.profile", version)
CollectInitialMetadata.print_metadata(workflow)

// MAIN PARAMETERS
//      Fastq
if (params.fastq instanceof Boolean){exit 1, "'fastq' must be a path pattern. Provide value:'$params.fastq'"}
if (!params.fastq){ exit 1, "'fastq' parameter missing"}
// size: -1 -> allows for single and paired-end files to be passed through. Change if necessary
IN_fastq_raw = Channel.fromFilePairs(params.fastq, size: -1).ifEmpty { exit 1, "No fastq files provided with pattern:'${params.fastq}'" }

//      Clear
clear = params.clearInput ? "true" : "false"
checkpointClear = Channel.value(clear)


// SET CHANNELS FOR ASSEMBLERS
IN_fastq_raw.into{ IN_BCALM2; IN_GATB_MINIA; IN_MEGAHIT; IN_METASPADES; IN_UNICYCLER; IN_IDBA; IN_SPADES; IN_SKESA; IN_VELOUR; IN_VELVETOPTIMIZER; IN_PANDASEQ; IN_BBAP }


// BCALM 2

if ( !params.bcalmKmerSize.toString().isNumber() ){
    exit 1, "'bcalmKmerSize' parameter must be a number. Provided value: '${params.bcalmKmes%rSize}'"
}

process BCALM2 {

    tag { sample_id }
    publishDir "results/bcalm2/"

    input:
    set sample_id, file(fastq) from IN_BCALM2
    val KmerSize from Channel.value(params.bcalmKmerSize)

    output:
    file "*_BCALM2.fasta"

    script:
    """
    ls -1 $fastq  > list_reads
    bcalm -in list_reads -out ${sample_id} -kmer-size $KmerSize

    # workdir cleanup
    rm list_reads
    rm *.h5
    rm *.glue.*

    mv ${sample_id}.unitigs.fa  ${sample_id}_BCALM2.fasta
    """
}


// GATB MINIA

IN_GATB_kmers = Channel.value(params.gatbkmer)
IN_GATB_besst_iter = Channel.value(params.gatb_besst_iter)
GATB_error_correction = params.GATB_error_correction ? "true" : "false"
IN_error_correction = Channel.value(GATB_error_correction)

process GATBMINIAPIPELINE {

    tag { sample_id }
    publishDir 'results/GATBMiniaPipeline/', pattern: '*.fasta'

    input:
    set sample_id, file(fastq_pair) from IN_GATB_MINIA
    val kmer_list from IN_GATB_kmers
    val do_error_correction from GATB_error_correction
    val besst_iter from IN_GATB_besst_iter

    output:
    file('*.fasta')

    script:
    """
    if [ $do_error_correction ];
    then
        gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} -o ${sample_id}_GATBMiniaPipeline
    else
        gatb -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} --kmer-sizes ${kmer_list} -o ${sample_id}_GATBMiniaPipeline --no-error-correction
    fi

    link=$(readlink ${sample_id}_GATBMiniaPipeline.fasta) && rm ${sample_id}_GATBMiniaPipeline.fasta && mv $link ${sample_id}_GATBMiniaPipeline.fasta

    # rm temp dirs
    rm -r ${sample_id}_GATBMiniaPipeline.lib* ${sample_id}_GATBMiniaPipeline_besst *.unitigs* *contigs.fa *.h5
    rm *list_reads*

    """
}


// MEGAHIT

IN_megahit_kmers = Channel.value(params.megahitKmers)

process MEGAHIT {

    tag { sample_id }
    publishDir 'results/MEGAHIT/', pattern: '*_megahit*.fasta'

    input:
    set sample_id, file(fastq_pair) from IN_MEGAHIT
    val kmers from IN_megahit_kmers

    output:
    set sample_id, file('*_megahit.fasta')

    script:
    """
    /NGStools/megahit/bin/megahit --num-cpu-threads $task.cpus -o megahit --k-list $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}
    mv megahit/final.contigs.fa ${sample_id}_megahit.fasta
    rm -r megahit
    """

}


// METASPADES

if ( params.metaspadesKmers.toString().split(" ").size() <= 1 ){
    if (params.metaspadesKmers.toString() != 'auto'){
        exit 1, "'metaspadesKmers' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.metaspadesKmers}"
    }
}
IN_metaspades_kmers = Channel.value(params.metaspadesKmers)

process METASPADES {

    tag { sample_id }
    publishDir 'results/MetaSPAdes/'

    input:
    set sample_id, file(fastq_pair) from IN_METASPADES
    val kmers from IN_metaspades_kmers

    output:
    set sample_id, file('*_metaspades.fasta')
    file('*_metaspades.fastg')

    script:
    """
    metaspades.py --only-assembler --threads $task.cpus -k $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o metaspades
    mv metaspades/contigs.fasta ${sample_id}_metaspades.fasta
    mv metaspades/assembly_graph.fastg ${sample_id}_metaspades.fastg
    rm -r metaspades
    """
}


// UNICYCLER
process UNICYCLER {

    tag { sample_id }
    publishDir 'results/unicycler'

    input:
    set sample_id, file(fastq_pair) from IN_UNICYCLER

    output:
    set sample_id, file('*_unicycler.*')
    file('*_unicycler.gfa')

    script:
    """
    unicycler -t $task.cpus -o . --no_correct --no_pilon -1 ${fastq_pair[0]} -2 ${fastq_pair[1]}
    mv assembly.fasta ${sample_id}_unicycler.fasta
    mv assembly.gfa ${sample_id}_unicycler.gfa
    rm *best_spades_graph* *overlaps_removed* *bridges_applied* *final_clean*
    """
}


// SPADES
if ( params.spadesKmers.toString().split(" ").size() <= 1 ){
    if (params.spadesKmers.toString() != 'auto'){
        exit 1, "'spadesKmers' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.spadesKmers}"
    }
}
IN_spades_kmers = Channel.value(params.spadesKmers)

process SPADES {

    tag { sample_id }
    publishDir 'results/SPAdes/', pattern: '*.fasta'

    input:
    set sample_id, file(fastq_pair) from IN_SPADES
    val kmers from IN_spades_kmers

    output:
    set sample_id, file('*_spades.fasta')
    file('*_spades.fastg')

    script:
    """
    spades.py --only-assembler --threads $task.cpus -k $kmers -1 ${fastq_pair[0]} -2 ${fastq_pair[1]} -o spades
    mv spades/contigs.fasta ${sample_id}_spades.fasta
    mv spades/assembly_graph.fastg ${sample_id}_spades.fastg
    rm -r spades
    """
}

// SKESA
process SKESA {

    tag { sample_id }
    publishDir 'results/skesa/'

    input:
    set sample_id, file(fastq_pair) from IN_SKESA

    output:
    set sample_id, file('*_skesa.fasta')

    script:
    """
    skesa --cores $task.cpus --memory $task.memory --use_paired_ends --contigs_out ${sample_id}_skesa.fasta --gfa ${sample_id}_skesa.fastg --fastq ${fastq_pair[0]} ${fastq_pair[1]}
    """
}


// PANDASEQ
process PANDASEQ {

    tag { sample_id }
    publishDir 'results/pandaseq/$sample_id', pattern: 'assembly.fasta'

    input:
    set sample_id, file(fastq_pair) from IN_PANDASEQ

    output:
    set sample_id, file('*pandaseq.fasta')

    script:
    """
    cp -r /NGStools/pandaseq pandaseq/

    ./pandaseq/pandaseq -T $task.cpus -w ${sample_id}_pandaseq.fasta -f ${fastq_pair[0]} -r ${fastq_pair[1]}
    """
}


// VelvetOptimizer
process VELVETOPTIMIZER {
    tag { sample_id }
    publishDir 'results/velvet_optimiser/', pattern: '*fasta'

    input:
    set sample_id, file(fastq_pair) from IN_VELVETOPTIMIZER

    output:
    file('*.fasta')

    script:
    """
    VelvetOptimiser.pl -v -s $params.velvetoptimizer_hashs -e $params.velvetoptimizer_hashe -t $task.cpus \
    -f '-shortPaired -fastq.gz -separate ${fastq_pair[0]} ${fastq_pair[1]}'
    """

}


// IDBA

process reformat_IDBA {
    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from IN_IDBA

    output:
    set sample_id, file('*.fasta') into REFORMAT_IDBA

    script:
    "reformat.sh in=${fastq_pair[0]} in2=${fastq_pair[1]} out=${sample_id}_reads.fasta"
}

process IDBA {

    tag { sample_id }
    publishDir 'results/idba/', pattern: '*fasta'

    input:
    set sample_id, file(fasta_reads_single) from  REFORMAT_IDBA

    output:
    set sample_id, file('*_idba_contig.fa')

    script:
    """
    idba_ud -l ${fasta_reads_single} --num_threads $task.cpus -o .
    mv contig.fa ${sample_id}_idba_contig.fa

    """

}


// VELOUR
IN_kmers_velour = Channel.value(params.velourKmer)

process reformat_VELOUR {
    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from IN_VELOUR

    output:
    set sample_id, file('*.fasta') into REFORMAT_VELOUR

    script:
    "reformat.sh in=${fastq_pair[0]} in2=${fastq_pair[1]} out=${sample_id}_1.fasta out2=${sample_id}_2.fasta"
}

process VELOUR {
    tag { sample_id }
    publishDir 'results/assembly/Velour_{{ pid }}/', pattern: 'out/*.fasta', mode: 'copy'

    input:
    set sample_id, file(fasta_reads_pair) from REFORMAT_VELOUR
    val kmer from IN_kmers_velour

    output:
    file('out/*.fasta')

    script:
    """
    velour out ${kmer} ${fasta_reads_pair}
    """
}
