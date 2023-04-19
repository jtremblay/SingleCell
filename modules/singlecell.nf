/*
    Single cell pipeline functions 
*/

process STAR {
    debug true
    publishDir "$params.DEFAULT.outdir/star/", mode: 'symlink'
    cpus params.star.cluster_cpus
    memory params.star.cluster_memory
    time params.star.cluster_time
    
    input:
        tuple val(sample_id), path(reads)

    output:
        //tuple val("$sample_id"), path("${sample_id}/.bam"), emit: bam
        path("${sample_id}_features.tsv"), emit: matrix
        path("${sample_id}_barcodes.tsv"), emit: features
        path("${sample_id}_matrix.mtx"), emit: barcodes
        path("${sample_id}_Features.stats"), emit: features_stats
        path("${sample_id}_Summary.csv"), emit: summary

    script:
        """
        module load ${params.modules.star} && \\
        STAR \\
          --genomeDir ${params.DEFAULT.ref_genome} \\
          --soloType Droplet \\
          --readFilesIn ${reads[1]} ${reads[0]} \\
          --soloCBwhitelist ${params.DEFAULT.whitelist} \\
          --soloUMIlen ${params.DEFAULT.umi_length} \\
          --runThreadN ${params.star.num_threads} && \\
        mv Solo.out/Gene/filtered/features.tsv ${sample_id}_features.tsv && \\
        mv Solo.out/Gene/filtered/barcodes.tsv ${sample_id}_barcodes.tsv && \\
        mv Solo.out/Gene/filtered/matrix.mtx ${sample_id}_matrix.mtx && \\
        mv Solo.out/Gene/Features.stats ${sample_id}_Features.stats && \\
        mv Solo.out/Gene/Summary.csv ${sample_id}_Summary.csv

        """
}

