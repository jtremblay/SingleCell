/*
    Shotgun metagenomics functions 
*/

process TRIMMOMATIC {
    publishDir "$params.DEFAULT.outdir/qced_reads/", mode: 'symlink'
    debug true
    cpus params.trimmomatic.cluster_cpus
    memory params.trimmomatic.cluster_memory
    time params.trimmomatic.cluster_time


    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_paired_R1.fastq.gz"), path("${sample_id}_paired_R2.fastq.gz"), emit: reads 
        path("${sample_id}_trimmomatic_log.txt"), emit: log
        path("${sample_id}_trimmomatic_stats.tsv"), emit: stats

    script:
        command = """
        module load ${params.modules.trimmomatic} $params.modules.java && \\
        java -XX:ParallelGCThreads=$params.trimmomatic.threads -Xmx2G -jar \$TRIMMOMATIC_JAR PE \\
            -threads ${params.trimmomatic.threads} \\
            -phred${params.trimmomatic.quality_offset} \\
            ${reads[0]} ${reads[1]} \\
            ${sample_id}_paired_R1.fastq.gz ${sample_id}_single_R1.fastq.gz ${sample_id}_paired_R2.fastq.gz ${sample_id}_single_R2.fastq.gz \\
            ILLUMINACLIP:${params.trimmomatic.adapter_fasta}${params.trimmomatic.illumina_clip_settings} \\
            TRAILING:${params.trimmomatic.trailing_min_quality} \\
            SLIDINGWINDOW:${params.trimmomatic.sliding_window1}:${params.trimmomatic.sliding_window2} \\
            MINLEN:${params.trimmomatic.min_length} \\
            HEADCROP:${params.trimmomatic.headcrop}"""
        if(params.trimmomatic.crop){
            command += """ CROP:${params.trimmomatic.crop}"""
        }
        command += """ 2> ${sample_id}_trimmomatic_log.txt && \\
        grep ^Input ${sample_id}_trimmomatic_log.txt | \\
        perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*\$/Raw Fragments,\\1\\nFragment Surviving,\\2\\nSingle Surviving,\\3/' >  ${sample_id}_trimmomatic_stats.tsv
        """
        return command

}

process BBDUK {
    publishDir "$params.DEFAULT.outdir/qced_reads/", mode: 'symlink'
    debug true
    cpus params.bbduk.cluster_cpus
    memory params.bbduk.cluster_memory
    time params.bbduk.cluster_time

    input:
        tuple val(sample_id), path(infile_R1), path(infile_R2)

    output:
        tuple val(sample_id), path("${sample_id}_paired_ncontam_R1.fastq.gz"), path("${sample_id}_paired_ncontam_R2.fastq.gz"), emit: reads
        tuple val(sample_id), path("${sample_id}_paired_contam_R1.fastq.gz"), path("${sample_id}_paired_contam_R2.fastq.gz"), emit: bad_reads
        path("${sample_id}_bbduk_log.txt"), emit: log

    script:
        """
        module load ${params.modules.bbmap} ${params.modules.java} && \\
        bbduk.sh \\
            in=${infile_R1} \\
            in2=${infile_R2} \\
            out=${sample_id}_paired_ncontam_R1.fastq.gz \\
            out2=${sample_id}_paired_ncontam_R2.fastq.gz \\
            outm=${sample_id}_paired_contam_R1.fastq.gz \\
            outm2=${sample_id}_paired_contam_R2.fastq.gz \\
            stats=${sample_id}_bbduk_log.txt \\
            k=${params.bbduk.k} \\
            minkmerhits=${params.bbduk.c} \\
            ref=${params.DEFAULT.contaminants} \\
            overwrite=true \\
            threads=1
        """
}

process BBMAP_SUBTRACT {
    publishDir "$params.DEFAULT.outdir/qced_reads/", mode: 'symlink'
    debug true
    cpus params.bbmap_sub.cluster_cpus
    memory params.bbmap_sub.cluster_memory
    time params.bbmap_sub.cluster_time

    input:
        tuple val(sample_id), path(infile_R1), path(infile_R2)

    output:
        tuple val(sample_id), path("${sample_id}_paired_ncontam_sub_R1.fastq.gz"), path("${sample_id}_paired_ncontam_sub_R2.fastq.gz"), emit: reads
        path("${sample_id}_bbduksub_log.txt"), emit: log

    script:
        """
        module load ${params.modules.bbmap} ${params.modules.java} && \\
        bbmap.sh \\
            -threads=${params.bbmap_sub.num_threads} \\
            -Xmx${params.bbmap_sub.ram} \\
            in=${infile_R1} in2=${infile_R2} \\
            outu=${sample_id}_paired_ncontam_sub_R1.fastq.gz outu2=${sample_id}_paired_ncontam_sub_R2.fastq.gz \\
            minid=${params.bbmap_sub.min_id} \\
            sam=1.3 nhtag=t mdtag=t xmtag=t amtag=t nmtag=t xstag=us \\
            maxindel=3 bwr=0.16 bw=12 fast=t overwrite=t minhits=2 \\
            path=${params.DEFAULT.ref_genome_to_subtract} \\
            statsfile=${sample_id}_bbduksub_log.txt
        """
}

process MEGAHIT {
    publishDir "$params.DEFAULT.outdir/assembly/", mode: 'symlink'
    debug true
    cpus params.megahit.cluster_cpus
    memory params.megahit.cluster_memory
    time params.megahit.cluster_time
    
    input:
        val(sample_id)
        path(reads1)
        path(reads2)

    output:
        path("megahit/contigs.fna"), emit: assembly
        path("megahit/log"), emit: log
        path("megahit/assembly_stats.txt"), emit: stats

    script:
        def input_R1 = "-1 " + reads1.join(",")
        def input_R2 = "-2 " + reads2.join(",")
        def requested_memory_in_bytes = params.megahit.cluster_memory.toBytes()

        """
        module load ${params.modules.megahit} ${params.modules.tools} && \\
        rm -rf megahit -rf && \\
        megahit -t ${params.megahit.num_threads} --k-min ${params.megahit.kmin} \\
            --k-max ${params.megahit.kmax} --k-step ${params.megahit.kstep} \\
            --min-contig-len ${params.megahit.min_contig_length} \\
            ${input_R1} \\
            ${input_R2} \\
            --memory ${requested_memory_in_bytes} \\
            --out-dir megahit && \\
            compileAssemblyResultsSingle.pl --infile megahit/final.contigs.fa > megahit/assembly_stats.txt && mv megahit/final.contigs.fa megahit/contigs.fna
        """
}

process METASPADES {
    publishDir "$params.DEFAULT.outdir/assembly/", mode: 'symlink'
    debug true
    
    input:
        val(sample_id)
        path(reads1)
        path(reads2)

    output:
        path("spades/scaffolds_gt${params.spades.length}.fasta"), emit: assembly
        path("spades/assembly_stats.txt"), emit: stats

    script:
        def input_R1 = " -1 " + reads1.join(" -1 ")
        def input_R2 = " -2 " + reads2.join(" -2 ")
        """
        module load ${params.modules.spades} ${params.modules.tools} && \\
        rm spades_outdir -rf && \\
        spades.py --meta -t {params.spades.num_threads} -k {params.spades.kmers} \\
            ${input_R1} \\
            ${input_R2} \\
            -o spades_outdir -m ${params.spades.memory} && \\
        filterFastaByLength.pl \\
            --infile spades_outdir/scaffolds.fasta \\
            --length ${params.spades.length} > \\
            spades/scaffolds_gt${params.spades.length}.fasta && \\
        compileAssemblyResultsSingle.pl --infile spades/scaffolds_gt${params.spades.length}.fasta > spades/assembly_stats.txt
        """
}

process PRODIGAL {
    publishDir "$params.DEFAULT.outdir/gene_prediction/", mode: 'symlink'
    debug true
    
    input:
        path(infile)

    output:
        path("contigs_renamed.fna"), emit: predicted_genes_fna
        path("contigs_renamed.faa"), emit: predicted_genes_faa
        path("contigs_renamed.gff"), emit: predicted_genes_gff
        path("contigs_raw.faa"), emit: predicted_genes_raw_faa

    script:
        """
        module load ${params.modules.prodigal} ${params.modules.tools} && \\
        prodigal -i ${infile} -f gff -p meta \\
            -o contigs_raw.gff \\
            -a contigs_raw.faa \\
            -d contigs_raw.fna && \\
        convertProdigalNames.pl \\
            --gff contigs_raw.gff \\
            --fna contigs_raw.fna \\
            --faa contigs_raw.faa \\
            --renamed_gff contigs_renamed.gff \\
            --renamed_faa contigs_renamed.faa \\
            > contigs_renamed.fna
        """
}

process EXONERATE_CONTIGS {
    debug true
    publishDir "$params.DEFAULT.outdir/", mode: 'symlink'
    
    input:
        path(infile)

    output:
        path("estimated_number_of_chunks_contigs.txt"), emit: chunks
        path("exonerate.done"), emit: done
        path("fasta_chunks"), emit: outdir

    script:
        """
        module load ${params.modules.exonerate} ${params.modules.tools} && \\
        estimateChunkFileSize.pl \\
            --infile ${infile} \\
            --targeted_chunk_size ${params.exonerate.targeted_chunk_file_size_contigs} > \\
            estimated_number_of_chunks_contigs.txt && \\
        rm -rf fasta_chunks && mkdir -p fasta_chunks && \\
        fastasplit -f ${infile} \\
            -o fasta_chunks \\
            -c `cat estimated_number_of_chunks_contigs.txt` && \\
        doublecheckExonerate.pl \\
            --indir fasta_chunks \\
            --prefix contigs.fna_chunk_ \\
            --no_chunks `cat estimated_number_of_chunks_contigs.txt` && \\
          touch exonerate.done && touch fasta_chunks
        """
}

process EXONERATE_GENES {
    debug true
    publishDir "$params.DEFAULT.outdir/gene_prediction", mode: 'symlink'
    
    input:
        path(infile)

    output:
        path("estimated_number_of_chunks_genes.txt"), emit: chunks
        path("exonerate.done"), emit: done
        path("fasta_chunks"), emit: outdir
        env(NUMCHUNKS), emit: num_chunks
        path("fasta_chunks/*"), emit: outfiles // Here the fact that we emit multiple files
                                               // (using the * wildcard) means that we will
                                               // generate a multi-file channel. Otherwise,
                                               // it would have been a single file channel.

    script:
        """
        module load ${params.modules.exonerate} ${params.modules.tools} && \\
        estimateChunkFileSize.pl \\
            --infile ${infile} \\
            --targeted_chunk_size ${params.exonerate.targeted_chunk_file_size_genes} > \\
            estimated_number_of_chunks_genes.txt && \\
        rm -rf fasta_chunks && mkdir -p fasta_chunks && \\
        fastasplit -f ${infile} \\
            -o fasta_chunks \\
            -c `cat estimated_number_of_chunks_genes.txt` && \\
        doublecheckExonerate.pl \\
            --indir fasta_chunks \\
            --prefix contigs_renamed.faa_chunk_ \\
            --no_chunks `cat estimated_number_of_chunks_genes.txt` && \\
          touch exonerate.done && touch fasta_chunks && \\
          NUMCHUNKS=\$(cat estimated_number_of_chunks_genes.txt)
        """
}

process BEDFILE_CONTIGS {
    debug true
    publishDir "$params.DEFAULT.outdir/assembly/megahit", mode: 'symlink'
    
    input:
        path(infile_contigs)

    output:
        path("contigs.bed"), emit: bed_contigs

    script:
        """
        module load ${params.modules.tools} && \\
        fastaToBed.pl --fasta ${infile_contigs} > contigs.bed
        """
}

process BEDFILE_GENES {
    debug true
    publishDir "$params.DEFAULT.outdir/gene_prediction", mode: 'symlink'
    
    input:
        path(infile_gff)

    output:
        path("contigs_genes.bed"), emit: bed_genes

    script:
        """
        module load ${params.modules.tools} && \\
        gffToBed.pl --infile ${infile_gff} > contigs_genes.bed
        """
}

process MAKE_BWA_INDEX {
    debug true
    publishDir "$params.DEFAULT.outdir/assembly", mode: 'symlink'
    cpus params.make_index.cluster_cpus
    memory params.make_index.cluster_memory
    time params.make_index.cluster_time
    
    input:
        path(infile)

    output:
        path("contigs.fna.amb"), emit: amb
        path("contigs.fna.ann"), emit: ann
        path("contigs.fna.bwt"), emit: bwt
        path("contigs.fna.fai"), emit: fai
        path("contigs.fna.pac"), emit: pac
        path("contigs.fna.sa"),  emit: sa

    script:
        """
        module load ${params.modules.bwa} && \\
        bwa index ${infile} && \\
        touch contigs.fna.amb && \\
        touch contigs.fna.ann && \\
        touch contigs.fna.bwt && \\
        touch contigs.fna.fai && \\
        touch contigs.fna.pac && \\
        touch contigs.fna.sa
        """
}

process BWAMEM_PE {
    debug true
    publishDir "$params.DEFAULT.outdir/contig_abundance/bams", mode: 'symlink'
    cpus params.bwa.cluster_cpus
    memory params.bwa.cluster_memory
    time params.bwa.cluster_time
    
    input:
        path(amb)
        path(ann)
        path(bwt)
        path(fai)
        path(pac)
        path(sa)
        path(infile_fasta)
        tuple val(sample_id), path(infile_R1), path(infile_R2)

    output:
        tuple val("$sample_id"), path("${sample_id}.bam"), emit: bam
        path("${sample_id}.bam.bai"), emit: bam_index

    script:
        """
        module load ${params.modules.bwa} ${params.modules.samtools} && \\
        bwa mem -M \\
          -t ${params.bwa.num_threads} \\
          ${infile_fasta} \\
          ${infile_R1} \\
          ${infile_R2} \\
          | samtools view -Sbh -F 0x100 - > ${sample_id}.bam.tmp && \\
          samtools sort -@ ${params.bwa.num_threads} -m ${params.samtools.mem_per_thread} ${sample_id}.bam.tmp -o ${sample_id}.bam && \\
          rm ${sample_id}.bam.tmp && \\
          samtools index ${sample_id}.bam && touch ${sample_id}.bam.bai
        """
}

process BEDTOOLS_COV_CONTIGS {
    debug true
    publishDir "$params.DEFAULT.outdir/contig_abundance/cov", mode: 'symlink'
    cpus params.bedtools.cluster_cpus
    memory params.bedtools.cluster_memory
    time params.bedtools.cluster_time
    
    input:
        tuple val(sample_id), path(bam)
        path(bed_contigs)

    output:
        tuple val("$sample_id"), path("${sample_id}.cov"), emit: cov

    script:
        """
        module load ${params.modules.bedtools} ${params.modules.samtools} && \\
        samtools view -b -f 0x2 ${bam} > ${bam}.tmp && \\
        coverageBed -abam ${bam}.tmp \\
            -b ${bed_contigs} \\
            -counts \\
            > ${sample_id}.cov && \\
        rm -f ${bam}.tmp
        """
}

process BEDTOOLS_COV_GENES {
    debug true
    publishDir "$params.DEFAULT.outdir/gene_abundance/cov", mode: 'symlink'
    cpus params.bedtools.cluster_cpus
    memory params.bedtools.cluster_memory
    time params.bedtools.cluster_time
    
    input:
        tuple val(sample_id), path(bam)
        path(bed_genes)

    output:
        tuple val(sample_id), path("${sample_id}.cov"), emit: cov

    script:
        """
        module load ${params.modules.bedtools} ${params.modules.samtools} && \\
        samtools view -b -f 0x2 ${bam} > ${bam}.tmp && \\
        coverageBed -abam ${bam}.tmp \\
            -b ${bed_genes} \\
            -counts \\
            > ${sample_id}.cov && \\
        rm -f ${bam}.tmp
        """
}

/*
    Can't use the same process twice?
    Process 'MERGE_COV' has been already used -- If you need to reuse the same component, include it with a different name or include it in a different workflow context
    Anyways, I'll leave it here if it becomes ever supported.
*/
process MERGE_COV {
    debug true
    publishDir "$params.DEFAULT.outdir/${type}_abundance/", mode: 'symlink'
    cpus params.merge_abundance.cluster_cpus
    memory params.merge_abundance.cluster_memory
    time params.merge_abundance.cluster_time
    
    input:
        path(cov_files)
        val(type)

    output:
        path("${type}_abundance.tsv"), emit: contig_abundance_matrix

    script:
        cov_files_comma = cov_files.toString().replaceAll(/ /, ",")
        """
        module load ${params.modules.perl} ${params.modules.tools} && \\
        mergeAbundance.pl \\
            --type contigs \\
            --infiles ${cov_files_comma} \\
            > ${type}_abundance.tsv
        """
}

process MERGE_COV_CONTIGS {
    debug true
    publishDir "$params.DEFAULT.outdir/contig_abundance/", mode: 'symlink'
    cpus params.merge_abundance.cluster_cpus
    memory params.merge_abundance.cluster_memory
    time params.merge_abundance.cluster_time
    
    input:
        path(cov_files)

    output:
        path("contig_abundance.tsv"), emit: contig_abundance_matrix

    script:
        cov_files_comma = cov_files.toString().replaceAll(/ /, ",")
        """
        module load ${params.modules.perl} ${params.modules.tools} && \\
        mergeAbundance.pl \\
            --type contigs \\
            --infiles ${cov_files_comma} \\
            > contig_abundance.tsv
        """
}

process MERGE_COV_GENES {
    debug true
    publishDir "$params.DEFAULT.outdir/gene_abundance/", mode: 'symlink'
    cpus params.merge_abundance.cluster_cpus
    memory params.merge_abundance.cluster_memory
    time params.merge_abundance.cluster_time
    
    input:
        path(cov_files)

    output:
        path("gene_abundance.tsv"), emit: gene_abundance_matrix

    script:
        cov_files_comma = cov_files.toString().replaceAll(/ /, ",")
        """
        module load ${params.modules.perl} ${params.modules.tools} && \\
        mergeAbundance.pl \\
            --type genes \\
            --infiles ${cov_files_comma} \\
            > gene_abundance.tsv
        """
}

process DIAMOND_BLASTP_NR {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/diamond_blastp_nr/", mode: 'symlink'
    cpus params.diamond_blastp.cluster_cpus
    memory params.diamond_blastp.cluster_memory
    time params.diamond_blastp.cluster_time

    input:
        /* 
            Here 'each' (instead of 'path') will launch each process in parallel. For instance, if your have
            a channel of 5 path in the infile channel (here the infile channel is :
            EXONERATE_GENES.out.outfiles, 5 diamond-blastp processes will be launch in
            parallel. Finally, I had to use tuple and not each.
        */
        tuple val(prefix), path(infile)

    output:
        path("${prefix}.tsv"), emit: diamond_blastp_outfiles 
                                                             
                                                             
                                                             

    script:
        """
        module load ${params.modules.diamond} && \\

        if [[ -s ${infile} ]] ; then
            diamond blastp --quiet \\
                -d ${params.diamond_blastp.db_nr} \\
                -q ${infile} \\
                -o ${prefix}.tsv \\
                -k 10 \\
                -e ${params.diamond_blastp.evalue} \\
                -p ${params.diamond_blastp.num_threads}
        else
            touch ${prefix}.tsv
        fi ;
        """
}

process HMMSEARCH_KEGG {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/hmmsearch_kegg/", mode: 'symlink'
    cpus params.hmmsearch.cluster_cpus
    memory params.hmmsearch.cluster_memory
    time params.hmmsearch.cluster_time

    input:
        tuple val(prefix), path(infile)

    output:
        path("${prefix}_hmmsearch_kegg_domtblout.tsv"), emit: domtblout_outfiles
        path("${prefix}_hmmsearch_kegg_tblout.tsv"), emit: tblout_outfiles

    script:
        """
        module load ${params.modules.hmmer} && \\
        if [[ -s ${infile} ]] ; then
            hmmsearch \\
                --tblout ${prefix}_hmmsearch_kegg_tblout.tsv \\
                --domtblout ${prefix}_hmmsearch_kegg_domtblout.tsv \\
                -E ${params.hmmsearch.evalue} \\
                --cpu ${params.hmmsearch.num_threads} \\
                ${params.kegg.KO_profiles} \\
                ${infile} > /dev/null
        else
            touch ${prefix}_hmmsearch_kegg_tblout.tsv
            touch ${prefix}_hmmsearch_kegg_domtblout.tsv
        fi ;
        """
}

process HMMSEARCH_PFAM {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/hmmsearch_pfam/", mode: 'symlink'
    cpus params.hmmsearch.cluster_cpus
    memory params.hmmsearch.cluster_memory
    time params.hmmsearch.cluster_time

    input:
        tuple val(prefix), path(infile)

    output:
        path("${prefix}_hmmsearch_pfam_domtblout.tsv"), emit: domtblout_outfiles
        path("${prefix}_hmmsearch_pfam_tblout.tsv"), emit: tblout_outfiles

    script:
        """
        module load ${params.modules.hmmer} && \\
        if [[ -s ${infile} ]] ; then
            hmmsearch \\
                --tblout ${prefix}_hmmsearch_pfam_tblout.tsv \\
                --domtblout ${prefix}_hmmsearch_pfam_domtblout.tsv \\
                -E ${params.hmmsearch.evalue} \\
                --cpu ${params.hmmsearch.num_threads} \\
                ${params.pfam.db} \\
                ${infile} > /dev/null
        else
            touch ${prefix}_hmmsearch_pfam_tblout.tsv
            touch ${prefix}_hmmsearch_pfam_domtblout.tsv
        fi ;
        """
}

process RPSBLAST_COG {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/rpsblast_cog/", mode: 'symlink'
    cpus params.rpsblast.cluster_cpus
    memory params.rpsblast.cluster_memory
    time params.rpsblast.cluster_time

    input:
        /*
            In the end, a tuple and not two each channels had to be used in here.
        */
        tuple val(prefix), path(infile)

    output:
        path("${prefix}_rpsblast_cog.tsv"), emit: outfiles

    script:
        """
        module load ${params.modules.blast} && \\
        if [[ -s ${infile} ]] ; then
            rpsblast \\
                -db ${params.cog.db} \\
                -query ${infile} \\
                -out ${prefix}_rpsblast_cog.tsv \\
                -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \\
                -max_target_seqs 1 \\
                -evalue 1e-10 \\
                -num_threads 2
        else
            touch ${prefix}_rpsblast_cog.tsv
        fi ;
        """
}

process COG_OVERREP {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/", mode: 'symlink'
    cpus params.cog_overrep.cluster_cpus
    memory params.cog_overrep.cluster_memory
    time params.cog_overrep.cluster_time

    input:
        path(rpsblast_cog)
        path(gene_abundance)

    output:
        path("COG_abundance.tsv"), emit: COG_abundance

    script:
        """
        module load ${params.modules.tools} ${params.modules.python} && \\
        getCOG.py \\
            --infile-blastp ${rpsblast_cog} \\
            --infile-gene-abundance ${gene_abundance} \\
            > COG_abundance.tsv
        """
}

process RPSBLAST_KOG {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/rpsblast_kog/", mode: 'symlink'
    cpus params.rpsblast.cluster_cpus
    memory params.rpsblast.cluster_memory
    time params.rpsblast.cluster_time

    input:
        each prefix
        each infile

    output:
        path("${prefix}_rpsblast_kog.tsv"), emit: outfile

    script:
        """
        module load ${params.modules.blast} && \\
        if [[ -s ${infile} ]] ; then
            rpsblast \\
                -db ${params.kog.db} \\
                -query ${infile} \\
                -out ${prefix}_rpsblast_kog.tsv \\
                -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \\
                -max_target_seqs 1 \\
                -evalue 1e-10 \\
                -num_threads 2
        else
            touch ${prefix}_rpsblast_kog.tsv
        fi ;
        """
}

process MERGE_DIAMOND_BLASTP_NR {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/", mode: 'symlink'
    cpus params.merge.cluster_cpus
    memory params.merge.cluster_memory
    time params.merge.cluster_time

    input:
        path(infiles)

    output:
        path("diamond_blastp_nr.tsv"), emit: tsv

    script:
        """
        cat ${infiles} > diamond_blastp_nr.tsv
        """
}

process MERGE_COG {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/", mode: 'symlink'
    cpus params.merge.cluster_cpus
    memory params.merge.cluster_memory
    time params.merge.cluster_time

    input:
        path(infiles)

    output:
        path("rpsblast_cog.tsv"), emit: tsv

    script:
        """
        cat ${infiles} > rpsblast_cog.tsv
        """
}

process MERGE_KEGG {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/", mode: 'symlink'
    cpus params.merge.cluster_cpus
    memory params.merge.cluster_memory
    time params.merge.cluster_time

    input:
        path(infile_tblout)
        path(infile_domtblout)

    output:
        path("hmmsearch_kofam_tblout.tsv"), emit: tblout
        path("hmmsearch_kofam_domtblout.tsv"), emit: domtblout

    script:
        """
        cat ${infile_tblout} > hmmsearch_kofam_tblout.tsv.tmp && \\
        cat ${infile_domtblout} > hmmsearch_kofam_domtblout.tsv.tmp && \\
        awk 'BEGIN{OFS=FS=\" \"} NR<=3{print}; NR>3{tmp=\$1; \$1=\$3; \$3=tmp; tmp=\$2; \$2=\$4; \$4=tmp; print}' hmmsearch_kofam_tblout.tsv.tmp > hmmsearch_kofam_tblout.tsv && \\
        awk 'BEGIN{OFS=FS=\" \"} NR<=3{print}; NR>3{tmp=\$1; \$1=\$4; \$4=tmp; tmp=\$2; \$2=\$5; \$5=tmp; print}' hmmsearch_kofam_domtblout.tsv.tmp > hmmsearch_kofam_domtblout.tsv
        """
}

process MERGE_PFAM {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/", mode: 'symlink'
    cpus params.merge.cluster_cpus
    memory params.merge.cluster_memory
    time params.merge.cluster_time

    input:
        path(infile_tblout)
        path(infile_domtblout)

    output:
        path("hmmsearch_pfam_tblout.tsv"), emit: tblout
        path("hmmsearch_pfam_domtblout.tsv"), emit: domtblout

    script:
        """
        cat ${infile_tblout} > hmmsearch_pfam_tblout.tsv.tmp && \\
        cat ${infile_domtblout} > hmmsearch_pfam_domtblout.tsv.tmp && \\
        awk 'BEGIN{OFS=FS=\" \"} NR<=3{print}; NR>3{tmp=\$1; \$1=\$3; \$3=tmp; tmp=\$2; \$2=\$4; \$4=tmp; print}' hmmsearch_pfam_tblout.tsv.tmp > hmmsearch_pfam_tblout.tsv && \\
        awk 'BEGIN{OFS=FS=\" \"} NR<=3{print}; NR>3{tmp=\$1; \$1=\$4; \$4=tmp; tmp=\$2; \$2=\$5; \$5=tmp; print}' hmmsearch_pfam_domtblout.tsv.tmp > hmmsearch_pfam_domtblout.tsv
        """
}


process PARSE_KEGG {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/", mode: 'symlink'
    cpus params.parse_kofam.cluster_cpus
    memory params.parse_kofam.cluster_memory
    time params.parse_kofam.cluster_time

    input:
        path(infile_tblout)

    output:
        path("KOs_parsed.tsv"), emit: KOs_parsed

    script:
        """
        module load ${params.modules.tools}
        parseKofam.pl \\
            --hmmsearch \\
            --infile ${infile_tblout} \\
            --ref_database ${params.parse_kofam.ref_database} \\
            > KOs_parsed.tsv
        """
}

process KO_OVERREP {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/", mode: 'symlink'
    cpus params.kegg_overrep.cluster_cpus
    memory params.kegg_overrep.cluster_memory
    time params.kegg_overrep.cluster_time

    input:
        path(KOs_parsed)
        path(gene_abundance)

    output:
        path("KO_abundance.tsv"), emit: KO_abundance

    script:
        """
        module load ${params.modules.tools} ${params.modules.python} && \\
        getKeggKO.py \
            --infile-blastp ${KOs_parsed} \
            --infile-gene-abundance ${gene_abundance} \
            > KO_abundance.tsv
        """
}

/* 
    I've tried, but I couldn't make array jobs work with current version of Nextflow. To 
    revisit in the future...
*/
process DIAMOND_BLASTP_NR_ARRAY_JOB {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/diamond_blastp_nr/", mode: 'symlink'
    clusterOptions = params.clusterOptions + " --array=0-1"
    cpus params.diamond_blastp.cluster_cpus
    memory params.diamond_blastp.cluster_memory
    time params.diamond_blastp.cluster_time
    afterScript 'touch alldone.done'

    input:
        path(indir)
        //path(number_of_chunks_file)

    output:
        //path("diamond_blastp_nr_out_\$SLURM_ARRAY_TASK_ID2.tsv"), emit: tsv
        env(SLURM_ARRAY_TASK_ID), emit: sati
        path("alldone.done"), emit: done

    script:
        
        """
        module load ${params.modules.diamond} && \\
        SLURM_ARRAY_TASK_ID2=\$(( SLURM_ARRAY_TASK_ID - 0  ))

        SLURM_ARRAY_TASK_ID2=\$(printf "%07d" \$SLURM_ARRAY_TASK_ID2)
        echo 'SLURM/SGE_TASK_ID2' \$SLURM_ARRAY_TASK_ID2

        if [[ -s ${indir}/contigs_renamed.faa_chunk_\$SLURM_ARRAY_TASK_ID2 ]] ; then
            diamond blastp --quiet \\
                -d ${params.diamond_blastp.db_nr} \\
                -q ${indir}/contigs_renamed.faa_chunk_\$SLURM_ARRAY_TASK_ID2 \\
                -o diamond_blastp_nr_out_\$SLURM_ARRAY_TASK_ID2.tsv \\
                -k 10 \\
                -e ${params.diamond_blastp.evalue} \\
                -p ${params.diamond_blastp.num_threads}
        else
            touch diamond_blastp_nr_out.tsv
        fi ;
        """

}

process CONVERT_IDS_FOR_CAT {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/taxonomy/", mode: 'symlink'
    cpus params.CAT.cluster_cpus
    memory params.CAT.cluster_memory
    time params.CAT.cluster_time

    input:
        path(infile_gff)
        path(infile_diamond)

    output:
        path("blastp_nr_alt_orf_ids.tsv"), emit: orf_ids

    script:
        """
        module load ${params.modules.tools} && \\
        convertDiamondBlastpORFIDForCAT.pl \\
            --infile_gff ${infile_gff} \\
            --infile_blastp ${infile_diamond} \\
            > blastp_nr_alt_orf_ids.tsv
        """
}

process CAT {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/taxonomy/", mode: 'symlink'
    cpus params.CAT.cluster_cpus
    memory params.CAT.cluster_memory
    time params.CAT.cluster_time

    input:
        path(contigs_fna)
        path(genes_faa)
        path(orf_ids)

    output:
        path("out.contig2classification_with_names.tsv"), emit: classification_with_names

    script:
        """
        module load ${params.modules.CAT} && \\
        CAT contigs \\
            -r ${params.CAT.r} -f ${params.CAT.f} \\
            -c ${contigs_fna} \\
            -p ${genes_faa} \\
            -a ${orf_ids} \\
            -d ${params.CAT.database_folder} \\
            -t ${params.CAT.taxonomy_folder} \\
            -o out --force && \\
        CAT add_names \\
            -i out.contig2classification.txt \\
            -o out.contig2classification_with_names.tsv \\
            -t ${params.CAT.taxonomy_folder} --force
        """
}

process GENERATE_FEATURE_TABLES {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/taxonomy/", mode: 'symlink'
    cpus params.DEFAULT.cluster_cpus
    memory params.DEFAULT.cluster_memory
    time params.DEFAULT.cluster_time

    input:
        path(classification_with_names)
        path(contig_abundance)

    output:
        path("taxonomy.tsv"), emit: taxonomy
        path("feature_table.tsv"), emit: feature_table
        path("feature_table_bacteriaArchaea.tsv"), emit: feature_table_ba
        path("feature_table_others.tsv"), emit: feature_table_others

    script:
        """
        module load ${params.modules.tools} ${params.modules.perl} && \\
        generateFeatureTableFromCAT.pl \\
            --infile_taxonomy ${classification_with_names} \\
            --infile_abundance ${contig_abundance} \\
            --outfile_taxonomy taxonomy.tsv \\
            > feature_table.tsv && \\
        splitFeatureTable.pl \\
            --infile feature_table.tsv \\
            --matched feature_table_bacteriaArchaea.tsv \\
            --unmatched feature_table_others.tsv \\
            --keepChloro yes \\
            --keepMito yes \\
            --select bacteriaArchaea
        """
}

process COG_MATRIX_RPOB{
    debug true
    publishDir "$params.DEFAULT.outdir/alpha_diversity/rpob/", mode: 'symlink'
    cpus params.DEFAULT.cluster_cpus
    memory params.DEFAULT.cluster_memory
    time params.DEFAULT.cluster_time

    input:
        path(gene_abundance)
        path(rpsblast_cog)

    output:
        path("rpob_abundance.tsv"), emit: rpob_abundance

    script:
        """
        module load ${params.modules.tools} && \\
        generateCOGMatrix.pl \\
            --infile ${gene_abundance} \\
            --cog_term COG0085 \\
            --cog_file ${rpsblast_cog} \\
            > rpob_abundance.tsv
        """
}

process COG_MATRIX_RECA{
    debug true
    publishDir "$params.DEFAULT.outdir/alpha_diversity/reca/", mode: 'symlink'
    cpus params.DEFAULT.cluster_cpus
    memory params.DEFAULT.cluster_memory
    time params.DEFAULT.cluster_time

    input:
        path(gene_abundance)
        path(rpsblast_cog)

    output:
        path("reca_abundance.tsv"), emit: reca_abundance

    script:
        """
        module load ${params.modules.tools} && \\
        generateCOGMatrix.pl \\
            --infile ${gene_abundance} \\
            --cog_term COG0468 \\
            --cog_file ${rpsblast_cog} \\
            > reca_abundance.tsv
        """
}

process METABAT_ABUNDANCE{
    debug true
    publishDir "$params.DEFAULT.outdir/mags/", mode: 'symlink'
    cpus params.metabat_abundance.cluster_cpus
    memory params.metabat_abundance.cluster_memory
    time params.metabat_abundance.cluster_time

    input:
        path(bams)

    output:
        path("abundance.txt"), emit: abundance

    script:
        """
        module load ${params.modules.metabat2} ${params.modules.tools} && \\
        export OMP_NUM_THREADS=${params.metabat_abundance.num_threads} && \\
        jgi_summarize_bam_contig_depths \\
            --outputDepth abundance.txt \\
            --pairedContigs paired_contigs.txt \\
            --minContigLength ${params.metabat_abundance.min_contig_length} \\
            --minContigDepth ${params.metabat_abundance.min_contig_depth} \\
            --percentIdentity ${params.metabat_abundance.perc_id} \\
            ${bams}
        """
}

process METABAT2{
    debug true
    publishDir "$params.DEFAULT.outdir/mags/metabat2", mode: 'symlink'
    cpus params.metabat2.cluster_cpus
    memory params.metabat2.cluster_memory
    time params.metabat2.cluster_time

    input:
        path(contigs)
        path(abundance_profile)

    output:
        path("./fna"), emit: fasta_dir

    script:
        """
        module load ${params.modules.metabat2} ${params.modules.tools} && \\
        metabat2 \\
            -i ${contigs} \\
            -a ${abundance_profile} \\
            -o ./fna/out \\
            --maxP ${params.metabat2.max_p} \\
            -m ${params.metabat2.min_contig} \\
            -t ${params.metabat2.num_threads} \\
            --saveTNF saved.tnf \\
            --saveDistance saved.dist
        """
}

process CHECKM_METABAT2{
    debug true
    publishDir "$params.DEFAULT.outdir/mags/metabat2/", mode: 'symlink'
    cpus params.checkm.cluster_cpus
    memory params.checkm.cluster_memory
    time params.checkm.cluster_time

    input:
        path(indir)

    output:
        path("out_checkm.tsv"), emit: tsv

    script:
    """
    module load ${params.modules.hmmer} ${params.modules.prodigal} ${params.modules.pplacer} ${params.modules.checkm} && \\
    bash -c '
        set +u && source \$CHECKM_HOME/venv_checkm/bin/activate && set -u && \\
        mkdir -p out_checkm && \\
        \$CHECKM_HOME/bin/checkm lineage_wf \\
            -f out_checkm.tsv \\
            -x fa \\
            -t 4 \\
            --tab_table \\
            ${indir} out_checkm'
    """
}
