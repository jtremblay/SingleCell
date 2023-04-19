/*
    Microbial ecology functions 
*/

process BETA_DIVERSITY_BACTARCH{
    debug true
    publishDir "$params.DEFAULT.outdir/beta_diversity/bacteria_archaea/", mode: 'symlink'
    cpus params.beta_diversity.cluster_cpus
    memory params.beta_diversity.cluster_memory
    time params.beta_diversity.cluster_time

    input:
        path(feature_table)

    output:
        path("bray_curtis_contig_bacteriaArchaea.tsv"), emit: distance_matrix
        path("bray_curtis_contig_bacteriaArchaea_coords.tsv"), emit: coords
        path("3d_bray_curtis_plot"), emit: interactive_plot

    script:
        """
        module load ${params.modules.tools} ${params.modules.python} && \\
        microbiomeutils.py betadiv \\
            -i ${feature_table} \\
            -m bray-curtis \\
            > bray_curtis_contig_bacteriaArchaea.tsv && \\
            sed -i 's/nan/0.0/g' bray_curtis_contig_bacteriaArchaea.tsv && \\
        microbiomeutils.py pcoa \\
            -i bray_curtis_contig_bacteriaArchaea.tsv \\
            > bray_curtis_contig_bacteriaArchaea_coords.tsv && \\
        microbiomeutils.py emperor \\
            -i bray_curtis_contig_bacteriaArchaea_coords.tsv \\
            -m ${params.DEFAULT.mapping_file} \\
            -o 3d_bray_curtis_plot
        """
}

process BETA_DIVERSITY_ALL{
    debug true
    publishDir "$params.DEFAULT.outdir/beta_diversity/all/", mode: 'symlink'
    cpus params.beta_diversity.cluster_cpus
    memory params.beta_diversity.cluster_memory
    time params.beta_diversity.cluster_time

    input:
        path(feature_table)

    output:
        path("bray_curtis_contig.tsv"), emit: distance_matrix
        path("bray_curtis_contig_coords.tsv"), emit: coords
        path("3d_bray_curtis_plot"), emit: interactive_plot

    script:
        """
        module load ${params.modules.tools} ${params.modules.python} && \\
        microbiomeutils.py betadiv \\
            -i ${feature_table} \\
            -m bray-curtis \\
            > bray_curtis_contig.tsv && \\
            sed -i 's/nan/0.0/g' bray_curtis_contig.tsv && \\
        microbiomeutils.py pcoa \\
            -i bray_curtis_contig.tsv \\
            > bray_curtis_contig_coords.tsv && \\
        microbiomeutils.py emperor \\
            -i bray_curtis_contig_coords.tsv \\
            -m ${params.DEFAULT.mapping_file} \\
            -o 3d_bray_curtis_plot
        """
}

process ALPHA_DIVERSITY_CONTIGS{
    debug true
    publishDir "$params.DEFAULT.outdir/alpha_diversity/contig_abundance/", mode: 'symlink'
    cpus params.alpha_diversity.cluster_cpus
    memory params.alpha_diversity.cluster_memory
    time params.alpha_diversity.cluster_time

    input:
        path(contig_abundance)

    output:
        path("*.tsv"), emit: diversity_index_files

    script:
        """
        module load ${params.modules.tools} ${params.modules.rtk} && \\
        generateDepthPointsForRarefaction.pl \\
            --remove_last_col false \\
            --infile ${contig_abundance} \\
            --number_of_points ${params.alpha_diversity.number_of_points} \\
            > depth_list.txt && \\
        sed '1 s/^\\S\\+\\t/\\t/' ${contig_abundance} > tmp.tsv && \\
        rtk swap \\
            -i tmp.tsv \\
            -o ./ \\
            -r ${params.alpha_diversity.perm} \\
            -d `cat depth_list.txt` && \\
        rm tmp.tsv && \\
        touch rtk.done
        """
}

process ALPHA_DIVERSITY_GENES{
    debug true
    publishDir "$params.DEFAULT.outdir/alpha_diversity/gene_abundance/", mode: 'symlink'
    cpus params.alpha_diversity.cluster_cpus
    memory params.alpha_diversity.cluster_memory
    time params.alpha_diversity.cluster_time

    input:
        path(gene_abundance)

    output:
        path("*.tsv"), emit: diversity_index_files

    script:
        """
        module load ${params.modules.tools} ${params.modules.rtk} && \\
        generateDepthPointsForRarefaction.pl \\
            --remove_last_col false \\
            --infile ${gene_abundance} \\
            --number_of_points ${params.alpha_diversity.number_of_points} \\
            > depth_list.txt && \\
        sed '1 s/^\\S\\+\\t/\\t/' ${gene_abundance} > tmp.tsv && \\
        rtk swap \\
            -i tmp.tsv \\
            -o ./ \\
            -r ${params.alpha_diversity.perm} \\
            -d `cat depth_list.txt` && \\
        rm tmp.tsv && \\
        touch rtk.done
        """
}

process SUMMARIZE_TAXONOMY {
    debug true
    publishDir "$params.DEFAULT.outdir/annotations/taxonomy/", mode: 'symlink'
    cpus params.summarize_taxonomy.cluster_cpus
    memory params.summarize_taxonomy.cluster_memory
    time params.summarize_taxonomy.cluster_time

    input:
        path(feature_table)

    output:
        path("feature_table_L{1,2,3,4,5,6,7}_relative.tsv"), emit: relative
        path("feature_table_L{1,2,3,4,5,6,7}_absolute.tsv"), emit: absolute

    script:
        """
        module load ${params.modules.tools} ${params.modules.python} && \\
        for i in {1..7}
        do
            microbiomeutils.py taxsum \\
                -i ${feature_table} \\
                -l \$i \\
                -t relative \\
                > feature_table_L\${i}_relative.tsv && \\
            microbiomeutils.py taxsum \\
                -i ${feature_table} \\
                -l \$i \\
                -t absolute \\
                > feature_table_L\${i}_absolute.tsv
        done
        """
}

process ALPHA_DIVERSITY_RPOB{
    debug true
    publishDir "$params.DEFAULT.outdir/alpha_diversity/rpob/", mode: 'symlink'
    cpus params.alpha_diversity.cluster_cpus
    memory params.alpha_diversity.cluster_memory
    time params.alpha_diversity.cluster_time

    input:
        path(gene_abundance)

    output:
        path("*.tsv"), emit: diversity_index_files

    script:
        """
        module load ${params.modules.tools} ${params.modules.rtk} && \\
        generateDepthPointsForRarefaction.pl \\
            --remove_last_col false \\
            --infile ${gene_abundance} \\
            --number_of_points ${params.alpha_diversity.number_of_points} \\
            > depth_list.txt && \\
        sed '1 s/^\\S\\+\\t/\\t/' ${gene_abundance} > tmp.tsv && \\
        rtk swap \\
            -i tmp.tsv \\
            -o ./ \\
            -r ${params.alpha_diversity.perm} \\
            -d `cat depth_list.txt` && \\
        rm tmp.tsv && \\
        touch rtk.done
        """
}

process ALPHA_DIVERSITY_RECA{
    debug true
    publishDir "$params.DEFAULT.outdir/alpha_diversity/reca/", mode: 'symlink'
    cpus params.alpha_diversity.cluster_cpus
    memory params.alpha_diversity.cluster_memory
    time params.alpha_diversity.cluster_time

    input:
        path(gene_abundance)

    output:
        path("*.tsv"), emit: diversity_index_files

    script:
        """
        module load ${params.modules.tools} ${params.modules.rtk} && \\
        generateDepthPointsForRarefaction.pl \\
            --remove_last_col false \\
            --infile ${gene_abundance} \\
            --number_of_points ${params.alpha_diversity.number_of_points} \\
            > depth_list.txt && \\
        sed '1 s/^\\S\\+\\t/\\t/' ${gene_abundance} > tmp.tsv && \\
        rtk swap \\
            -i tmp.tsv \\
            -o ./ \\
            -r ${params.alpha_diversity.perm} \\
            -d `cat depth_list.txt` && \\
        rm tmp.tsv && \\
        touch rtk.done
        """
}
