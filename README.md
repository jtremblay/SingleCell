# SingleCell
Single cell bioinformatics pipeline

```
                       _____ _             _       _____     _ _ 
                      / ____(_)           | |     / ____|   | | |
                     | (___  _ _ __   __ _| | ___| |     ___| | |
                      \___ \| | '_ \ / _\'| |/ _ \ |    / _ \ | |
                      ____) | | | | | (_| | |  __/ |___|  __/ | |
                     |_____/|_|_| |_|\__, |_|\___|\_____\___|_|_|
                                        __/ |                      
                                       |___/   for  N E X T F L O W 
                                              
                     Support: jtremblay514@gmail.com
                   Home page: github.com/jtremblay/SingleCell
```

## Background
This is a simple pipeline to process single cell sequencing libraries of the 10X genomics type. 
Briefly, each library has to be demultiplexed so that we have one set of R1 and R2 fastqs per library.
Then reads 3' reads (i.e. .R2) are mapped (STAR aligner) against the corresponding reference genome. 
Reference genome has to be formatted before mapping the reads. For instance for the mouse genome, we'd have the following command:
```
module load STAR/2.7.10b
STAR \
    --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir ./ \
    --genomeFastaFiles ../refdata-gex-mm10-2020-A/fasta/genome.fa \
    --sjdbGTFfile ../refdata-gex-mm10-2020-A/genes/genes.gtf
```

Then make sure that the singlecell.config file is all setup properly. For instance:
```
         raw_reads = "$projectDir/raw_reads/*_R{1,2}_001.fastq"
         outdir = "$projectDir/output/"
         ref_genome = "$INSTALL_HOME/databases/single_cell/refdata-gex-mm10-2020-A_alt"
         whitelist = "$INSTALL_HOME/databases/single_cell/3M-february-2018.txt"
         umi_length = 12
         # The next parameters are essentially for the generation of the .h5 (hdf5) file generated
         # at the end of the pipeline.
         gene_info = "$INSTALL_HOME/databases/single_cell/refdata-gex-mm10-2020-A_alt/geneInfo.tab"
         chemistry = "Chromium V3"
         genome_info = "gencode_m10_GRCm38p4"
```

## Execution
Nextflow can then be run:
```
module load nextflow/22.10.7.5853 java/17.0.2
nextflow run -c ./singlecell.config ./singlecell.nf -resume
```

## Results
The results are located in the `./ouput/star/` folder. The important files are the .h5 files which contains the abundance of features (rows) x cells (columns).
These files can be read for example using the `Seurat::Read10X_h5()` function (Seurat library).



