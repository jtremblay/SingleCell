/*
    parameters are defined in the singlecell.config file
*/


log.info """
###############################################################################
                 _____ _             _       _____     _ _ 
                / ____(_)           | |     / ____|   | | |
               | (___  _ _ __   __ _| | ___| |     ___| | |
                \\___ \\| | '_ \\ / _\' | |/ _ \\ |    / _ \\ | |
                ____) | | | | | (_| | |  __/ |___|  __/ | |
               |_____/|_|_| |_|\\__, |_|\\___|\\_____\\___|_|_|
                                  __/ |                      
                                 |___/   for  N E X T F L O W                       
                                        

               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/pipelines.html
               Version: 1.4.0-beta
###############################################################################

 reads   : ${params.DEFAULT.raw_reads}
 outdir  : ${params.DEFAULT.outdir}
    """.stripIndent()
/*
  Import processes from external files
  It is common to name processes with UPPERCASE strings, to make
  the program more readable
*/
include { STAR } from './modules/singlecell'
//include { SUMMARIZE_TAXONOMY; BETA_DIVERSITY_BACTARCH; BETA_DIVERSITY_ALL; ALPHA_DIVERSITY_CONTIGS; ALPHA_DIVERSITY_GENES; ALPHA_DIVERSITY_RPOB; ALPHA_DIVERSITY_RECA } from './modules/microbial_ecology'
raw_reads_channel = Channel.fromFilePairs(params.DEFAULT.raw_reads, checkIfExists:true)

workflow {
    /*
        STAR
    */
    STAR(raw_reads_channel)

    // If host/contaminant DNA is to be removed
}

