#!/usr/bin/env nextflow

import groovy.json.JsonOutput // used for parameter output 

nextflow.enable.dsl=2

//import processes for 
include{
    quality_control
    quality_control_2
    adapter_removal
    quality_filter
} from './modules/read_processing.nf'

include{
    collect_metadata
    get_md5sum
} from './modules/default_processes.nf'

/*
 * Prints help and exits workflow afterwards when parameter --help is set to true
 */

if ( params.help ) {
    help = """main.nf: Collects processing data from preprocessing steps (adapter trimming and quality filtering) as well as from alignment steps and returns them. Provides an overview of the general quality of the data.
                |Required arguments:
            --reads         Location of the input file file (FASTQ).
                |
                |Optional arguments:
                |  --min_length    Minimum length for reads after adapter trimming.
                |                  [default: ${params.min_length}]
                |  --min_qual      Minimum base quality.
                |                  [default: ${params.min_qual}]
                |  --min_percent_qual_filter   Minimum percentage of bases within a read that need to be above the quality threshold
                |                              [default: ${params.min_percent_qual_filter}]
                |  -w              The NextFlow work directory. Delete the directory once the process
                |                  is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

//essential input files
input_reads     = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
}
annotation = file(params.annotation)

/*
 * preparation for workflow
 */

/*
 * Welcome log to be displayed before workflow
 */
log.info """\
        ${params.manifest.name} v${params.manifest.version}
        ==========================
        input reads  : ${params.reads}
        output to    : ${params.output_dir}
        --
        run as       : ${workflow.commandLine}
        started at   : ${workflow.start}
        config files : ${workflow.configFiles}
        """
        .stripIndent()


/*
 * Starting subworkflow descriptions
 */

/*
 * Actual workflow connecting subworkflows
 */
workflow {
    quality_control(input_reads)
    adapter_removal(input_reads)
    quality_filter(adapter_removal.out.fastq_trimmed)
    quality_control_2(quality_filter.out.fastq_quality_filtered)

    // Collect metadata
    collect_metadata()
    get_md5sum(input_reads.concat(annotation))
}


/*
 * Prints complection status to command line
 */
workflow.onComplete{
    // Create JSON file with all parameters -> will be saved in metadata output
    jsonStr = JsonOutput.toJson(params)
	file("${params.output_dir}/metadata/params.json").text = JsonOutput.prettyPrint(jsonStr)

	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}