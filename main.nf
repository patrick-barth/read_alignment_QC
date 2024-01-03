#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//import processes for 
include{

} from './modules/< name module to import processes from >'

//essential input files
input_reads     = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
}
annotation = file(params.annotation)

//preparation for workflow

/*
 * Welcome log to be displayed before workflow
 */
log.info """\
         ${params.manifest.name} v${params.manifest.version}
         ==========================
         input from   : ${params.input_file}
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
workflow preprocessing {
    take: 
        input_reads
    main:
        quality_control(input_reads)
        adapter_removal(input_reads)
        quality_filter(adapter_removal.out.fastq_trimmed)
        quality_control_2(quality_filter.out.fastq_quality_filtered)

    emit:
        //data for multiqc
        multiqc_quality_control                     = quality_control.out
        multiqc_quality_control_post_preprocessing  = quality_control_2.out
        multiqc_adapter_removal                     = adapter_removal.out.report_trimming
        multiqc_quality_filter                      = quality_filter.out.report_quality_filter

        // data for downstream processes
        fastq_reads_quality_filtered                = quality_filter.out.fastq_quality_filtered
}

/*
 * 
 */

if ( params.help ) {
    help = """your_script.nf: A description of your script and maybe some examples of how
             |                to run the script
             |Required arguments:
             |  --input_file  Location of the input file file.
             |                [default: ${params.input_file}]
             |
             |Optional arguments:
             |  --use_thing   Do some optional process.
             |                [default: ${params.use_thing}]
             |  -w            The NextFlow work directory. Delete the directory once the process
             |                is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

/*
 * Actual workflow connecting subworkflows
 */
workflow {
    preprocessing(input_reads)
}


/*
 * Prints complection status to command line
 */
workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}