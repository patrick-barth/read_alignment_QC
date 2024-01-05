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
    build_index_bowtie
    mapping_bowtie
} from './modules/alignment.nf'

include{
    collect_metadata
    get_md5sum
    collect_versions
} from './modules/default_processes.nf'

/*
 * Prints help and exits workflow afterwards when parameter --help is set to true
 */

if ( params.help ) {
    help = """main.nf: Collects processing data from preprocessing steps (adapter trimming and quality filtering) as well as from alignment steps and returns them. Provides an overview of the general quality of the data.
                |Required arguments:
                |   --reads         Location of the input file file (FASTQ).
                |
                |Optional arguments:
                |   --aligner       States which alignment tool is used. Currently available are: 
                |                   'bowtie2' and 'star'
                |                   [default: ${params.aligner}]
                |   --min_length    Minimum length for reads after adapter trimming.
                |                   [default: ${params.min_length}]
                |   --min_qual      Minimum base quality.
                |                   [default: ${params.min_qual}]
                |   --min_percent_qual_filter   Minimum percentage of bases within a read that need to
                |                               be above the quality threshold
                |                               [default: ${params.min_percent_qual_filter}]
                |   -w              The NextFlow work directory. Delete the directory once the process
                |                   is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

//essential input files
input_reads     = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
reference       = Channel.fromPath( params.reference )
//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
}
annotation = file(params.annotation)

/*
 * preparation for workflow
 */

// Collect all input files
input_files = input_reads.concat(Channel.of(annotation))
                    .concat(reference)
                    .flatten().toList()


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

workflow alignment {
    take:
        reference
        annotation
        reads

    main:
        if(params.aligner == "bowtie2"){
            build_index_bowtie(reference)
            mapping_bowtie(build_index_bowtie.out.index.first(),
                            reads)

            alignments_tmp          =   mapping_bowtie.out.bam_alignments
            version_index_tmp       =   build_index_bowtie.out.version
            version_align_tmp   =   mapping_bowtie.out.version
            report_tmp              =   mapping_bowtie.out.report
        } else if (params.aligner == "star"){
            //TODO: STAR logic
        } 

    emit:
        version_index   =   version_index_tmp
        version_align   =   version_align_tmp
        reports         =   report_tmp

        alignments      =   alignments_tmp
}

/*
 * Actual workflow connecting subworkflows
 */
workflow {
    // Preprocessing
    quality_control(input_reads)
    adapter_removal(input_reads)
    quality_filter(adapter_removal.out.fastq_trimmed)
    quality_control_2(quality_filter.out.fastq_quality_filtered)

    alignment(reference,
            annotation,
            quality_filter.out.fastq_quality_filtered
    )

    // Collect metadata
    collect_metadata()
    get_md5sum(input_files)
    collect_versions(collect_metadata.out.version
                        .concat(get_md5sum.out.version)
                        .concat(quality_control.out.version)
                        .concat(quality_control_2.out.version)
                        .concat(adapter_removal.out.version)
                        .concat(quality_filter.out.version)
                        .concat(alignment.out.version_index)
                        .concat(alignment.out.version_align)
                        .unique()
                        .flatten().toList()
                    )
}

/*
 * Check for parameter correctness
 */
def validAligners = ["bowtie2","star"]
if(params.aligner !in validAligners){
    log.info """\
            Stated aligner (${params.aligner}) is not a valid option.
            Available aligners are: "bowtie2" and "star".
            Please use one of them via the parameter --aligner and restart the pipeline.
            """
            .stripIndent()
            exit(1)
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