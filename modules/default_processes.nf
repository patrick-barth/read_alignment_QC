process collect_metadata {
	publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "pipeline_metadata.txt"

    output:
    path("pipeline_metadata.txt")

    script:
    """
    cat <<EOF > pipeline_metadata.txt
    Author: ${params.manifest.author}
    Pipeline version: ${params.manifest.version}
    Working directory: ${workflow.workDir}
    User name: ${workflow.userName}
    Command line: ${workflow.commandLine}
    Container engine: ${workflow.containerEngine}    
    Containers used: ${workflow.container}
    Git repository: ${workflow.repository}
    Repository revision: ${workflow.revision}
    """
}