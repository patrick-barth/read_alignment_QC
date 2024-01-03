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

process get_md5sum {
    publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "md5sums.txt"

    input:
    path(query)

    output:
    path("md5sums.txt")

    script:
    """
    for i in ${query}
	do
		if test -f \$i; then
			md5sum \$i >> md5sums.txt
		fi
	done
    """
}