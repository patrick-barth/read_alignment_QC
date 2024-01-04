process collect_metadata {
	publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "pipeline_metadata.txt"

    output:
    path("pipeline_metadata.txt"), emit: metadata_output
    //tuple val('collect_metadata'), val('cat'), cmd("cat --version | head -1 | rev | cut -f 1 -d' ' | rev"), emit: version
    tuple val('collect_metadata'), val('cat'), env(VERSION), emit: version

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
    EOF

    VERSION=\$(cat --version | head -1 | rev | cut -f 1 -d' ' | rev)
    """
}

process get_md5sum {
    publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "md5sums.txt"

    input:
    path(query)

    output:
    path("md5sums.txt"), emit: metadata_output
    tuple val('collect_metadata'), val('cat'), env(VERSION), emit: version

    script:
    """
    for i in ${query}
	do
		if test -f \$i; then
			md5sum \$i >> md5sums.txt
		fi
	done

    VERSION=\$(md5sum --version | head -1 | rev | cut -f 1 -d' ' | rev)
    """
}

process collect_versions {
    publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "tool_versions.txt"

    input:
    tuple val(process), val(tool), val(version)

    output:
    path('tool_versions.txt')

    script:
    """
    echo -e "Process\ttool\tversion" > tool_versions.txt
    for i in ${query}
	do
		echo -e "${process}\t${tool}\t${version}" >> tool_versions.txt
	done
    """
}