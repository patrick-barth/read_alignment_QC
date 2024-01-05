/*
 * Checks reads for general metrics
 * Input: [FASTQ] Unpreprocessed reads 
 * Output: [HTML] General report for unpreprocessed reads 
 */
process quality_control {
	tag {query.simpleName}
	publishDir "${params.output_dir}/statistics/qc-preprocessing", mode: 'copy', pattern: "${query.baseName}_fastqc.{html,zip}"

	
	input:
	path query

	output:
	path "${query.baseName}_fastqc.{html,zip}", 				emit: output
	path("${task.process}.version.txt"), 	emit: version

	"""
	fastqc ${query} -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

/*
 * Checks preprocessed reads for general metrics
 * Input: [FASTQ] Preprocessed reads 
 * Output: [HTML] General report for preprocessed reads
 */
process quality_control_2 {
	tag {query.simpleName}
	publishDir "${params.output_dir}/statistics/qc-postprocessing", mode: 'copy', pattern: "${query.simpleName}_2_fastqc.{html,zip}"

	
	input:
	path query

	output:
	path "${query.simpleName}_2_fastqc.{html,zip}", emit: output
	path("${task.process}.version.txt"), 			emit: version

	"""
	cat ${query} > ${query.simpleName}_2.fastq
	fastqc ${query.simpleName}_2.fastq -o .

	echo -e "${task.process}\tFastQC\t\$(fastqc --version | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

/*
 * Removes adapters from reads
 * Input: [FASTQ] Read file
 * Params: params.min_length -> Minimum length reads need after trimming to not be omitted
 * Output:  fastq_trimmed 	-> [FASTQ] Read file with afdapters trimmed
 *			report_trimming -> [TXT] Report adapter trimming
 */
process adapter_removal {
	tag {query.simpleName}
	publishDir "${params.output_dir}/statistics", mode: 'copy', pattern: "${query}_trimming_report.txt"


	input:
	path query

	output:
	path "${query}_trimmed.fq", 														emit: fastq_trimmed 
	path "${query}_trimming_report.txt", 												emit: report
	tuple path("${task.process}.version.txt"), path("${task.process}.version2.txt"), 	emit: version

	"""
	trim_galore --cores ${task.cpus} --basename ${query} -o . --length ${params.min_length} ${query} --quality 0

	echo -e "${task.process}\ttrim_galore\t\$(trim_galore -v | head -4 | tail -1 | sed -e 's/^[ \t]*//' | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tcutadapt\t\$(cutadapt --version)" > ${task.process}.version2.txt
	"""
}

/*
 * Removes bases with low quality from reads
 * Input: [FASTQ] Read file
 * Params: 	params.min_qual					-> Bases below this threshold are omitted 
 *			params.min_percent_qual_filter	-> Minimum percentage of bases of a read need to be above this threshold to keep the it 
 * Output: 	fastq_quality_filtered 	-> [FASTQ] Read file with low quality bases filtered out
 *			report_quality_filter 	-> [TXT] Report of quality filtering
 */
process quality_filter {
	tag {query.simpleName}
	publishDir "${params.output_dir}/statistics", mode: 'copy', pattern: "summary-quality-filter.txt"
	publishDir "${params.output_dir}/preprocessed-reads", mode: 'copy', pattern: "${query.baseName}.qual-filter.fastq"

	input:
	path query

	output:
	path "${query.baseName}.qual-filter.fastq", emit: fastq_quality_filtered 
	path 'summary-quality-filter.txt', 			emit: report
	path("${task.process}.version.txt"), 		emit: version


	"""
	fastq_quality_filter -v -q ${params.min_qual} -p ${params.min_percent_qual_filter} -i ${query} -o ${query.baseName}.qual-filter.fastq > summary-quality-filter.txt

	echo -e "${task.process}\tfastq_quality_filter\t\$(fastq_quality_filter -h | head -2 | tail -1 | rev | cut -d ' ' -f 5- | rev)" > ${task.process}.version.txt
	"""
}