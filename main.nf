#! /usr/bin/env nextflow

log.info """
	Visium HD pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()

/* input files:	
 * samplesheet
 * fastq files for each sample
 * spaceranger ref
 * probeset
 * slide serial number
 * cytassist image
 * optional HE or IF image
 */
 
 
  /*
 * Run spaceranger
 */

process spaceranger_count {
 	tag "$meta.sample"
 	label "process_high"
	publishDir "${params.results_dir}/spaceranger/${meta.sample}", mode: 'copy'
	container = "cumulusprod/spaceranger:3.1.2"
	
	input:
	tuple val(meta), path(fastq), path(image), path(cytaimage), path(alignment_file)
    path(reference)
    path(probeset)	
    
	output:
	tuple val(meta), path("outs")
	
	script:
	def alignment_param = alignment_file.baseName != 'NO_FILE' ? "--loupe-alignment=${alignment_file}" : ''
    """
	spaceranger count \
		--id=${meta.sample} \
        --transcriptome=${reference} \
        --probe-set=${probeset} \
        --fastqs=. \
        --cytaimage=${cytaimage} \
        --image ${image} \
        --slide=${meta.slide_serial} \
        --area=${meta.slide_area} \
        --create-bam=false \
        $alignment_param \
        --nosecondary \
        --custom-bin-size=4 \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()}
        
	mv ${meta.sample}/outs outs
	"""
    
    stub:
    def alignment_param = alignment_file.baseName != 'NO_FILE' ? "--loupe-alignment=${alignment_file}" : ""
    """
    mkdir outs
    touch outs/test_file.txt
    """	
} 
 
 
  /*
 * Perform cell segmentation
 */
 
  /*
 * Perform square bin clustering at all bin sizes
 */
 
process cluster_bins {
	tag "$meta.sample"
	label "process_medium"
	publishDir "${params.results_dir}/clusters/${bin_size}_um_square_bins", mode: 'copy'
 	container = "tjmgison/seurat_v5:latest"
 	
 	input:
 	tuple val(meta), path("outs"), val(bin_size)
 	val(n_sketch_cells)
 	val(cluster_res)
 	val(cluster_npcs)
 	
 	output:
	tuple val(meta), path("${meta.sample}_${bin_size}um_clusters.csv.gz"), emit: clusters
 	path("${meta.sample}_${bin_size}um_clusters_UMAP.pdf"), emit: UMAP
 	
 	script:
 	"""
 	./cluster_square_bins.R $bin_size $n_sketch_cells $cluster_res $cluster_npcs $meta.sample
 	"""
 	
 	stub:
 	"""
 	touch "${meta.sample}_${bin_size}um_clusters.csv.gz"
 	touch "${meta.sample}_${bin_size}um_clusters_UMAP.pdf"
 	"""
}
 
  /*
 * Perform cell-based clustering
 */
 
  /*
 * create base spatialdata object
 */
//  process create_sdata {
// 	tag "$meta.sample"
// 	publishDir "${params.results_dir}/spatialdata_objects/", mode: 'copy'
//  	container = "tjmgison/seurat_v5:latest"
//  	
//  	input:
//  	tuple val(meta), path("outs"), path("*.clusters.csv.gz")
//  	
//  	output:
// 	tuple val(meta), path("${meta.sample}_${bin_size}um_clusters.csv.gz"), emit: clusters
//  	path("${meta.sample}_${bin_size}um_clusters_UMAP.pdf"), emit: UMAP
//  	
//  	script:
//  	"""
//  	./cluster_square_bins.R $bin_size $n_sketch_cells $cluster_res $cluster_npcs $meta.sample
//  	"""
//  	
//  	stub:
//  	"""
//  	touch "${meta.sample}_${bin_size}um_clusters.csv.gz"
//  	touch "${meta.sample}_${bin_size}um_clusters_UMAP.pdf"
//  	"""
// }
//  
  /*
 * Add additional information to spatialdata object:
 * segmentation masks and cell-based counts
 * cell-based clusters
 * square bin clusters
 * metadata for all bin sizes and cell bins: total counts, n genes, etc
 */
 
  /*
 * generate reports for each set of clusters:
 * UMAP embeddings
 * clusters overlaid on spatial image
 * lists of marker genes for clusters
 */
 
 
 /*
 * Run workflow
 */
 
workflow {
	
	fastq_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
	| splitCsv( header:true )
    | map { row ->
        fastq_meta = row.subMap('sample', 'slide_serial', 'slide_area')
        [
        	fastq_meta, 
        	file(row.fastq, checkIfExists: true),
            file(row.cytaimage, checkIfExists: true),
            file(row.image, checkIfExists: true),
            file(row.alignment_file, checkIfExists: true)
            ]
    }
    
    spaceranger_ch = spaceranger_count(
    fastq_ch,
    file("${params.reference}", checkIfExists: true),
	file("${params.probeset}", checkIfExists: true)
	)
	
	
	bin_size_ch = Channel.fromList(params.bin_sizes)

	
	bin_ch = spaceranger_ch
		.combine(bin_size_ch) 
	
	cluster_ch = cluster_bins(
	bin_ch,
 	params.n_sketch_cells,
 	params.cluster_res,
 	params.cluster_npcs
	)
	.clusters
	| groupTuple
	
	spatialdata_input = spaceranger_ch
	.mix(cluster_ch)
	| groupTuple

}
 