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
	publishDir "${params.results_dir}/spaceranger/${meta.sample}", mode: 'copy'
	container = "wailz/spaceranger:4.0.1"
	
	input:
	tuple val(meta), path(fastq), path(cytaimage), path(image), path(alignment_file)
    path(reference)
    path(probeset)	
    
	output:
	tuple val(meta), path("outs"), path(image)
	
	script:
	def alignment_param = alignment_file.baseName != 'NO_FILE' ? "--loupe-alignment=${alignment_file}" : ''
	
	if(meta.image_type == "brightfield"){
		image_params = "--darkimage ${image} --dapi-index ${meta.DAPI_index}"
	} else if(meta.image_type == "darkimage") { 
		image_params = "--darkimage ${image} --dapi-index ${meta.DAPI_index}"
	} else if (meta.image_type == "colorizedImage") {
		image_params = "--colorizedImage ${image} --dapi-index ${meta.DAPI_index}"
	}

    """
	spaceranger count \
		--id=${meta.sample} \
        --transcriptome=${reference} \
        --probe-set=${probeset} \
        --fastqs=. \
        --cytaimage=${cytaimage} \
        $image_params \
        --slide=${meta.slide_serial} \
        --area=${meta.slide_area} \
        --create-bam=false \
        $alignment_param \
        --custom-bin-size=4 \
        --localcores=${task.cpus} \
        --localmem=${task.memory.toGiga()}
        
	mv ${meta.sample}/outs outs
	mv outs/feature_slice.h5 outs/${meta.sample}_feature_slice.h5
	mv outs/web_summary.html outs/${meta.sample}_web_summary.html
	mv outs/cloupe_008um.cloupe outs/${meta.sample}_008um.cloupe
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
 	container = "tjmgibson/seurat_spatial:latest"
 	
 	input:
 	tuple val(meta), path("outs"), path(image), val(bin_size)
 	val(n_sketch_cells)
 	val(cluster_res)
 	val(cluster_npcs)
 	
 	output:
	tuple val(meta), path("${meta.sample}_${bin_size}um_clusters.csv.gz"), emit: clusters
 	path("${meta.sample}_${bin_size}um_clusters_UMAP.pdf"), emit: UMAP
 	
 	script:
 	"""
 	cluster_square_bins.R $bin_size $n_sketch_cells $cluster_res $cluster_npcs $meta.sample
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
 process create_sdata {
	tag "$meta.sample"
	label "process_low" 
	publishDir "${params.results_dir}/spatialdata_objects/", mode: 'copy'
 	container = "ghcr.io/scverse/spatialdata:spatialdata0.3.0_spatialdata-io0.1.7_spatialdata-plot0.2.9"
 	
 	input:
 	tuple val(meta), path("outs"), path(image), path("*")
 	val(bin_sizes)
 	
 	output:
	tuple val(meta), path("${meta.sample}_sdata.zarr.gz")
 	
 	script:
 	def bin_sizes_str = bin_sizes.join(',')
 	"""
 	create_spatialdata_object.py outs/ ${image} ${meta.sample} $bin_sizes_str
	gzip '${meta.sample}_sdata.zarr'
 	"""
 	
 	stub:
 	"""
	touch '${meta.sample}_sdata.zarr'
	gzip '${meta.sample}_sdata.zarr'
 	"""
}
 
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
        fastq_meta = row.subMap('sample', 'slide_serial', 'slide_area', 'image_type', "DAPI_index")
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
	.join(cluster_ch)
	
	
	sdata_ch = create_sdata(
	spatialdata_input,
	params.bin_sizes
	)

	
	
}