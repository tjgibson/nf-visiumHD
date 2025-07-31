## Overview
This is a minimal nextflow pipeline for processing data generated with 10X genomics probe-based visium HD assay. 
The pipeline performs the following analysis steps:
1. Process raw Visium HD data using Spaceranger (currently v 4.0.1)
2. Use Seurat to perform clustering
3. Generate a spatialdata object that can be used for downstream processing with the spatialdata python package

## Inputs
This pipeline requires the following input files:
1. A samplesheet describing the input samples (more details below)
2. A probeset csv file. Standard probesets are avaiable from 10X Genomics, as are instructions for making custom probesets.
3. A Spaceranger reference genome. Standard references are also available from 10X genomics, or you can build a custom reference using Spaceranger mkref

## Usage
The pipeline can be run with the following command:
```
nextflow run \
    https://github.com/tjgibson/nf-visiumHD \
    --samplesheet samplesheet.csv \
    --results_dir <RESULTSDIR> \
    -profile <docker/singularity/...> \
    --probeset probeset.csv \
     --reference spaceranger_reference_path

```
## The samplesheet
The samplesheet should be a CSV file with the following columns:
1. **sample**: sample name
2. **slide_serial**: the unique slide serial number for each VisiumHD slide. This can be found in the file name for the cytassist images
3. **fastq**: path to the fastq files for a given sample. e.g. "data/fastq/sample_1/*.fastq.gz"
4. **cytaimage**: cytassist image taken on the Cytassist instrument during sample processing.
5. **image**: file path for an H&E or fluorescent image for each sample
6. **image_type**: should be one of the following: 'brightfield' for H&E images, 'darkimage' for grayscale tiff, 'colorizedImage' for colorized tiff. See [here](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/image-image-recommendation) for more details
7. **DAPI_index**: which channel of the fluorescent image file has DAPI, 1-indexed.
8. **alignment_file**: path to json file from Loupe Browser if using custom alignment or tissue detection. If not performing manual image alignment, enter 'assets/NO_FILE' to proceed with automatic image processing. Create an empty file 'assets/NO_FILE' from the directory where you will launch the nextflow pipeline.

See [samplesheet.csv](samplesheet.csv) for an example.

## Pipeline output
The pipeline produces the following output subdirectories:
- clusters: PDFs of UMAP embeddings and CSV files of cluster assignments for each bin size
- reports: nextflow reports
- spaceranger: outs/ folder from Spaceranger for each sample
- spatialdata_objects/ zarr store containing spatialdata object for each sample

## To-do:
- Add integration and integrated clusters
- Add integrated clustering of segmentation-based bins for H&E images