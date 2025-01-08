#!/usr/bin/env python3

# import libraries =======================================================================
from spatialdata_io import visium_hd
import spatialdata as sd
import sys
import anndata as ad

# parse command line arguments ===========================================================
spaceranger_results_dir = sys.argv[1]
fullres_image_fn = sys.argv[2]
sample_name = sys.argv[3]
bin_sizes = sys.argv[4].strip("[]").replace(" ", "").split(",")

## create spatialdata object =============================================================
print("importing Visium HD data")
sdata = visium_hd(
        spaceranger_results_dir,
        dataset_id = sample_name,
        fullres_image_file = fullres_image_fn,
        load_all_images = True
        )
    
# make the var names unique ==============================================================
print("making variable names unique")
for table in sdata.tables.values():
	table.var_names_make_unique()

# rasterize expression data ==============================================================
print("rasterizing transcriptome data")
for i in bin_sizes:
	bin_size = f"{i:0>3}"
	# rasterize_bins() requires a compresed sparse column (csc) matrix
	sdata.tables[f"square_{bin_size}um"].X = sdata.tables[f"square_{bin_size}um"].X.tocsc()
	rasterized = sd.rasterize_bins(
		sdata,
		f"{sample_name}_square_{bin_size}um",
		f"square_{bin_size}um",
		"array_col",
		"array_row",
	)
	sdata[f"rasterized_{bin_size}um"] = rasterized

# add clusters to object =================================================================
import pandas as pd

for i in bin_sizes:
	bin_size = f"{i:0>3}"

	clusters_fn = f"{sample_name}_{i}um_clusters.csv.gz"
	clusters = pd.read_csv(clusters_fn)
	
	
	clusters["cluster_full"] = clusters["cluster_full"].astype("category")
	clusters.set_index("spot_id", inplace=True)
	
	join_data = sdata[f"square_{bin_size}um"].obs.merge(clusters, how="left", left_index=True, right_index=True)
	sdata[f"square_{bin_size}um"].obs = join_data
	
	sdata[f"rasterized_clusters_{bin_size}um"] = sd.rasterize_bins(
		sdata,
		f"{sample_name}_square_{bin_size}um",
		f"square_{bin_size}um",
		"array_col",
		"array_row",
		"clusters_full"
	)
	
# write sdata object to zarr file ========================================================
    # Write sdata object to file
    print("writing sdata object to zarr store")
out_path = f"{sample_name}_sdata.zarr"

out_data = sdata

images = [
f"{sample_name}_cytassist_image",
f"{sample_name}_full_image"
]

padded_bin_sizes = [f"{i:0>3}" for i in bin_sizes]
rasterized_expression = ["rasterized_{}um".format(bin_size) for bin_size in padded_bin_sizes]
rasterized_clusters = ["rasterized_clusters_{}um".format(bin_size) for bin_size in padded_bin_sizes]

keep_layers = images + rasterized_expression + rasterized_clusters

out_data = sdata.subset(keep_layers)

print(f"file_name: {out_path}")
out_data.write(out_path)