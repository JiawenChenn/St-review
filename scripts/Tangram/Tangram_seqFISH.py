import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import os
import time

# import tangram for spatial deconvolution
import tangram as tg

stdir = "/processed_data/MOB/"
stfile = 'seqfish.cnt.genexrow.tsv'
cellfile = "ob_seqfish_cell_pro.tsv"

scdir = "/processed_data/MOB/ref/internal/"
reffile = "OB.sc.cnt.tsv"
metafile = "OB.sc.mta.tsv"

resultdir = "/scripts/Tangram/MOB/"
resultfile = "MOB.Seqfish.own_all.tangram.txt"
result_cellfile = "MOB.Seqfish.own_all.tangram_cell.txt"
if not os.path.exists(resultdir):
    os.mkdir(resultdir)

# read in ST data
st = pd.read_csv(stdir + stfile, sep='\t', index_col=0)
st = st.transpose()
adata_st = AnnData(st)
adata_st

## append cell count for each spot
cell = pd.read_csv(stdir + cellfile, sep='\t', index_col=0)
adata_st.obs = adata_st.obs.merge(cell.sum(axis = 1).to_frame(name="cell_count"), how = 'outer', left_index = True, right_index = True)

## create spatial coordinates information
spatial_coord = adata_st.obs.reset_index()['index'].str.split('_', expand = True).to_numpy()
spatial_coord[:,0] = spatial_coord[:,0] + spatial_coord[:,2]
spatial_coord = spatial_coord[:, 0:2]
adata_st.obsm['spatial'] = spatial_coord

centroid = pd.Series(index = adata_st.obs.index, dtype = "object")
for i in range(len(centroid)):
    centroid[i] = np.tile(spatial_coord[i], (adata_st.obs.cell_count[i],1))

adata_st.obsm['image_features'] = cell.sum(axis = 1).to_frame(name="segmentation_label").merge(centroid.to_frame(name = "segmentation_centroid"),left_index = True, right_index = True)

# Read in scRNA-seq data
scdat = pd.read_csv(scdir + reffile, sep='\t', index_col=0)
adata_sc = AnnData(scdat)

sc_meta = pd.read_csv(scdir + metafile, sep='\t')
sc_meta.set_index('cell', inplace = True)
sc_meta.index = sc_meta.index.astype('str')
adata_sc.obs = adata_sc.obs.merge(sc_meta, how = 'left', left_index=True, right_index=True)
adata_sc.obs["bio_celltype"] = pd.Categorical(adata_sc.obs['bio_celltype'])
adata_sc

start_time = time.time()
# preprocessing: find the common genes between sc and st
tg.pp_adatas(adata_sc, adata_st, genes=None)


# Deconvolution
ad_map = tg.map_cells_to_space(
    adata_sc,
    adata_st,
    mode="constrained",
    target_count=adata_st.obs.cell_count.sum(),
    density_prior=np.array(adata_st.obs.cell_count) / adata_st.obs.cell_count.sum(),
    num_epochs=1000,
    device="cuda:0",
    #device='cpu',
)

# Gather deconvolution results
## map the cell type information to the st AnnData object
## The output created is the unnormalized probability matrix
tg.project_cell_annotations(ad_map, adata_st, annotation="bio_celltype")
end_time = time.time()
print("--- %.2f seconds ---" % (time.time() - start_time))

## normalize the probability matrix and save as csv
prob_mat = adata_st.obsm["tangram_ct_pred"]
prob_mat = prob_mat.div(prob_mat.sum(axis=1), axis=0)
prob_mat.to_csv(resultdir + resultfile, sep = '\t')


## create cell-level mapping assignments
tg.create_segment_cell_df(adata_st)
tg.count_cell_annotations(
    ad_map,
    adata_sc,
    adata_st,
    annotation="bio_celltype",
)
adata_st.obsm["tangram_ct_count"].drop(columns = ['centroids']).to_csv(resultdir + result_cellfile, sep = '\t')

