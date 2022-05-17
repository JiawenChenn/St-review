
import os.path as osp
import argparse as arp
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from os import mkdir, makedirs,getcwd

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns


prs = arp.ArgumentParser()
parser = arp.ArgumentParser()

parser.add_argument('-scc','--sc_cnt',
                    required = True,
                    type = str)

parser.add_argument('-scl','--sc_labels',
                    required = True,
                    type = str)


parser.add_argument('-stc','--st_cnt',
                    required = True,
                    type = str)

parser.add_argument('-out','--out_dir',
                    required = True,
                    type = str)

parser.add_argument('-pre','--perfix',
                    required = True,
                    type = str)

args = parser.parse_args()


#pip install git+https://github.com/BayraktarLab/cell2location.git#egg=cell2location[tutorials]


scRNA_path = args.sc_cnt
scRNA_anno_path =  args.sc_labels

output_pre = args.perfix


results_folder = args.out_dir+output_pre+"/"
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
makedirs(results_folder,exist_ok=True)

# cell * gene
scRna = pd.read_csv(scRNA_path,sep="\t",index_col=0)
#scRna=scRna.transpose() # use this step when the input is gene*cnt
anno = pd.read_csv(scRNA_anno_path,sep="\t")
anno.index=anno.iloc[:,0].values
intersect = anno.index.intersection(scRna.index)
scRna = scRna.loc[intersect,]
anno = anno.loc[intersect,]

gene_name=pd.DataFrame(scRna.columns.values,columns=['gene'],index=scRna.columns.values)


adata_ref1 = anndata.AnnData(X=scRna,obs= anno,var=gene_name)
adata_ref1.X
adata_ref1.var
adata_ref1.obs['Sample']="scRNA"

from cell2location.utils.filtering import filter_genes
selected=adata_ref1.var.index
# filter the object (not filter in our analysis)
adata_ref1 = adata_ref1[:, selected].copy()

scvi.data.setup_anndata(adata=adata_ref1,
                        batch_key='Sample',
                        labels_key=anno.columns.values[1],
                    )
scvi.data.view_anndata_setup(adata_ref1)

from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref1)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=6000, batch_size=1024, train_size=1, lr=0.002, use_gpu=True)

# plot ELBO loss history during training, removing first 20 epochs from the plot
plt.clf()
mod.plot_history(0)
plt.savefig(results_folder+"sc_loss.png")
plt.clf()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref1 = mod.export_posterior(
    adata_ref1, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref1.write(adata_file)
adata_file

adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref1 = sc.read_h5ad(adata_file)


run_name = f'{results_folder}/cell2location_map/'
st = pd.read_csv(args.st_cnt,sep="\t",index_col=0)
st=st.transpose()
anno=pd.DataFrame(st.index.values,columns=['spots'],index=st.index.values)

gene_name=pd.DataFrame(st.columns.values,columns=['gene'],index=st.columns.values)
adata_vis1 = anndata.AnnData(X=st,obs= anno,var=gene_name)
adata_vis1.obs['sample'] = "st"

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref1.varm.keys():
    inf_aver = adata_ref1.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_ref1.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref1.var[[f'means_per_cluster_mu_fg_{i}' for i in adata_ref1.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref1.uns['mod']['factor_names']
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis1.var_names, inf_aver.index)
adata_vis1 = adata_vis1[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
# prepare anndata for cell2location model
scvi.data.setup_anndata(adata=adata_vis1, batch_key="sample")
scvi.data.view_anndata_setup(adata_vis1)

sc.pp.filter_genes(adata_vis1, min_counts=1)
sc.pp.filter_cells(adata_vis1,min_genes=3)
# create and train the model
mod2 = cell2location.models.Cell2location(
    adata_vis1, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)
mod2.train(max_epochs=30000,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
        use_gpu=True)
# plot ELBO loss history during training, removing first 100 epochs from the plot
plt.clf()
mod2.plot_history(1000)
plt.legend(labels=['full data training'])
plt.savefig(f'{results_folder}mod2_train.png')
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis1 = mod2.export_posterior(
    adata_vis1, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)
# Save model
mod2.save(f"{run_name}", overwrite=True)
# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis1.write(adata_file)
adata_file
plt.clf()
mod2.plot_QC()
plt.savefig(f'{results_folder}/plot_QC.png')
plt.clf()


