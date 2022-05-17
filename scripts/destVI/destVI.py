import sys

#if True, will install via pypi, else will install from source
stable = False
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sns
from matplotlib.lines import Line2D
import umap

import torch
import anndata
import time
import argparse as arp
import scvi
from scvi.model import CondSCVI, DestVI

from os import mkdir, makedirs,getcwd

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

import sys

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


outdir=args.out_dir+args.perfix+"/"

# process scRNA data
scRna = pd.read_csv(args.sc_cnt,sep="\t",index_col=0)
#scRna=scRna.transpose() # use this step when the input is gene*cnt

anno = pd.read_csv(args.sc_labels,sep="\t",index_col=0)
intersect = anno.index.intersection(scRna.index)
scRna = scRna.loc[intersect,]
anno = anno.loc[intersect,]
gene_name=pd.DataFrame(scRna.columns.values,columns=['gene'],index=scRna.columns.values)

adata_ref1 = anndata.AnnData(X=scRna,obs= anno,var=gene_name)
adata_ref1.obs['Sample']="scRNA"
sc.pp.filter_genes(adata_ref1, min_counts=1)
sc.pp.filter_cells(adata_ref1,min_genes=3)
adata_ref1.layers["counts"] = adata_ref1.X.copy()
sc.pp.normalize_total(adata_ref1, target_sum=10e4)
sc.pp.log1p(adata_ref1)
sc_adata = adata_ref1.copy()

# fit scRNA model
CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key=anno.columns.values[0])
sc_model = CondSCVI(sc_adata, weight_obs=True)
sc_model.train(max_epochs=500)

# save training figures
sc_model.history["elbo_train"].plot()
plt.savefig(outdir+"/sc_train.png")

# process st data
st = pd.read_csv(args.st_cnt,sep="\t",index_col=0)
st=st.transpose()
anno=pd.DataFrame(st.index.values,columns=['spots'],index=st.index.values)
anno['in_tissue']=1
gene_name=pd.DataFrame(st.columns.values,columns=['gene'],index=st.columns.values)
adata_vis1 = anndata.AnnData(X=st,obs= anno,var=gene_name)
adata_vis1.obs['sample'] = 'ST'
intersect = np.intersect1d(sc_adata.var_names, adata_vis1.var_names)
st_adata = adata_vis1[:, intersect].copy()
st_adata.layers["counts"] = st_adata.X.copy()
sc.pp.filter_genes(st_adata, min_counts=1)
sc.pp.filter_cells(st_adata,min_genes=3)

# fit st model
DestVI.setup_anndata(st_adata, layer="counts")
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.train(max_epochs=4000)

# save training figures
st_model.history["elbo_train"].plot()
plt.savefig(outdir+"/st_train.png")

# save results
st_adata.obsm["proportions"] = st_model.get_proportions()
st_adata.obsm["proportions"].to_csv(outdir+"/destvi.csv",sep='\t')
