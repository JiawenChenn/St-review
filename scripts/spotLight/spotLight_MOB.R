# r/3.6.0
#devtools::install_git("https://github.com/MarcElosua/SPOTlight")
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(tidyverse)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(R.matlab)
options(stringsAsFactors = F,check.names=F)

scRNA_path = "./processed_data/MOB/ref/internal/ob.sc_cnt.1950cell.tsv"

scRNA_anno_path = "./processed_data/MOB/ref/internal/ob.sc.1950cell.mta.tsv"

output_dir = "int_default"

###### Input processing #####

scRna = fread(scRNA_path,sep="\t") %>% as.data.frame(check.names=F)
rownames(scRna) = scRna$V1
scRna=scRna[,-1]
scRna = scRna%>% t()  %>%  as.data.frame(check.names=F) # do not use t() when the input is gene*cnt

cell_type = fread(scRNA_anno_path,header=TRUE)%>% as.data.frame(check.names=F)
cell_type=cell_type%>% filter(cell_type[,1] %in% colnames(scRna))
scRna = scRna[,as.character(cell_type[,1])]
sc = CreateSeuratObject(scRna, project = "SeuratProject", assay = "RNA",
  min.cells = 1, min.features = 0, names.field = 1, meta.data = NULL)

identical(as.character(cell_type$cell),sc@assays$RNA@counts@Dimnames[[2]])
sc@meta.data$subclass = cell_type[,2]
Seurat::Idents(object = sc) <- sc@meta.data$subclass
sc=subset(sc,nFeature_RNA>1)

cluster_markers_all <- Seurat::FindAllMarkers(object = sc, verbose = TRUE, only.pos = TRUE)
saveRDS(object = cluster_markers_all,
       file = here::here(paste0("./scripts/spotLight/out/MOB.seqfish.markers_sc.",output_dir,".RDS")))

set.seed(123)

st_count = data.frame(fread("./processed_data/MOB/seqfish.cnt.genexrow.tsv"),check.names=F)
st_count=na.omit(st_count)
rownames(st_count) = st_count[,1]
st_count=st_count[,-1]

for(j in 1:ncol(st_count)){
    st_count[,j] = as.numeric(st_count[,j])
}

##### start deconvolution #####

spotlight_ls <- spotlight_deconvolution(
  se_sc = sc,
  counts_spatial = st_count,
  clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 300, # number of cells per cell type to use
  hvg =500, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
  )

saveRDS(object = spotlight_ls, file = here::here(paste0("./scripts/spotLight/out/MOB.seqfish.",output_dir,".spotlight.rds")))

decon_mtrx <- spotlight_ls[[2]] %>% as.data.frame(check.names=F)
cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
decon_mtrx=decon_mtrx[,cell_types_all]
rownames(decon_mtrx)=colnames(st_count)

fwrite(decon_mtrx,paste0("./scripts/spotLight/out/MOB.seqfish.",output_dir,".spotlight.txt"),
col.names=T,row.names=T,sep="\t",quote=F)

print(timestamp())
