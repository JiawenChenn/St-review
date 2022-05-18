library(RCTD)
library(Matrix)
library(tidyverse)
library(data.table)

stdir = "/processed_data/MOB/"
stfile = "seqfish.cnt.genexrow.tsv"

scdir = "/processed_data/MOB/ref/internal/"
reffile = "OB.sc.cnt.tsv"
metafile = "OB.sc.mta.tsv"

resultdir = "/scripts/RCTD/MOB/"
resultfile = "MOB.Seqfish.own_all.rctd.txt"
if (!dir.exists(file.path(resultdir))) {
  dir.create(file.path(resultdir))
}

options(check.name = F) 
##### read in single cell and ST data ####
ref = fread(file.path(scdir, reffile))
ref = as.data.frame(ref, check.names = FALSE)
rownames(ref) = ref[,1]
ref[,1] = NULL

cluster = fread(file.path(scdir, metafile))
spots <- fread(file.path(stdir, stfile)) 
spots = as.data.frame(spots, check.names = FALSE)
rownames(spots) = spots[,1]
spots[,1] = NULL

# keep only overlap genes with minimum count > 0 
spots = spots[rowSums(spots) > 0,]
overlap_gene = intersect(rownames(spots), rownames(ref))
ref_overlap = ref[overlap_gene,]
spots = spots[overlap_gene,]

# clean the cluster labels
cell_types <- cluster$bio_celltype
names(cell_types) <- cluster$cell
cell_types[cell_types == "Neuron_M/TC"] = "Neuron_M.TC"

cell_types <- as.factor(cell_types)

nUMI <- colSums(ref_overlap)
names(nUMI) <- colnames(ref_overlap)

## create reference data object
reference <- Reference(ref_overlap, cell_types, nUMI)

#### Spatial RNA files ####
coords <- tibble(barcodes = colnames(spots))

#### for ST data
coords = coords %>% 
  separate(barcodes, into = c("xcoord", "ycoord"), sep = "x", remove = FALSE)

coords = as.data.frame(coords)
rownames(coords) <- coords$barcodes
coords$barcodes <- NULL
coords$xcoord = as.numeric(coords$xcoord)
coords$ycoord = as.numeric(coords$ycoord)
nUMI_spot <- colSums(spots)

### Create SpatialRNA object
puck <- SpatialRNA(coords, spots, nUMI_spot)

## create RCTD object
myRCTD <- create.RCTD(puck, reference, max_cores = 1, CELL_MIN_INSTANCE = 0)

#### run RCTD ####
myRCTD_full <- run.RCTD(myRCTD, doublet_mode = 'full')

#### Summarize results ####
RCTD_results_full <- myRCTD_full@results
RCTD_norm_weights = as.matrix(sweep(RCTD_results_full$weights, 1, rowSums(RCTD_results_full$weights), '/'))
RCTD_norm_weights_tibb = as_tibble(RCTD_norm_weights, rownames = "barcodes")
write_tsv(RCTD_norm_weights_tibb, file.path(resultdir, resultfile))

