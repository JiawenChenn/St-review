# seqFISH+ MOB default gene set
library(reticulate)
library(Giotto)
options(stringsAsFactors = F,check.names=F)

# setting the environment
# use python/3.6.6 and unset PYTHONPATH
Sys.setenv(RETICULATE_PYTHON = "~/.local/share/r-miniconda/envs/giotto_env/bin/python")
use_python("~/.local/share/r-miniconda/envs/giotto_env/bin/python", required = T)

# script for running SpatialDWLS
run("MOB.Seqfish.own_all","./processed_data/MOB/seqfish.cnt.genexrow.tsv","./processed_data/MOB/ref/internal/OB.sc.cnt.tsv","./processed_data/MOB/ref/internal/OB.sc.mta.tsv")

run = function(dat, st_dir, sc_expression_dir, sc_metadata_dir, feat_min=10, min_feat=100) {
    
    ## Generate Giotto objects 
    st = createGiottoObject(expression = st_dir)
    
    cat("Dimension of ST data:", dim(st@expression[[1]]$rna$raw), "\n")
    
    ## filter
    st = filterGiotto(gobject = st,
                              expression_threshold = 1,
                              feat_det_in_min_cells = feat_min,
                              min_det_feats_per_cell = min_feat, # lower for few number of genes
                              expression_values = c('raw'),
                              verbose = T)
    
    ## normalize
    st = normalizeGiotto(gobject = st, scalefactor = 6000)
    
    ## add gene & cell statistics
    st = addStatistics(gobject = st, expression_values = 'raw')
    
    ## PCA ##
    st = calculateHVF(gobject = st)
    st = runPCA(gobject = st, center = TRUE, scale_unit = TRUE)
    
    ## sNN network (default)
    st = createNearestNetwork(gobject = st,
                               dim_reduction_to_use = 'pca', dim_reduction_name = 'pca',
                               dimensions_to_use = 1:10, k = 10)
    
    # Leiden clustering
    st = doLeidenCluster(gobject = st, resolution = 0.4, n_iterations = 1000)
    
    ## scRNAseq data
    sc_expression = paste0(datdir, sc_expression_dir)
    sc_metadata = paste0(datdir, sc_metadata_dir)
    
    giotto_SC = createGiottoObject(
      expression = sc_expression
    )
    giotto_SC = normalizeGiotto(giotto_SC)
    giotto_SC = addCellMetadata(giotto_SC, new_metadata = data.table::fread(sc_metadata,sep="\t"))
    
    markers_scran = findMarkers_one_vs_all(gobject=giotto_SC, method="scran",
                                           expression_values="normalized", cluster_column='bio_celltype', min_feats=3)
    
    top_markers = markers_scran[, head(.SD, 10), by="cluster"]
    celltypes = levels(factor(markers_scran$cluster))
    sign_list = list()
    for (i in 1:length(celltypes)){
      sign_list[[i]] = top_markers[which(top_markers$cluster == celltypes[i]),]$feats
    }
    
    # SpatialDWLS
    DWLS_matrix = makeSignMatrixDWLSfromMatrix(matrix = as.matrix(get_expression_values(giotto_SC,values = "normalized")),
                                              cell_type = pDataDT(giotto_SC)$bio_celltype,
                                              sign_gene = top_markers$feats)
    st = runDWLSDeconv(gobject = st, sign_matrix = DWLS_matrix)
    
    # results at testcombo@spatial_enrichment[[1]]$DWLS; spot by cell type (first column is cell ID)
    results = st@spatial_enrichment[[1]]$DWLS
    rowname = results$cell_ID
    results = as.data.frame(results[,-1])
    rownames(results) = rowname
    write.table(results, paste0("results/", dat, ".dwls.txt"), quote=F, row.names=T, col.names=T, sep="\t")
}