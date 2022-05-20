# seqFISH+ MOB default gene set
library(STdeconvolve)

st = read.table("./processed_data/MOB/seqfish.cnt.genexrow.tsv", row.names = 1, header = T, check.names=FALSE)

## remove pixels with too few genes
st_counts = cleanCounts(counts = as.matrix(st),
                  min.lib.size = 100,
                  min.reads = 1,
                  min.detected = 1)
                  
## feature select for genes
corpus = restrictCorpus(st_counts,
                     removeAbove=1.0,
                     removeBelow = 0.02, 
                     alpha = 0.05,
                     plot = F,
                     verbose = TRUE)

## Note: the input corpus needs to be an integer count matrix of pixels x genes
corpus = corpus[,Matrix::colSums(corpus)!=0] 

pdf("Seqfish+ OB perplexity.pdf")
ldas = fitLDA(t(as.matrix(corpus)), Ks = 6:20, plot=T, verbose=TRUE)
dev.off()

## select model with minimum perplexity
optLDA = optimalModel(models = ldas, opt = "kneed") # can be min or a given input

## extract pixel cell-type proportions (theta) and cell-type gene expression profiles (beta) for the given dataset
results = getBetaTheta(optLDA, corpus = t(as.matrix(corpus)))

## deconvolution results
deconProp = results$theta
write.table(deconProp, "results/MOB.Seqfish.all.txt", quote=F)

# cell type annotation
# predicted transcriptional profiles of the deconvolved cell-types as the beta matrix
deconGexp = results$beta*1000

sc = read.table("./processed_data/MOB/ref/internal/OB.sc.cnt.tsv", row.names = 1, header = T, check.names=F)
label = read.table("./processed_data/MOB/ref/internal/OB.sc.mta.tsv", sep="\t", header=T)

mobProxyTheta2 = model.matrix(~ 0 + label$bio_celltype)
rownames(mobProxyTheta2) = colnames(sc)
library(stringr)
colnames(mobProxyTheta2) = substring(colnames(mobProxyTheta2), 19)
mobProxyGexp2 = sc %*% mobProxyTheta2

corMtx_beta = getCorrMtx(m1 = as.matrix(deconGexp), # the deconvolved cell-type `beta` (celltypes x genes)
                          m2 = t(as.matrix(mobProxyGexp2)), # the reference `beta` (celltypes x genes)
                          type = "b") # "b" = comparing beta matrices, "t" for thetas
write.table(corMtx_beta, "seqfish+_results/celltype_annotation.txt", quote=F, col.names=T, row.names=T)

pdf("results/celltype_annotation.pdf")
correlationPlot(mat = corMtx_beta,
                colLabs = "Deconvolved cell-types", # aka x-axis, and rows of matrix
                rowLabs = "Ground truth cell-types", # aka y-axis, and columns of matrix
                title = "Transcriptional correlation", annotation = TRUE)
dev.off()



