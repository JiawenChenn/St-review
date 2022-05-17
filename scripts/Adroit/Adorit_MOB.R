#r/4.0.1

options(stringsAsFactors = F,check.names=F)
#devtools::install_github("TaoYang-dev/AdRoit/AdRoit")
library(data.table)
library(AdRoit)
library(parallel)
library(doParallel)
library(dplyr)

print(timestamp())

############################################################################################################

# Seqfish+ internal
###### Input path #####
scRNA_path = "./processed_data/MOB/ref/internal/ob.sc_cnt.hvg2000.1950cell.tsv"

scRNA_anno_path = "./processed_data/MOB/ref/internal/ob.sc.1950cell.mta.tsv"

output_dir = "int_hvg2000"

###### Input processing #####

scRna = fread(scRNA_path,sep="\t") %>% as.data.frame(check.names=F)
rownames(scRna) = scRna$V1
scRna=scRna[,-1]
scRna = scRna%>% t()  %>%  as.data.frame(check.names=F) # do not use t() when the input is gene*cnt

#row gene column cell
st_count = data.frame(fread("./processed_data/MOB/seqfish.cnt.genexrow.tsv"),check.names=F)
rownames(st_count) = st_count[,1]
st_count=st_count[,-1]


anno<-fread(scRNA_anno_path,header=TRUE)%>% as.data.frame(check.names=F)
rownames(anno) = anno[,1]
anno=anno[colnames(scRna),]
#head(anno)

anno = anno[,2]
#print(length(intersect(rownames(st_count),rownames(scRna))))
for(i in 1:ncol(scRna)){
    scRna[,i]=as.numeric(scRna[,i])
}
scRna =as.matrix(scRna)

###### scRNA model #####
single.ref = ref.build(scRna,
                       anno,
                       genes=intersect(rownames(st_count),rownames(scRna)),
                       samples = rep("MOB",ncol(scRna)),
                       normalize = "None",
                       multi.sample.bulk = T,
                       multi.sample.single = T,
                       nbootsids=5, 
                       minbootsize=50,
                       silent = F)

saveRDS(single.ref, file = paste0("./scripts/Adroit/result/MOB.Seqfish.",output_dir,".rds"))

print("single cell ref done")

###### st model #####
st_count_ma=as.matrix(st_count)
st_count_ma = st_count_ma[intersect(rownames(st_count_ma),rownames(scRna)),]

result=AdRoit.est(
       st_count_ma,
       single.ref,
       use.refvar = FALSE,
       per.sample.adapt = FALSE,
       silent = TRUE)

result=t(result%>% as.data.frame(check.names=F))%>% as.data.frame(check.names=F)

fwrite(result,paste0("./scripts/Adroit/result/MOB.Seqfish.",output_dir,".adroit.txt"),
col.names=T,row.names=T,sep="\t",quote=F)

print(timestamp())
