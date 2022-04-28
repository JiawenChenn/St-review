# Evaluation metric
library(energy)
library(data.table)
library(dplyr)
library(funfun)

evaluation_metric = function(truth,data){
  result_all=list()
  intersect_spot=intersect(rownames(truth),rownames(data))
  truth=truth[intersect_spot,]
  data=data[intersect_spot,]

  data=data.frame(data,check.names = F)
  truth=data.frame(truth,check.names=F)
  
  intersect_celltype=intersect(colnames(truth),colnames(data))
  truth=truth[,intersect_celltype]
  data=data[,intersect_celltype]
  
  mse=calc.mse(t(truth), t(data), rsq = FALSE)
  dcor=rep(0,ncol(data))
  diff_all = c()
  for(i in 1:ncol(data)){
    dcor[i]=as.numeric(dcor.test(truth[,i], data[,i], R=200)$statistic)
    diff_all=rbind(diff_all,data.frame(V1=rownames(truth),value=data[,i]-truth[,i],cell_type=colnames(data)[i]))
  }
  result_all[['mse']]=data.table(V1=rownames(truth),value=mse)
  result_all[['rmse']]=data.table(V1=rownames(truth),value=sqrt(mse))
  result_all[['dcor']]=data.table(cell_type=colnames(data),value=dcor)
  result_all[['diff']]=diff_all
  return(result_all)
}

# Simulate cell proportion and estimates
set.seed(1)
make_sim_proportion=function(){
data=matrix(runif(1000,min=0.01,max=0.99),nrow=100,ncol=10)
rownames(data)=paste0("spot",1:100);colnames(data)=paste0("celltype",1:10)
data=sweep(data,1,rowSums(data),FUN = "/")
}

truth=make_sim_proportion()
method1=make_sim_proportion()

result_method1=evaluation_metric(truth,method1)

# mse
head(result_method1[['rmse']])
# distance correlation
head(result_method1[['dcor']])
# difference
head(result_method1[['diff']])






