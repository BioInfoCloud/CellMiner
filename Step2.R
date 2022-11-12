


rm(list = ls())
library(impute)
library(limma)

#读取药物输入文件
drugDat <- read.table("output/drug.txt",sep="\t",header=T,check.names=F)
drugDat <- as.matrix(drugDat)
rownames(drugDat) <- drugDat[,1]
drug <- drugDat[,2:ncol(drugDat)]
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)



# 考虑到药物敏感性数据中存在部分NA缺失值，通过impute.knn()函数来评估并补齐药物数据。其中，impute.knn()函数是一个使用最近邻平均来估算缺少的表达式数据的函数。
mat <- impute.knn(data)
drug <- mat$data
drug <- avereps(drug) %>% t() %>% as.data.frame()

colnames(drug)[1:12]



# 读取表达输入文件
exp <- read.table("output/geneExp.txt", sep="\t", header=T, row.names = 1, check.names=F)
dim(exp)
exp[1:4, 1:4]


# 提取特定基因表达
library(WGCNA)
library(tidyr)
inputgene <- c("TP53","PTEN","BCAT2","EGFR","TMEM178A")
gl <- intersect(inputgene,row.names(exp))
exp <- exp[gl,] %>% t() %>% as.data.frame()

identical(rownames(exp),rownames(drug))


##======药物敏感性计算
corTab <-cor(drug,exp,method="pearson")
corPval <- corPvalueStudent(corTab,nSamples = nrow(drug))


##=========筛选有显著性的
##以EGFR为例

fitercor <- lapply(gl, function(g){
  index <- abs(corTab[,g])> 0.3 & corPval[,g] < 0.05
  df <- cbind(corTab[index,g],corPval[index,g])
  colnames(df) <- c("pearson","Pvalue")
  write.csv(df,file = paste0("output/",g,"-cor.csv"))
  df
})
length(fitercor) == length(gl)
names(fitercor) <- gl

save(fitercor,drug,exp,file = "output/fitercor.Rdata")











