

> 作者：DoubleHelix
>
> 微信公众号：生物信息云

**CellMiner数据库，主要是通过国家癌症研究所癌症研究中心(NCI)所列出的60种癌细胞为基础而建立的。该数据库最初发表于2009年，后于2012年在Cancer Research杂志上进行了更新，题目为“CellMiner: a web-based suite of genomic and pharmacologic tools to explore transcript and drug patterns in the NCI-60 cell line set”。大家后期在使用该数据库记得应用相关文献。**

数据库地址：[CellMiner - Analysis Tools | Genomics and Pharmacology Facility (nih.gov)](https://discover.nci.nih.gov/cellminer/)

视频讲解：[微信](https://mp.weixin.qq.com/s/eR9wmWJzPj6Oih717qWToA)   [B站](https://www.bilibili.com/video/BV1kP411c7Rs/?vd_source=ad7486d1c0a79f7459d20781ce805fbc)

## 1.数据下载

[Download Data Sets](https://discover.nci.nih.gov/cellminer/loadDownload.do) 》》》》Processed Data Set: 

勾选RNA: RNA-seq和Compound activity: DTP NCI-60，点击下载即可

## 2.读入药物数据

```R

library(readxl)
dat1 <- read_excel(path = "data/DTP_NCI60_ZSCORE.xlsx", skip = 7)


colnames(dat1) <- dat1[1,]
dat1 <- dat1[-1,-c(67,68)]

# 筛选药物标准

table(dat1$`FDA status`)

# 选取经过临床试验（Clinical trial）和FDA批准（FDA approved）的药物结果
dat1 <- dat1[dat1$`FDA status` %in% c("FDA approved", "Clinical trial"),]
dat1 <- dat1[,-c(1, 3:6)]

ifelse(dir.exists("output"),FALSE,dir.create("output"))
write.table(dat1, file = "output/drug.txt",sep = "\t",row.names = F,quote = F)
```

## 3.读入表达数据

```R
###==============读入表达数据
dat2 <- read_excel(path = "data/RNA__RNA_seq_composite_expression.xls", skip = 9)
colnames(dat2) <- dat2[1,]
dat2 <- dat2[-1,-c(2:6)]

write.table(dat2, file = "output/geneExp.txt",sep = "\t",row.names = F,quote = F)
```

## 4.整理数据

```R
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

```

## 5.药敏相关性分析

```R
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
```

## 6.相关性拟合曲线

```R
###======可视化========
rm(list = ls())
library(ggplot2)
library(ggpubr)

load("output/fitercor.Rdata")
names(fitercor)


ifelse(dir.exists("opFig"),FALSE,dir.create("opFig"))
g <- names(fitercor)[1]
lapply(names(fitercor), function(g){
  data <- fitercor[[g]]
  data <- na.omit(data)
  if(!is.null(data)){
    #dr <- rownames(data)[1]
    for(dr in rownames(data)){
      df <- data.frame(exp = exp[,g],dr = drug[,dr])
      tit <- paste0("R:",round(data[dr,1],2),",p value = ",round(data[dr,2],3))
      
      p <- ggplot(data = df, aes(x = exp, y = dr)) + #数据映射
        geom_point(alpha = 0.6,shape = 19,size=3,color="#DC143C") +#散点图，alpha就是点的透明度
        #geom_abline()+
        labs(title = tit)+
        geom_smooth(method = lm, formula = y ~ x,aes(colour = "lm"), size = 1.2,se = T)+
        scale_color_manual(values = c("#808080")) + #手动调颜色c("#DC143C","#00008B", "#808080")
        theme_bw() +#设定主题
        theme(axis.title=element_text(size=15,face="plain",color="black"),
              axis.text = element_text(size=12,face="plain",color="black"),
              legend.position =  "none",
              panel.background = element_rect(fill = "transparent",colour = "black"),
              plot.background = element_blank(),
              plot.title = element_text(size=15, lineheight=.8,hjust=0.5, face="plain"),
              legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
        ylab(paste0("Activity z scores of ",dr)) + #expression的作用就是让log10的10下标
        xlab(paste0("The expression of ",g))
      ggsave(filename = paste0("opFig/",g,"-",dr,"-cor.pdf"),plot = p,width = 5,height = 5)
    }
  }

})

```

## 7.基因高低表达分组小提琴图

```R
###======可视化========
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(RColorBrewer)
load("output/fitercor.Rdata")
names(fitercor)


ifelse(dir.exists("opFig"),FALSE,dir.create("opFig"))
g <- names(fitercor)[1]
lapply(names(fitercor), function(g){
  data <- fitercor[[g]]
  data <- na.omit(data)
  if(!is.null(data)){
    #dr <- rownames(data)[1]
    for(dr in rownames(data)){
      df <- data.frame(exp = exp[,g],dr = drug[,dr])
      med <- median(df$exp)
      df$group <- ifelse(df$exp > med,"High","Low")
      head(df)

      p = ggplot(df, aes(group,dr,fill= group))+ 
        geom_violin(aes(fill = group),trim = FALSE)+
        geom_signif(comparisons = list(c("High","Low")),
                    step_increase = 0.1,
                    map_signif_level = T,
                    margin_top=0.2,
                    tip_length =0.02,
                    test = "t.test")+
        geom_boxplot(width = 0.1,fill = "white")+
        scale_fill_manual(values=c(brewer.pal(2,"Dark2")))+
        theme_classic()+
        labs(y= paste0("Activity z scores of ",dr),title= g)+
        theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
              plot.title = element_text(hjust = 0.5),
              axis.line=element_line(colour="black",size=0.25),
              axis.title.x = element_blank(),
              axis.text.x = element_text(face = "plain",colour = "black"),
              axis.text = element_text(size=12,face="plain",color="black"),
              legend.position="none"
        )
      p
      ggsave(filename = paste0("opFig/",g,"-",dr,"-violin.pdf"),plot = p,
             width = 3,height = 4)
    }
  }

})

```

