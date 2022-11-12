

# 读取药物相关数据
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

###==============读入表达数据
dat2 <- read_excel(path = "data/RNA__RNA_seq_composite_expression.xls", skip = 9)
colnames(dat2) <- dat2[1,]
dat2 <- dat2[-1,-c(2:6)]

write.table(dat2, file = "output/geneExp.txt",sep = "\t",row.names = F,quote = F)



