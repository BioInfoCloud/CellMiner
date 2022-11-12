###======可视化========
rm(list = ls())
library(ggplot2)
library(ggpubr)

load("output/fitercor.Rdata")
names(fitercor)


ifelse(dir.exists("opFig"),FALSE,dir.create("opFig"))
# g <- names(fitercor)[1]
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



