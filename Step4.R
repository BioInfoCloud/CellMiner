###======可视化========
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(RColorBrewer)
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



