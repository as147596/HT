library(ggplot2)
library(ggside)
library(tidyverse)
library(ggforce)
library(patchwork)
data<-read.csv("data/test_cohort/test_16s/genus_16s.csv",row.names = 1)
meta<-read.csv("data/test_cohort/test_16s/SraRunTable.csv",row.names = 1)
meta$disease<-ifelse(meta$diastolic_bp>=90|meta$Systolic_BP>=140,"hypertension","healthy")
samples<-unique(meta$Sample.Name)
net_list<-list()

data<-data[,meta$Sample.Name]

genus<-rownames(data)
data$genus<-genus
data<-aggregate(.~genus,data,sum)
rownames(data)<-data$genus
data<-data[,-1]
genus_count<-apply(data,2,function(x){
  sum(x!=0)
})

#rownames(data)<-gsub("g__","",rownames(data))

ag<-read.csv("D:/master/project/ratio/agora.csv",row.names = 1)
ag_count<-apply(data[intersect(ag$x,rownames(data)),],2,function(x){
  ind<-which(x!=0)
  length(intersect(rownames(data)[ind],ag$x))
  #sum(x[ind])
})
ag_ratio<-apply(data[intersect(ag$x,rownames(data)),],2,function(x){
  ind<-which(x!=0)
  #length(intersect(rownames(data)[ind],ag$x))
  sum(x[ind])
})
ag_count/genus_count


bardata<-data.frame(sampleid=meta$Sample.Name,
                    ratio=ag_ratio,
                    group=c("AGORA Recall Species Percentage"))
bardata$group1<-meta$disease
d<-bardata[bardata$group1=="healthy",]
h<-bardata[bardata$group1=="hypertension",]

bardata$sampleid<-factor(bardata$sampleid,
                         levels = c(d$sampleid[order(d$ratio,decreasing = T)],
                                    h$sampleid[order(h$ratio,decreasing = T)]))



pd<-ggplot()+
  geom_bar(data=bardata,mapping=aes(x=sampleid,y=ratio,fill=group1),
           stat = "identity", alpha=0.6)+
  scale_x_discrete(breaks = c())+
  #geom_ysidedensity(data = bardata_d,mapping = aes(y=count,fill=group),alpha=0.6)+
  theme_bw()+
  scale_y_continuous(breaks = seq(0,100,20))+
  scale_fill_manual(name="",values = c("#5b94c2","#ea5b57"))+
  theme(#axis.text.x = element_text(angle = 45,hjust=1,size = 15),
    #axis.text.x = element_blank(),
    #axis.title = element_text(size = 18),
    axis.title = element_text(size = 14),
    legend.text = element_text(size=18),
    legend.title = element_text(size=20),
    plot.margin = margin(l=10,b=3,r=5,t=5),
    plot.title = element_text(size=20),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    strip.text = element_text(size = 20))+
  ylab("# AGORA2 Recall\nSpecies Percentage")+
  xlab("Samples")+
  geom_hline(yintercept = 80,linetype="dashed")

pd_side<-ggplot(bardata,aes(y=ratio,fill=group1))+
  geom_density(alpha=0.6)+
  scale_x_continuous(breaks = c(0.02,0.06))+
  theme_bw()+
  scale_fill_manual(name="",values = c("#5b94c2","#ea5b57"))+
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.ticks.y= element_blank(),
        plot.margin = margin(l=0.1,r=5)
  )
p<-wrap_plots(list(pd,pd_side),widths = c(0.8,0.1))
saveRDS(p,file = "output/p_filter_sample.rds")
keep<-bardata$ratio>80
keep_sample<-bardata$sampleid[keep]|>as.character()
write.csv(keep_sample,"data/res_grow_wd_pFBA/selected_sample.csv",row.names = F,quote = F)






















