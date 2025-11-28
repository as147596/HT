library(tidyverse)
library(ggplot2)
library(reticulate)
library(cowplot)
library(patchwork)
#use_condaenv("F:/software/anaconda/envs/doubleml/")
#dml<-import("doubleml")
#np=import("numpy")
#pd=import("pandas")
#sklearn=import("sklearn")
#RandomForestRegressor=sklearn$ensemble$RandomForestRegressor
flux<-read.csv("D:/master/micom/res_grow_wd/exchanges.csv",header = T)
flux<-flux[,c(2,7,5)]
meta<-read.csv("data/meta.csv",row.names = 1)
meta$sample_id<-gsub("-|\\.","_",meta$sample_id)

flux<-flux %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux<-flux%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
flux[is.na(flux)]<-0

colnames(flux)<-paste("m",colnames(flux),sep = "_")
colnames(flux)<-gsub("\\[e]$","_e",colnames(flux))


flux<-flux[meta$sample_id,]

flux<-flux[rownames(flux)!="SZAXPI015258_57",]
keep<-apply(flux,2,function(x){
  mean(x!=0)>0.1&length(unique(x))>3
})

flux<-flux[,keep]
flux<-scale(flux)|>as.data.frame()
flux$y<-as.numeric(as.factor(meta$disease[meta$sample_id!="SZAXPI015258_57"]))-1
flux<-flux[,grep("_m$",colnames(flux))]
predictors<-colnames(flux)[-ncol(flux)]
meta<-meta[meta$sample_id%in%rownames(flux),]

#group<-meta$study_condition|>factor(levels = c("control","pre-hypertension","hypertension"))|>as.numeric()
#flux$group<-group
#cf_y = 0.1
#cf_d = 0.1
#theta = 5.0
#set.seed(123)
#flux1<-flux[,c('m_26dap_M_m','group',"m_oxa_m")]
#dml_data = dml$DoubleMLData(flux, 'group', 'm_26dap_M_m')
#
#set.seed(123)
#dml_obj = dml$DoubleMLPLR(dml_data,
#                          ml_l=RandomForestRegressor(),
#                          ml_m=RandomForestRegressor(),
#                          n_folds=as.integer(5),
#                          score='partialling out')
#dml_obj$fit()
#print(dml_obj)
#
#dml_obj$sensitivity_analysis()
#print(dml_obj$sensitivity_summary)
#
#dml_obj$sensitivity_params
#
#dml_obj$sensitivity_plot()

ATE<-read.csv("output/doubleml_systolic.csv")
ATE$padj<-p.adjust(ATE$pvalue,method = "BH")
sign<-sign(ATE[,c("ATE","CI2.5","CI97.5","theta.lower","theta.upper")])
keep<-apply(sign, 1, function(x){
  length(unique(x))==1
})
sig<-ATE[keep&ATE$padj<0.005&ATE$RV>0.05,]
sig$color<-ifelse(sig$ATE>0,"R","P")
sig$metabolite<-factor(sig$X,levels = sig$X[order(sig$ATE)])
sig_SBP<-sig

ATE<-read.csv("output/doubleml_diastolic.csv")
ATE$padj<-p.adjust(ATE$pvalue,method = "BH")
sign<-sign(ATE[,c("ATE","CI2.5","CI97.5","theta.lower","theta.upper")])
keep<-apply(sign, 1, function(x){
  length(unique(x))==1
})
sig<-ATE[keep&ATE$padj<0.005&ATE$RV>0.05,]
sig$color<-ifelse(sig$ATE>0,"R","P")
sig$metabolite<-factor(sig$X,levels = sig$X[order(sig$ATE)])
sig_DBP<-sig

sig_SBP$class<-"SBP"
sig_DBP$class<-"DBP"
sig<-rbind(sig_DBP,sig_SBP)
#sig$metabolite<-factor(sig$metabolite,levels = sig$metabolite[order(sig$ATE)])

sig$x=ifelse(sig$ATE<0,0.05,min(sig$CI2.5))
sig$padj1<-ifelse(sig$padj<0.001,"<0.001",round(sig$padj,3))
sig$label<-paste(sprintf("%.3f", sig$RV),sig$padj1,sep = "      ")

p1=ggplot(data=sig,mapping = aes(y=reorder(metabolite,ATE)))+
  geom_text(aes(x=0,label=label),hjust=0)+
  xlim(-0.01,0.05)+xlab("")+ylab("")+
  facet_grid(class~., scale='free',space = 'free')+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.margin = ggplot2::margin(r=-80,t=10),
    plot.subtitle = element_text(hjust = 0.92,size = 12,face="bold"),
    plot.title.position = "plot",
    plot.title = element_text(size=14,vjust = 2,,face = "bold"),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
  )+
  labs(subtitle = "RV        Padj",title = "B.")
p1

p2=ggplot(data=sig, aes(x = ATE, y = reorder(metabolite,ATE))) +
  #annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,ymax=14.5,fill="green4",alpha=0.05)+
  #annotate("rect",xmin=-Inf,xmax=Inf,ymin=14.5,ymax=25,fill="red4",alpha=0.05)+
  geom_point(data=sig, aes(x = ATE, y = metabolite, color = color),size=3) +
  scale_color_manual(values = c("blue4","red4"))+
  geom_errorbar(data=sig,mapping = aes(xmin = theta.lower, xmax = theta.upper,y=metabolite),color="black",width=0.3)  +
  geom_vline(xintercept = 0, color = "grey") +
  theme_test()+ ylab("") +
  xlab("ATE") +
  theme(legend.position = "none",
        text = element_text(size = 16),axis.title.x = element_text(size = 14,face = "bold"),
        axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        plot.margin = ggplot2::margin(10, 10, 10, -30),
        strip.text = element_text(face = "bold",size = 16),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))+
  #geom_rect(data=sig[sig$ATE>0,],mapping = aes(xmin=x-0.01,xmax=x+0.04,y=metabolite))
  #geom_text(data=sig[sig$ATE>0,],mapping=aes(x=x,y=metabolite,label=label),hjust=0)+
  #geom_text(data=sig[sig$ATE<0,],mapping=aes(x=x,y=metabolite,label=label),hjust=0)+
  facet_grid(class~., scale='free',space = 'free')
p_double<-wrap_plots(p1,p2,nrow = 1,widths = c(0.3,0.5))
saveRDS(p_double,file = "output/doubleml_plot.rds")
ggsave("output/doubleml_plot.pdf",width = 8,height = 6)
anno<-read.csv("data/res_grow_wd_pFBA/annotations.csv",header = T)
anno_select<-anno[match(sig$X,anno$metabolite),]
sig<-cbind(sig,anno_select)
write.csv(sig,"output/doubleml_res.csv")


filter<-readRDS("output/p_filter_sample.rds")
filter[[1]]=filter[[1]]+
  ggtitle("A.")+theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size=14,face = "bold"),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.key.size = unit(0.5,"cm"),
    strip.text = element_blank(),
    plot.margin = ggplot2::margin(r=-80),
    plot.subtitle = element_text(hjust = 0.9,size = 12,face="bold"),
    plot.title.position = "plot",
    plot.title = element_text(size=14,vjust = 0,face = "bold"),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )
filter  

flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv")
selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
flux<-flux[flux$sample_id%in%selected$x,]
flux$flux<-flux$flux*flux$abundance
flux<-flux[,c(2,1,7,5)]

colnames(flux)[2:3]<-c("source","target")
index<-which(flux$flux<0)
flux[index,2:3]<-flux[index,3:2]
meta<-read.csv("data/test_cohort/test_16s/SraRunTable.csv",row.names = 1)
meta<-meta[meta$Sample.Name%in%selected$x,]
meta$disease<-ifelse(meta$diastolic_bp>=90|meta$Systolic_BP>=140,"hypertension","healthy")
flux<-flux[flux$sample_id%in%meta$Sample.Name,]
flux$flux<-abs(flux$flux)
importance<-read.csv("output/doubleml_res.csv",row.names = 1)
#importance<-importance[importance$ATE<0,]
importance<-gsub("_m$","[e]",importance$X)
res<-list()
topo<-data.frame()

for(i in meta$Sample.Name){
  sample_net<-flux[flux$sample_id==i,]
  sample_net1<-sample_net[sample_net$target%in%importance|sample_net$source%in%importance,]
  sample_net_g<-graph_from_data_frame(sample_net1[,2:4],directed = T)
  density<-edge_density(sample_net_g)
  diam<-diameter(sample_net_g,directed = T)
  md<-mean_distance(sample_net_g,directed = T,weights = log10((1/E(sample_net_g)$flux)+1))
  CC<-transitivity(sample_net_g,type = "global")
  comp<-components(sample_net_g)$no
  topo<-rbind(topo,data.frame(density=density,diam=diam,md=md,CC=CC,comp=comp,meanwei=sum(E(sample_net_g)$flux)))
  eigen<- eigen_centrality(sample_net_g,directed=T,scale = T)$vector
  bet<-betweenness(sample_net_g,directed = T,normalized = T)
  deg_in<-igraph::degree(sample_net_g,mode="in",normalized = T)
  deg_out<-igraph::degree(sample_net_g,mode="out",normalized = T)
  close<-closeness(sample_net_g,normalized = T)
  res[[i]]<-data.frame(degree_in=deg_in,degree_out=deg_out,betweenness=bet,closeness=close,eigen_centrality=eigen)
}

boxdata<-data.frame(value=c(topo$density,topo$md),
                    var=rep(c("density","mean\ndistance"),each=nrow(topo)),
                    group=c(meta$disease,meta$disease))

p_box_metabolite<-ggplot(boxdata,aes(x=group,y=value,fill=group))+
  geom_violin(width=0.7,alpha=0.5)+
  geom_boxplot(width=0.2,alpha=0.8,size=1)+
  #geom_point(size=2,alpha=0.2)+
  geom_line(aes(group=group),linewidth=1,alpha=0.1)+
  facet_wrap(.~var,nrow = 1,scales = "free_y")+
  scale_fill_manual(values = c(hypertension="#ea5b57",healthy="#5b94c2"))+#不显示箱式图的图例
  #scale_color_manual(values = c("#e54445","#5b94c2"))+
  xlab("")+ggtitle("C.")+
  #scale_y_continuous(limits = c(0.05,0.12))+
  theme_bw()+
  scale_y_continuous(breaks = c(0.05,0.08,0.11,5,10,15))+
  ggsignif::geom_signif(data = boxdata[boxdata$var=="density",],comparisons = list(c("healthy","hypertension")),
                        test = "wilcox.test",
                        size = 1,y_position = 0.1)+
  ggsignif::geom_signif(data = boxdata[boxdata$var=="mean\ndistance",],comparisons = list(c("healthy","hypertension")),
                        test = "wilcox.test",
                        size = 1,y_position = 15.5)+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title.y = element_text(size = 14,face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none",
        plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = 0,face = "bold"),
        panel.border = element_rect(linewidth = 1),
        panel.grid.major = element_line(linewidth = 1),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 16),
        legend.key.size = unit(0.8,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = ggplot2::margin(b=0,t=10))
tmp<-wrap_plots(filter,p_box_metabolite,ncol=1,heights = c(7,5))
dev.off()
pdf("output/fig2.pdf",width = 11,height = 7)
grid.arrange(plot_grid(tmp),plot_grid(p_double),nrow = 1,widths = c(4.5,6))
dev.off()
