library(ggplot2)
library(tidyverse)
library(tidyr)
library(reshape2)
library(igraph)
library(limma)
library(dplyr)
library(ggpubr)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(ggthemes)
library(ggforce)
library(rstatix)
library(ALDEx2)
library(cowplot)
library(grid)
library(gground)


flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv")
flux<-flux[flux$taxon!="medium",]
flux$flux<-flux$flux*flux$abundance
selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
flux<-flux[flux$sample_id%in%selected$x,]
flux<-flux[,c(2,1,7,5)]

colnames(flux)[2:3]<-c("source","target")
index<-which(flux$flux<0)
flux[index,2:3]<-flux[index,3:2]
meta<-read.csv("data/test_cohort/test_16s/SraRunTable.csv",row.names = 1)
meta<-meta[meta$Sample.Name%in%selected$x,]
meta$disease<-ifelse(meta$diastolic_bp>=90|meta$Systolic_BP>=140,"hypertension","healthy")
flux<-flux[flux$sample_id%in%meta$Sample.Name,]
flux$flux<-abs(flux$flux)
#cluster1<-read.csv("output/cluster_1.csv")
res<-list()
topo<-data.frame()

for(i in meta$Sample.Name){
  sample_net<-flux[flux$sample_id==i,]
  #sample_net1<-sample_net[sample_net$target%in%cluster1$x|sample_net$source%in%cluster1$x,]
  sample_net_g<-graph_from_data_frame(sample_net[,2:4],directed = T)
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
  strength_in<-strength(sample_net_g,mode = "in",weights = E(sample_net_g)$flux)
  strength_out<-strength(sample_net_g,mode = "out",weights = E(sample_net_g)$flux)
  res[[i]]<-data.frame(degree_in=deg_in,degree_out=deg_out,betweenness=bet,
                       closeness=close,eigen_centrality=eigen,strength_in=strength_in,strength_out=strength_out)
}

boxdata_all<-data.frame(value=c(topo$density,topo$md),
                    var=rep(c("density","mean_distance"),each=nrow(topo)),
                    group=c(meta$disease,meta$disease))
ggplot(boxdata_all,aes(x=group,y=value,fill=group))+
  geom_violin(width=0.7,alpha=0.5)+
  geom_boxplot(width=0.2,alpha=0.8,size=1)+
  geom_point(size=2,alpha=0.2)+
  geom_line(aes(group=group),linewidth=1,alpha=0.1)+
  facet_wrap(.~var,nrow = 1,scales = "free_y")+
  scale_fill_manual(values = c(hypertension="#ea5b57",healthy="#5b94c2"))+#不显示箱式图的图例
  #scale_color_manual(values = c("#e54445","#5b94c2"))+
  xlab("")+
  #scale_y_continuous(limits = c(0,6))+
  theme_bw()+
  theme(axis.text = element_text(color = "black",face = "italic",size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 1),
        panel.grid.major = element_line(linewidth = 1),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 16),
        legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold.italic"),
        legend.text = element_text(size = 12),
        plot.margin = ggplot2::margin(b=0,t=10))+
  ggsignif::geom_signif(comparisons = list(c("healthy","hypertension")),
                        test = "wilcox.test",
                        size = 1)


res_metabolite<-lapply(res,function(x){
  index<-grep("_m$|\\[e\\]$",rownames(x))
  x[index,]
})

metab<-unique(unlist(lapply(res_metabolite,rownames)))

res_metabolite<-lapply(res_metabolite,function(x){
  tmp<-x[metab,]
  tmp[is.na(tmp)]<-0
  rownames(tmp)<-metab
  tmp
})




deg_inmat<-lapply(res_metabolite,function(x){
  x[,6]
})
deg_inmat<-as.data.frame(t(as.data.frame(deg_inmat)))
colnames(deg_inmat)<-metab


deg_outmat<-lapply(res_metabolite,function(x){
  x[,7]
})
deg_outmat<-as.data.frame(t(as.data.frame(deg_outmat)))
colnames(deg_outmat)<-metab

ml_g = lrn("regr.ranger", num.trees = 100, mtry = 2, min.node.size = 2, max.depth = 5)
ml_m = lrn("regr.ranger", num.trees = 100, mtry = 2, min.node.size = 2, max.depth = 5)
meta1<-meta[,c('AGE', 'Adiponectin', 'BMI', 'Body_Fat_Percentage', 'Glucose', 'Hb1Ac',
               'HDL', 'hs.CRP', 'Insulin', 'Total_Cholesterol', 'Triglycerides', 'VLDL', 
               'Waist_circumference', 'city', 'sex', 'smoker', 'Medicament_use')]
meta2<-model.matrix(~.-1,meta1)
double_bp<-function(mat,d_cols){
  res<-data.frame()
  for(i in 1:ncol(mat)){
    double_data<-data.frame(mat[,i,drop=F],y=meta[,d_cols],
                            meta2)
    double_data<-scale(double_data)|>as.data.frame()
    double_data<-na.omit(double_data)|>as.data.frame()
    model_data<-DoubleMLData$new(double_data,
                                 y_col = "y",
                                 d_cols = colnames(double_data)[1],
                                 x_cols = c('AGE','Adiponectin','BMI','Body_Fat_Percentage','Glucose','Hb1Ac',
                                            'HDL','hs.CRP','Insulin','Total_Cholesterol','Triglycerides',
                                            'VLDL','Waist_circumference','cityBarranquilla','cityBogota',
                                            'cityBucaramanga','cityCali','cityMedellin','sexmale','smokerYes',
                                            'Medicament_useYes'))
    set.seed(123)
    model = DoubleMLPLR$new(model_data, ml_g, ml_m)
    model$fit()
    tmp<-data.frame(ATE=model$coef,lower=model$confint()[1],
                    upper=model$confint()[2],pvalue=model$pval)
    res<-rbind(res,tmp)
  }
  res$metabolite<-rownames(res)
  res
}
keep<-apply(deg_inmat, 2, function(x){
  length(unique(x))>30
})
sbp_ind<-double_bp(mat = deg_inmat[,keep],d_cols = "Systolic_BP")
sbp_ind$fdr<-p.adjust(sbp_ind$pvalue,method = "fdr")
sbp_ind$y="Systolic_BP"
sbp_ind$x="Production\nStrength"
dbp_ind<-double_bp(mat = deg_inmat[,keep],d_cols = "diastolic_bp")
dbp_ind$fdr<-p.adjust(dbp_ind$pvalue,method = "fdr")
dbp_ind$y="Diastolic_BP"
dbp_ind$x="Production\nStrength"
keep<-apply(deg_outmat, 2, function(x){
  length(unique(x))>30
})
sbp_out<-double_bp(mat = deg_outmat[,keep],d_cols = "Systolic_BP")
sbp_out$fdr<-p.adjust(sbp_out$pvalue,method = "fdr")
sbp_out$y="Systolic_BP"
sbp_out$x="Consumption\nStrength"
dbp_out<-double_bp(mat = deg_outmat[,keep],d_cols = "diastolic_bp")
dbp_out$fdr<-p.adjust(dbp_out$pvalue,method = "fdr")
dbp_out$y="Diastolic_BP"
dbp_out$x="Consumption\nStrength"

sbp_all<-rbind(sbp_ind,sbp_out)
dbp_all<-rbind(dbp_ind,dbp_out)
bp_all<-rbind(sbp_all,dbp_all)
colors<- c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")
p_vol_all<-ggplot(bp_all,aes(x=ATE,y= -log10(fdr)))+
  facet_grid(x~y)+
  geom_point(aes(color = ATE, size = -log10(fdr)), alpha =0.5) +
  scale_color_gradientn(colors = colors) +
  theme_clean()+ggtitle("A.")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",size=1)+
  theme(legend.key = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = 'transparent'),
        axis.title = element_text(size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(r=10,b=10,t=10,l=10))+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = -2,hjust = -0.01,face = "bold"),
        strip.text = element_text(face = "bold",size = 16),
        legend.key.size = unit(0.5,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        plot.background = element_blank(),
        legend.text = element_text(size = 12))
p_vol_all
bp_all1<-bp_all
bp_all1$x<-gsub("\n","_",bp_all1$x)
write.table(bp_all1,file = "output/ATE of metabolic Strength.txt",quote = F,row.names = F,sep = "\t") 

anno<-read.csv("data/res_grow_wd_pFBA/annotations.csv")
quan<-function(x,anno){
  x$metabolite<-gsub("\\.e\\.","\\[e\\]",x$metabolite)
  anno1<-anno[anno$metabolite%in%x$metabolite,]
  tmp<-x$metabolite[x$fdr<0.05]
  anno1$group<-ifelse(anno1$metabolite%in%tmp,1,0)
  #q <- quantile(x$MW, probs = c(0.1, 0.9),na.rm=T)  
  #x<-na.omit(x)
  #x_filtered <- x[x >= q[1] & x <= q[2],]
  #x
  anno1
}

dbp_out1<-quan(dbp_out,anno)
dbp_ind1<-quan(dbp_ind,anno)
sbp_ind1<-quan(sbp_ind,anno)
sbp_out1<-quan(sbp_out,anno)

dbp_out1$x="Diastolic_BP"
dbp_ind1$x="Diastolic_BP"
dbp_out1$y="Consumption"
dbp_ind1$y="Production"
sbp_out1$x="Systolic_BP"
sbp_ind1$x="Systolic_BP"
sbp_out1$y="Consumption"
sbp_ind1$y="Production"

bp_all<-rbind(sbp_ind,sbp_out,dbp_ind,dbp_out)
bp_all$sign<-ifelse(bp_all$fdr<0.05,"Significant","Non\nSignificant")
#bp_all$MW<-anno$molecular_weight[match(bp_all$metabolite,anno$metabolite)]



ggplot(bp_all,aes(x=sign,y=MW,fill=sign))+
  geom_violin(width=0.7,alpha=0.5)+
  geom_boxplot(width=0.2,alpha=0.8,size=1)+
  geom_point(size=2,alpha=0.2)+
  scale_y_continuous(limits = c(0,800))+
  facet_wrap(.~x+y,nrow = 1,scales = "free_y")+
  #scale_fill_manual(values = c(hypertension="#ea5b57",healthy="#5b94c2"))+#不显示箱式图的图例
  #scale_color_manual(values = c("#e54445","#5b94c2"))+
  xlab("")+
  #scale_y_continuous(limits = c(0,6))+
  theme_bw()+
  theme(axis.text = element_text(color = "black",face = "italic",size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        panel.border = element_rect(linewidth = 1),
        panel.grid.major = element_line(linewidth = 1),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 16),
        #legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold.italic"),
        legend.text = element_text(size = 12),
        plot.margin = ggplot2::margin(b=0,t=10))+
  ggsignif::geom_signif(comparisons = list(c("Significant","Non\nSignificant")),
                        test = "wilcox.test",
                        size = 1)


pathway<-read.delim("data/ReactionDatabase.txt")
pathway<-pathway[pathway$Subsystem!="",]
subsys<-unique(pathway$Subsystem)
flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv")
background<-unique(gsub("^[0-9]+|_m$|\\[e\\]","",unique(flux$metabolite)))
source("script/enrichment_analysis.R")
pathwaylist1<-data.frame()
for(i in 1:length(subsys)){
  tmp<-pathway[pathway$Subsystem==subsys[i],]
  metabolite<-unlist(lapply(strsplit(tmp$Formula," -> | <=> | <=>|->"),function(x) x[1]))
  metabolite<-unlist(strsplit(metabolite," \\+ "))
  metabolite<-gsub("\\s*$","",unique(gsub("^[0-9.]+\\s*","",metabolite)))
  metabolite<-gsub("\\[.\\]","",metabolite)
  path_tmp<-data.frame(pathway=subsys[i],metabolite=unique(metabolite))
  pathwaylist1<-rbind(pathwaylist1,path_tmp)
}

strength_out<-rbind(sbp_out,dbp_out)

feeding_out_up<-gsub("\\.e\\.$","",strength_out$metabolite[strength_out$fdr<0.05&strength_out$ATE>0])
feeding_out_down<-gsub("\\.e\\.$","",strength_out$metabolite[strength_out$fdr<0.05&strength_out$ATE<0])

feeding_overpresentation_out_up<-enrichment_analysis(unique(feeding_out_up),background = background,pathway_data = pathwaylist1)
feeding_overpresentation_out_up$padj<-p.adjust(feeding_overpresentation_out_up$P_Value,method = "BH")
feeding_overpresentation_out_up<-feeding_overpresentation_out_up[feeding_overpresentation_out_up$padj<0.05,]
feeding_overpresentation_out_up$change<-"UP"
feeding_overpresentation_out_down<-enrichment_analysis(unique(feeding_out_down),background = background,pathway_data = pathwaylist1)
feeding_overpresentation_out_down$padj<-p.adjust(feeding_overpresentation_out_down$P_Value,method = "BH")
feeding_overpresentation_out_down<-feeding_overpresentation_out_down[feeding_overpresentation_out_down$padj<0.05,]
feeding_overpresentation_out_down$change<-"Down"

feeding_overpresentation_out_up$group<-"Out-degree Differences in the GMxN"

pathwaylist2<-data.frame()
for(i in 1:length(subsys)){
  tmp<-pathway[pathway$Subsystem==subsys[i],]
  metabolite<-unlist(lapply(strsplit(tmp$Formula," -> | <=> | <=>|->"),function(x) x[2]))
  metabolite<-unlist(strsplit(metabolite," \\+ "))
  metabolite<-gsub("\\s*$","",unique(gsub("^[0-9.]+\\s*","",metabolite)))
  metabolite<-gsub("\\[.\\]","",metabolite)
  path_tmp<-data.frame(pathway=subsys[i],metabolite=unique(metabolite))
  pathwaylist2<-rbind(pathwaylist2,path_tmp)
}
strength_in<-rbind(sbp_ind,dbp_ind)
feeding_in_up<-gsub("\\.e\\.$","",strength_in$metabolite[strength_in$fdr<0.05&strength_in$ATE>0])
feeding_in_down<-gsub("\\.e\\.$","",strength_in$metabolite[strength_in$fdr<0.05&strength_in$ATE<0])
feeding_overpresentation_in_up<-enrichment_analysis(unique(feeding_in_up),background = background,pathway_data = pathwaylist2)
feeding_overpresentation_in_up$padj<-p.adjust(feeding_overpresentation_in_up$P_Value,method = "BH")
feeding_overpresentation_in_up<-feeding_overpresentation_in_up[feeding_overpresentation_in_up$padj<0.05,]
feeding_overpresentation_in_up$change<-"UP"
feeding_overpresentation_in_up$group<-"In-degree Differences in the GMxN"

feeding_overpresentation_in_down<-enrichment_analysis(unique(feeding_in_down),background = background,pathway_data = pathwaylist2)
feeding_overpresentation_in_down$padj<-p.adjust(feeding_overpresentation_in_down$P_Value,method = "BH")
feeding_overpresentation_in_down<-feeding_overpresentation_in_down[feeding_overpresentation_in_down$padj<0.05,]
feeding_overpresentation_in_down$change<-"Down"
feeding_overpresentation_in_down$group<-"In-degree Differences in the GMxN"



feeding_overpresentation<-feeding_overpresentation_in_up[feeding_overpresentation_in_up$pathway%in%
                                                           feeding_overpresentation_out_up$pathway,]
feeding_overpresentation<-feeding_overpresentation[-grep("Transport",feeding_overpresentation$pathway),]
feeding_overpresentation$Ratio<-feeding_overpresentation$Overlap_Size/feeding_overpresentation$Pathway_Size
p_enrich_all<-ggplot(data = feeding_overpresentation,aes(reorder(pathway, Ratio), y = Ratio,fill=-log(padj)))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.7) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) +
  coord_flip() +
  labs(title = "",
       x = "",
       y = "Ratio") +theme(legend.key.size = unit(0.5,"cm"))+
  scale_fill_gradientn(values=seq(0,1,0.2),colors=c("#f9ed36","#f38466","#b81f25"),
                       limits = c(0, 10))+
  theme_cowplot()+
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(size=12,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.title.x = element_text(size = 14),
        #legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12)
  )+
  ggtitle("B.")+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        plot.title.position = "plot",
        plot.title = element_text(size=14,face = "bold",vjust = -4,hjust = 0.1),
        strip.text = element_text(face = "bold",size = 16),
        legend.title = element_text(size = 14,face = "bold"),
        plot.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.position = "top",
        plot.margin = margin(r=20))
pdf("output/fig4_cf.pdf",width = 12,height = 6)
grid.arrange(p_vol_all,p_enrich_all,nrow=1)
dev.off()
feeding_overpresentation<-feeding_overpresentation_in_up[feeding_overpresentation_in_up$pathway%in%
                                                           feeding_overpresentation_out_up$pathway,]
feeding_overpresentation<-feeding_overpresentation[-grep("Transport",feeding_overpresentation$pathway),]
feeding_overpresentation$Ratio<-feeding_overpresentation$Overlap_Size/feeding_overpresentation$Pathway_Size
feeding_overpresentation$category<-pathway$Subsystem_general[match(feeding_overpresentation$pathway,pathway$Subsystem)]

feeding_overpresentation_out<-feeding_overpresentation_out_up[feeding_overpresentation_out_up$pathway%in%
                                                                feeding_overpresentation_in_up$pathway,]
feeding_overpresentation_out<-feeding_overpresentation_out[-grep("Transport",feeding_overpresentation_out$pathway),]
feeding_overpresentation_out$Ratio<-feeding_overpresentation_out$Overlap_Size/feeding_overpresentation_out$Pathway_Size
feeding_overpresentation_out$category<-pathway$Subsystem_general[match(feeding_overpresentation_out$pathway,pathway$Subsystem)]

feeding_overpresentation$class<-"Production"
feeding_overpresentation_out$class<-"Consumption"
feeding_overpresentation<-rbind(feeding_overpresentation,feeding_overpresentation_out)
metabolite_name<-strsplit(feeding_overpresentation$Overlap,",")
metabolite_name<-sapply(metabolite_name,function(x){
  anno$name[match(paste0(x,"[e]"),anno$metabolite)]|>paste(collapse = ";")
})
feeding_overpresentation$metabolite_name<-metabolite_name
write.table(feeding_overpresentation,"output/enrichment.txt",row.names = F,quote = F,sep="\t")
order_path<-
  
  use_pathway <- group_by(feeding_overpresentation, category) %>%
  group_by(padj) %>%
  top_n(2, wt = Overlap_Size) %>%
  ungroup() %>%
  mutate(category = factor(category, 
                           levels = rev(c("Amino acid metabolism","Carbohydrate metabolism", "Central metabolism",
                                          "Inorganic nutrient metabolism","Polysaccharide metabolism", "Others")))) %>%
  dplyr::arrange(category, desc(padj)) %>%
  mutate(pathway = factor(pathway, levels = unique(pathway)))

width <- 0.5

xaxis_max <- max(-log10(use_pathway$padj)) + 1

rect.data <- group_by(use_pathway, category) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -2 * width,
    xmax = -1.5 * width,
    ymax = cumsum(n)/2,
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )


pal <-c('#7bc4e2', '#acd372', '#fbb05b', '#ed6ca4', '#74d4b3', '#c18be0')
names(pal)<-rev(c("Amino acid metabolism","Carbohydrate metabolism", "Central metabolism",
                  "Inorganic nutrient metabolism","Polysaccharide metabolism", "Others"))
p_enrich_all<-ggplot(use_pathway,aes(-log10(padj), y = pathway, fill = category)) +
  facet_grid(.~class)+
  geom_col(
    aes(y = pathway), width = 0.6, alpha = 0.8
  ) +labs(x="-Log(FDR)",title = "C.")+
  geom_text(
    aes(x = 0.05, label = pathway),
    hjust = 0, size = 5
  ) +
  geom_text(
    aes(x = 0.1, label = Overlap, colour = category), 
    hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic', 
    show.legend = FALSE
  ) +
  # 基因数量
  geom_point(
    aes(x = -0.4, size = Overlap_Size),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.4, label = Overlap_Size)
  ) +
  scale_size_continuous(name = 'Count', range = c(6, 12),breaks = c(6,8,10,12)) +
  # 分类标签
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = category),
    data = rect.data,
    radius = unit(2, 'mm'),
    inherit.aes = FALSE
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  labs(y = NULL) +
  scale_fill_manual(name = 'Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max, 2), 
    expand = expansion(c(0, 0))
  ) +
  theme_test() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    legend.position = c(0.9,0.5),
    strip.text=element_text(size=16,face = "bold"),
    #plot.title.position = "plot",
    legend.text = element_text(size = 12),
    plot.title = element_text(size=14,vjust = -1,hjust = 0,face = "bold")
  )



test<-data.frame(bp=c(meta$Systolic_BP,meta$diastolic_bp),
                 class=rep(c("Systolic_BP","Diastolic_BP"),each=312),
                 but=log(deg_inmat$`but[e]`+1))
but_p<-ggplot(test, aes(x = bp, y = but, color = class)) +
  annotate("rect",xmin = -Inf,xmax = 90,ymin = 0.5,ymax = 0.55,fill="green4",alpha=.5)+
  annotate("rect",xmin = 90,xmax = 99,ymin = 0.5,ymax = 0.55,fill="yellow4",alpha=.5)+
  annotate("rect",xmin = 99,xmax = 109,ymin = 0.5,ymax = 0.55,fill="#fbb05b",alpha=.5)+
  annotate("rect",xmin = 109,xmax = max(meta$diastolic_bp),ymin = 0.5,ymax = 0.55,fill="red3",alpha=.5)+
  annotate("rect",xmin = min(meta$Systolic_BP,na.rm = T),xmax = 140,ymin = 0.55,ymax = 0.6,fill="green4",alpha=.5)+
  annotate("rect",xmin = 140,xmax = 160,ymin = 0.55,ymax = 0.6,fill="yellow4",alpha=.5)+
  annotate("rect",xmin = 160,xmax = 180,ymin = 0.55,ymax = 0.6,fill="#fbb05b",alpha=.5)+
  annotate("rect",xmin = 180,xmax = Inf,ymin = 0.55,ymax = 0.6,fill="red3",alpha=.5)+
  annotate("rect",xmin = 150,xmax = 182,ymin = 0.42,ymax = 0.38,fill="green4",alpha=.5)+
  annotate("rect",xmin = 150,xmax = 182,ymin = 0.38,ymax = 0.34,fill="yellow4",alpha=.5)+
  annotate("rect",xmin = 150,xmax = 182,ymin = 0.34,ymax = 0.3,fill="#fbb05b",alpha=.5)+
  annotate("rect",xmin = 150,xmax = 182,ymin = 0.3,ymax = 0.26,fill="red3",alpha=.5)+
  annotate("text",x = 165,y = 0.4,label="healthy")+
  annotate("text",x = 165,y = 0.36,label="Grade1")+
  annotate("text",x = 165,y = 0.32,label="Grade2")+
  annotate("text",x = 165,y = 0.28,label="Grade3")+
  geom_point(alpha = 0.9) +
  scale_color_manual(name="",values = c('#acd372', '#fbb05b'))+
  geom_smooth(method = "loess",se = F) +
  theme_clean()+ylim(c(0,0.7))+labs(title = "D.",y="Butyrate production",x="BP")+
  theme(legend.key = element_rect(fill = 'transparent'),
        strip.text=element_text(size=16),
        axis.title = element_text(size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(r=10,b=10,t=10,l=10))+
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 14,face = "bold"),
        plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = -2,hjust = -0.015,face = "bold"),
        strip.text = element_text(face = "bold",size = 16),
        legend.key.size = unit(0.5,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        plot.background = element_blank(),legend.position = c(0.8,0.35),
        legend.text = element_text(size = 12))

but_p

sig_metab<-gsub(".e.$","[e]",bp_all$metabolite[bp_all$fdr<0.05e-5])|>unique()
sig_metab<-gsub("^X","",sig_metab)
flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv")
flux<-flux[flux$taxon!="medium",]
flux$flux<-flux$flux*flux$abundance
selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
flux<-flux[flux$sample_id%in%selected$x,]
flux<-flux[,c(2,1,7,5)]

#colnames(flux)[2:3]<-c("source","target")
meta<-read.csv("data/test_cohort/test_16s/SraRunTable.csv",row.names = 1)
meta<-meta[meta$Sample.Name%in%selected$x,]
meta$disease<-ifelse(meta$diastolic_bp>=90|meta$Systolic_BP>=140,"hypertension","healthy")
flux<-flux[flux$sample_id%in%meta$Sample.Name,]
flux<-flux %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux<-flux%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
flux<-flux[selected$x,]
mean_flux<-apply(flux,1,function(x){
  mean(abs(x),na.rm=T)
})

meta3<-cbind(meta2,mean_flux=mean_flux)|>as.data.frame()
fit<-lm(mean_flux~.,meta3)
mean_flux<-resid(fit)
alpha_df<-read.csv("output/alpha_df.csv",row.names = 1)
alpha_df<-alpha_df[rownames(flux),]
alpha_df<-prop.table(as.matrix(alpha_df),2)|>as.data.frame()
alpha_df<-reshape2::melt(as.matrix(alpha_df*100))
alpha_df$mean_flux<-mean_flux
Mload<-read.delim("data/MicrobialLoad_rdp_train_set_16_load.tsv",row.names = 1)
Mload<-Mload[rownames(flux),,drop=F]
Mload$mean_flux<-mean_flux
Mload$load<-Mload$load/2.5e11
Mload$Var2<-"Microbial load"


p_corline<-ggplot()+
  geom_point(data=alpha_df,mapping = aes(y=value,x=mean_flux,color=Var2),alpha = 0.5) +  # 添加散点图层，点的大小表示体重
  ggpmisc::stat_poly_line(data=alpha_df,mapping = aes(y=value,x=mean_flux,color=Var2),
                          formula = y ~ x,alpha=0.1) +  # 添加线性回归线
  stat_cor(data=alpha_df,mapping = aes(y=value,x=mean_flux,color=Var2),
           method = "spearman",size=4.5,label.x = 0.10,label.y = seq(0.25,0.1,len=6),alpha=1)+
  geom_point(data=Mload,mapping = aes(y=load,x=mean_flux,color=Var2),alpha = 0.5) +
  ggpmisc::stat_poly_line(data=Mload,mapping = aes(y=load,x=mean_flux,color=Var2),
                          formula = y ~ x,alpha=0.1) +  # 添加线性回归线
  stat_cor(data=Mload,mapping = aes(y=load,x=mean_flux,color=Var2),
           method = "spearman",size=4.5,label.x = 0.10,label.y = 0.07,alpha=1)+
  scale_y_continuous(    
    name = "scaled Alpha index",  # 主Y轴的名称  
    sec.axis = sec_axis(~.*2.5e11, name = "Microbial load",breaks = c(2e10,5e10,8e10,1.1e11)))+
  theme_bw()  +xlim(-0.2,0.4)+labs(x="Avg flux adj",title = "B.")+
  scale_color_manual(name="",values=(c('#7bc4e2', '#acd372', '#fbb05b', '#ed6ca4', '#74d4b3', '#c18be0','#35a9b1')))+
  theme(
    plot.title = element_text(size = 14, face = "bold",vjust = -2.5,hjust = 0.01),
    plot.title.position = "plot",
    legend.text = element_text(size = 12),
    axis.text = element_text(color = "black",size = 12),
    axis.title = element_text(size = 14,face = "bold")
  )

p_top<-grid.arrange(p_vol_all,p_corline,but_p,nrow=1,widths=c(4,5,2.8))
dev.off()
pdf("output/fig4_cf.pdf",width = 15.5,height = 11)
grid.arrange(p_top,p_enrich_all,nrow=2,heights=c(4,6))
dev.off()













hyper_sample<-meta$Sample.Name[meta$disease=="hypertension"]
deg_inmat_hy<-deg_inmat[hyper_sample,]
deg_outmat_hy<-deg_outmat[hyper_sample,]
wilcox.test(deg_inmat_hy$`but[e]`,deg_outmat_hy$`but[e]`,paired = T)






MGlcn<-deg_inmat[,grep("^MGlcn",colnames(deg_inmat))]
index<-rowSums(MGlcn)!=0
mglcn_sample<-rownames(MGlcn)[index]
bp_mglcn<-meta
bp_mglcn$group<-ifelse(bp_mglcn$Sample.Name%in%mglcn_sample,"MGlcn","Non MGlcn")
bp_mglcn<-na.omit(bp_mglcn)|>as.data.frame()
bp_mglcn<-bp_mglcn[,c('AGE', 'Adiponectin', 'BMI', 'Body_Fat_Percentage', 'Glucose', 'Hb1Ac',
                      'HDL', 'hs.CRP', 'Insulin', 'Total_Cholesterol', 'Triglycerides', 'VLDL', 
                      'Waist_circumference',"Systolic_BP","diastolic_bp","group")]
#bp_mglcn[,-ncol(bp_mglcn)]<-scale(bp_mglcn[,-ncol(bp_mglcn)])|>as.data.frame()

lon_data <- gather(bp_mglcn, key = "Variable", value = "value", -group)
lon_data$Variable[lon_data$Variable=="diastolic_bp"]<-"Diastolic_BP"
lon_data$value<-log(lon_data$value+1)
myt_test1 <- wilcox_test(group_by(lon_data, Variable), value ~ group) 
myt_test1<-myt_test1[order(myt_test1$p),]
lon_data$Variable<-factor(lon_data$Variable,levels = myt_test1$Variable)
myt_test1 <- wilcox_test(group_by(lon_data, Variable), value ~ group) 
myt_test11 <- add_significance(myt_test1, 'p') 
my.test <- add_xy_position(myt_test11, x = 'Variable', dodge = 0.8)
my.test$y.position<-6
p_bar_MG<-ggbarplot(lon_data, x = 'Variable', y = 'value',   
          fill = 'group', add = 'mean_sd',       
          color = 'gray30', position = position_dodge(0.6),  
          width = 0.6, size = 0.2, legend = 'right')+
  scale_fill_manual(values = c("#f9d580","#7ac7e2")) +
  labs(x = '', y = 'log(Value)', 
       fill = '')+ylim(0,6.5)+ggtitle("A.")+
  stat_pvalue_manual(my.test, label = 'p.signif', tip.length = 0.02)+
  theme(axis.text.x = element_text(size=12,angle = 30,hjust=1),
        axis.text.y=element_text(size=12),
        axis.title.y = element_text(size=16,face = "bold"),
        legend.text = element_text(size=14,face="bold"),
        plot.title.position = "plot",
        plot.title = element_text(size=14,face = "bold",vjust = 0.2,hjust = -0.01),
        plot.margin = margin(l=15))
  

flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv")
flux<-flux[flux$sample_id%in%selected$x,]
flux<-flux[flux$sample_id%in%meta$Sample.Name,]
flux<-flux[,c(2,7,5)]
flux<-flux[grep("_m$",flux$metabolite),]
flux_m<-flux %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux_m<-flux_m%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
flux_m[is.na(flux_m)]<-0
flux_m<-flux_m[meta$Sample.Name,]
meta$mg_group<-ifelse(meta$Sample.Name%in%mglcn_sample,1,0)
pvalue<-c()
logfc<-c()
for(i in 1:ncol(flux_m)){
  p<-wilcox.test(flux_m[meta$mg_group==1,i],flux_m[meta$mg_group==0,i])$p.val
  logfctmp<-mean(flux_m[meta$mg_group==1,i])-mean(flux_m[meta$mg_group==0,i])
  pvalue<-c(pvalue,p)
  logfc<-c(logfc,logfctmp)
}
dif_res<-data.frame(metabolite=colnames(flux_m),logfc=logfc,pvalue)
dif_res$fdr<-p.adjust(dif_res$pvalue,method = "fdr")


flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv")
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
res<-list()
for(i in meta$Sample.Name){
  sample_net<-flux[flux$sample_id==i,]
  sample_net_g<-graph_from_data_frame(sample_net[,2:4],directed = T)
  density<-edge_density(sample_net_g)
  eigen<- eigen_centrality(sample_net_g,directed=T,scale = T)$vector
  bet<-betweenness(sample_net_g,directed = T,normalized = T)
  deg_in<-igraph::degree(sample_net_g,mode="in",normalized = T)
  deg_out<-igraph::degree(sample_net_g,mode="out",normalized = T)
  
  close<-closeness(sample_net_g,normalized = T)
  strength_in<-strength(sample_net_g,mode = "in",weights = E(sample_net_g)$flux)
  strength_out<-strength(sample_net_g,mode = "out",weights = E(sample_net_g)$flux)
  res[[i]]<-data.frame(degree_in=deg_in,degree_out=deg_out,betweenness=bet,
                       closeness=close,eigen_centrality=eigen,strength_in=strength_in,strength_out=strength_out)
}

res_metabolite<-lapply(res,function(x){
  index<-grep("_m$|\\[e\\]$",rownames(x))
  x[index,]
})

metab<-unique(unlist(lapply(res_metabolite,rownames)))

res_metabolite<-lapply(res_metabolite,function(x){
  tmp<-x[metab,]
  tmp[is.na(tmp)]<-0
  rownames(tmp)<-metab
  tmp
})

deg_inmat<-lapply(res_metabolite,function(x){
  x[,6]
})
deg_inmat<-as.data.frame(t(as.data.frame(deg_inmat)))
colnames(deg_inmat)<-metab

deg_outmat<-lapply(res_metabolite,function(x){
  x[,7]
})

deg_outmat<-as.data.frame(t(as.data.frame(deg_outmat)))
colnames(deg_outmat)<-metab
meta$mg_group<-ifelse(meta$Sample.Name%in%mglcn_sample,1,0)
dif_mg<-function(mat){
  mat<-mat[meta$Sample.Name,]
  pvalue<-c()
  logfc<-c()
  for(i in 1:ncol(mat)){
    if(any(is.na(mat[meta$mg_group==1,i]))){
      pvalue<-c(pvalue,1)
      logfc=c(logfc,0)
    }
    else{
      p<-wilcox.test(mat[meta$mg_group==1,i],mat[meta$mg_group==0,i])$p.val
      logfctmp<-mean(mat[meta$mg_group==1,i])-mean(mat[meta$mg_group==0,i])
      pvalue<-c(pvalue,p)
      logfc<-c(logfc,logfctmp)
    }
  }
  dif_res<-data.frame(metabolite=colnames(mat),logfc=logfc,pvalue)
  dif_res$fdr<-p.adjust(dif_res$pvalue,method = "fdr")
  dif_res
}
mg_in<-dif_mg(deg_inmat)
mg_out<-dif_mg(deg_outmat)
mg_in$group<-"in"
mg_out$group<-"out"
mg_in$`-log10(fdr)`<- -log10(mg_in$fdr)
mg_out$`-log10(fdr)`<- log10(mg_out$fdr)
mg_all<-rbind(mg_in,mg_out)
mg_all$sign<-ifelse(mg_all$fdr<0.05&abs(mg_all$logfc)>0.5,"Ture","False")
mg_all$`-log10(fdr)`[mg_all$`-log10(fdr)`>10]<-10
mg_all$`-log10(fdr)`[mg_all$`-log10(fdr)`< -10]<- -10

in_table<-as.data.frame(table(mg_in$fdr<0.05&abs(mg_in$logfc)>0.5))|>column_to_rownames("Var1")|>t()
out_table<-as.data.frame(table(mg_out$fdr<0.05&abs(mg_out$logfc)>0.5))|>column_to_rownames("Var1")|>t()
colnames(in_table)<-colnames(out_table)<-c("Not\nSignificant","Significant")
rownames(in_table)<-rownames(out_table)<-""
p_vol<-ggplot(mg_all,aes(x=logfc,y= `-log10(fdr)`))+
  #facet_grid(x~y)+
  geom_point(aes(color = sign), alpha =0.5) +
  ylim(-10,10)+xlim(-6,6)+
  scale_y_continuous(breaks = seq(-10,10,2),labels = c(">=10",abs(seq(-8,8,2)),">=10"))+
  #scale_x_continuous(breaks = c(-3,-1.5,0,1.5,3,4.5,6))+
  scale_color_manual(name="Significant",values = c("black","red3")) +
  ggtitle("B.")+
  theme_few()+
  theme(legend.key = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = 'transparent'),
        strip.text=element_text(size=16),
        axis.title = element_text(size = 16),
        plot.title.position = "plot",
        plot.title = element_text(size=14,face = "bold",vjust = 0,hjust = -0.005))

tb_in<-gridExtra::tableGrob(in_table,theme = ttheme_default(base_size=10))
tb_out<-gridExtra::tableGrob(out_table,theme = ttheme_default(base_size=10))
p_vol_MG<-p_vol+annotation_custom(tb_in,xmin = -5.5,xmax = -1.5,ymin = 7,ymax = 9)+
  annotation_custom(tb_out,xmin = -5.5,xmax = -1.5,ymax = -7,ymin = -9)+
  annotate("rect",xmin = 5,xmax = 6,ymin=-10.5,ymax = 0,alpha = .2,fill="red4")+
  annotate("rect",xmin = 5,xmax = 6,ymin=0,ymax = 10.5,alpha = .2,fill="blue4")+
  annotate("text",x = 5.5,5,label="Production",angle=90,size=5)+
  annotate("text",x = 5.5,-5,label="Consumption",angle=90,size=5)+
  labs(y="-Log10(FDR)",x="LogFC")+
  geom_hline(yintercept = log10(0.05),linetype = "dashed",)+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed",)+
  geom_vline(xintercept = -0.5,linetype="dashed")+
  geom_vline(xintercept = 0.5,linetype="dashed")+
  theme(axis.title = element_text(size=14,face = "bold"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14,face = "bold"),
        plot.title.position = "plot",
        plot.title = element_text(size=14,face = "bold",vjust = -0,hjust = 0),
        legend.position = "top")



feeding_out_up<-gsub("\\[e\\]$","",mg_out$metabolite[mg_out$fdr<0.05&mg_out$logfc>0.5])
feeding_out_down<-gsub("\\[e\\]$","",mg_out$metabolite[mg_out$fdr<0.05&mg_out$logfc< -0.5])

feeding_overpresentation_out_up<-enrichment_analysis(unique(feeding_out_up),background = background,pathway_data = pathwaylist1)
feeding_overpresentation_out_up$padj<-p.adjust(feeding_overpresentation_out_up$P_Value,method = "BH")
feeding_overpresentation_out_up<-feeding_overpresentation_out_up[feeding_overpresentation_out_up$padj<0.05,]
feeding_overpresentation_out_up$change<-"UP"

feeding_in_up<-gsub("\\[e\\]$","",mg_in$metabolite[mg_in$fdr<0.05&mg_in$logfc>0.5])
feeding_in_down<-gsub("\\[e\\]$","",mg_in$metabolite[mg_in$fdr<0.05&mg_in$logfc< -0.5])
feeding_overpresentation_in_up<-enrichment_analysis(unique(feeding_in_up),background = background,pathway_data = pathwaylist2)
feeding_overpresentation_in_up$padj<-p.adjust(feeding_overpresentation_in_up$P_Value,method = "BH")
feeding_overpresentation_in_up<-feeding_overpresentation_in_up[feeding_overpresentation_in_up$padj<0.05,]
feeding_overpresentation_in_up$change<-"UP"
feeding_overpresentation_in_up$Ratio<-feeding_overpresentation_in_up$Overlap_Size/feeding_overpresentation_in_up$Pathway_Size
feeding_overpresentation_in_up<-feeding_overpresentation_in_up[-grep("Transport|Exchange",feeding_overpresentation_in_up$pathway),]

feeding_overpresentation_in_down<-enrichment_analysis(unique(feeding_in_down),background = background,pathway_data = pathwaylist2)
feeding_overpresentation_in_down$padj<-p.adjust(feeding_overpresentation_in_down$P_Value,method = "BH")
feeding_overpresentation_in_down<-feeding_overpresentation_in_down[feeding_overpresentation_in_down$padj<0.05,]
feeding_overpresentation_in_down$change<-"Down"
feeding_overpresentation_in_down$Ratio<-feeding_overpresentation_in_down$Overlap_Size/feeding_overpresentation_in_down$Pathway_Size
feeding_overpresentation_in_down<-feeding_overpresentation_in_down[feeding_overpresentation_in_down$Pathway_Size>=2&
                                                                     feeding_overpresentation_in_down$Ratio>=0.2&
                                                                     feeding_overpresentation_in_down$Overlap_Size>=2,]
sig_pathway<-rbind(feeding_overpresentation_in_up,feeding_overpresentation_in_down[1:5,])

sig_pathway$`-log10(padj)`<- -log10(sig_pathway$padj)
sig_pathway$`-log10(padj)`[sig_pathway$change=="Down"]<- -sig_pathway$`-log10(padj)`[sig_pathway$change=="Down"]

p_enrich_MG<-ggplot(data = sig_pathway,aes(reorder(pathway, `-log10(padj)`), y = `-log10(padj)`,fill=change))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.7) +
  #scale_y_continuous(labels = function(y) ifelse(y == 0, 0, ifelse(y<0,10^(-y),10^(y))),breaks = c(-log10(10),-log10(5),0,log10(5),log10(10),log10(20),log(40)))+
  #facet_wrap(~group,ncol = 1,scales = "free",space="free") +
  #facet_col(vars(group), scales = "free_y", space = "free")+
  scale_y_continuous(breaks = c(-5,0,5,10),labels = c(5,0,5,10))+
  coord_flip() +
  labs(title = "",
       x = "",
       y = "-Log10(FDR)") +theme(legend.key.size = unit(0.5,"cm"))+
  scale_fill_manual(name="Change",values = c(UP="firebrick",Down="dodgerblue4"))+
  theme_minimal()+
  theme(strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(size=12,face="bold"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12,face = "bold"),
        axis.title.x = element_text(size = 14),
        #legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12)
  )+
  geom_text(data = sig_pathway[sig_pathway$change == "UP",], aes(x = pathway, y = -0.01, label = pathway),
            hjust = 1, size = 4) +
  geom_text(data = sig_pathway[sig_pathway$change == "Down",], aes(x = pathway, y = 0.01, 
                                                                   label = str_wrap(pathway, width = 30)),
            hjust = 0, size = 4,lineheight=0.7)+expand_limits(y=c(-4,6))+
  ggtitle("C.")+theme(
    plot.title = element_text(size = 14, face = "bold",vjust = -0.5,hjust = 0.01),
    plot.title.position = "plot",legend.position = "top"
  )

tmp<-grid.arrange(p_vol_MG,p_enrich_MG,nrow=1,widths=c(6,5))
dev.off()
pdf("output/Fig_MG.pdf",width = 9.5,height = 7.5)
grid.arrange(p_bar_MG,tmp,nrow=2,heights=c(7,9))
dev.off()







bet_mat<-lapply(res_metabolite,function(x){
  x[,3]
})
bet_mat<-as.data.frame(t(as.data.frame(bet_mat)))
colnames(bet_mat)<-metab
rownames(meta)<-meta$Sample.Name

cor_ph<-function(mat,bp){
  keep<-apply(mat, 2, function(x){
    mean(x!=0)>0.1
  })
  mat<-mat[,keep]
  tmp<-cbind(mat,bp)
  r<-Hmisc::rcorr(as.matrix(tmp),type = "spearman")
  index<-which(colnames(r$r)%in%c("diastolic_bp","Systolic_BP"))
  P<-r$P[-index,index]
  R<-r$r[-index,index]
  net<-reshape2::melt(R)
  P<-reshape2::melt(P)
  net<-cbind(net,P[,3])
  colnames(net)<-c("metabolite","phenotype","R","pvalue")
  net#<-net[which(abs(net$R)>0.1&net$pvalue<0.01),]
}
bp<-meta[,c("diastolic_bp","Systolic_BP")]
ind_ph<-cor_ph(deg_inmat,bp)
ind_ph<-ind_ph[abs(ind_ph$R)>0.15&ind_ph$pvalue<0.05,]
out_ph<-cor_ph(deg_outmat,bp)
out_ph<-out_ph[abs(out_ph$R)>0.15&out_ph$pvalue<0.05,]
ind_ph$metabolite<-as.character(ind_ph$metabolite)
ind_ph$phenotype<-as.character(ind_ph$phenotype)
out_ph$metabolite<-paste0(as.character(out_ph$metabolite)," ")
out_ph$phenotype<-paste0(as.character(out_ph$phenotype)," ")



innode<-data.frame(node=as.character(c(unique(ind_ph$phenotype),unique(ind_ph$metabolite))),
                   nodetype=rep(c("phenotype","metabolite"),
                                c(length(unique(ind_ph$phenotype)),length(unique(ind_ph$metabolite)))),
                   class="InDegree")
outnode<-data.frame(node=c(unique(out_ph$phenotype),unique(out_ph$metabolite)),
                   nodetype=rep(c("phenotype","metabolite"),
                                c(length(unique(out_ph$phenotype)),length(unique(out_ph$metabolite)))),
                   class="OutDegree")
nodes<-rbind(innode,outnode)
edges<-rbind(ind_ph,out_ph)
cor_graph<-function(edges,nodes){
  cor_network<-graph_from_data_frame(edges,directed = F,vertices = nodes)
  V(cor_network)$nodes<-nodes$node
  ggraph(cor_network, layout = "dh") + #layout_tbl_graph_igraph unrooted stress
    geom_edge_link(aes(colour = R),
                   alpha = 1) +
    
    scale_edge_colour_gradientn(name = "Spearman\nCorrelation",
                                colors = c("aquamarine4","grey93","mediumpurple"),
                                limits = c(-0.25, 0.25),
                                space = "Lab",
                                na.value = "grey50", 
                                guide = "edge_colourbar") +
    
    geom_node_point(aes(
      colour = nodetype,
      stroke = 1.15),size=6,alpha=0.8) +
    
    scale_colour_manual(name = "node type",
                        values = c("metabolite"="#0f8096","phenotype"="#f09d69")
    )+
    geom_node_label(aes(label = nodes),max.overlaps = Inf,fill = scales::alpha("white", 0.6),
                    colour="black", repel = T,size=3) +
    #theme_void() +#ggtitle("C.")+
    guides(color = guide_legend(override.aes = list(size = 6)))+
    theme(plot.title = element_text(size = 16, face = "bold",vjust = 0,hjust = 0.01),
          #legend.key.size = unit(1,'cm'), 
          legend.title = element_text(size = 14,face = "bold"),
          legend.text = element_text(size = 12),
          panel.border = element_rect(color = "black", fill = NA),      # 面板边框
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.background = element_rect(color = "black", fill = "grey90"),  # facet 标题边框
          strip.text = element_text(color = "black",size = 16) ,
          plot.title.position = "plot",
          plot.margin = ggplot2::margin(t=0))+
    facet_graph(~class)
}






ind_cor_g<-cor_graph(edges = ind_ph,nodes = innode)+
  theme(legend.position = "none")
out_cor_g<-cor_graph(edges = out_ph,nodes = outnode)
patchwork::wrap_plots(list(ind_cor_g,out_cor_g))



