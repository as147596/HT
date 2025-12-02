library(NetCoMi)
library(tidyverse)
library(tidyr)
library(reshape2)
library(ggraph)
library(igraph)
library(colorspace)
library(RColorBrewer)
library(graphlayouts)
library(ggvenn)
library(data.table)
library(fgsea)
library(grDevices)
library(ggheatmap)
library(gggsea)
library(ggrepel)
library(showtext)
library(ggforce)
library(clusterProfiler)
library(enrichplot)
library(GseaVis)
library(aplot)
library(patchwork)
library(gridExtra)
library(cowplot)
library(curatedMetagenomicData)
library(ggpubr)
#showtext_auto()
ab<-read.csv("data/discovery_cohort/16s/genus_16s.csv",row.names = 1,header = T)
#samplename<-as.character(read.csv("data/abundance.csv",row.names = 1,header = F)[1,])
#colnames(ab)<-samplename
##taxa<-strsplit(rownames(ab),"\\|")
#genus<-unlist(lapply(taxa, function(x){x[6]}))
#genus<-gsub("g__","",genus)
#ab<-data.frame(genus=genus,ab)
#ab_genus<-aggregate(.~genus,data=ab,sum)
#rownames(ab_genus)<-ab_genus$genus
#ab_genus<-ab_genus[,-1]
meta<-read.csv("data/discovery_cohort/16s/SraRunTable.csv")
selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
meta<-meta[meta$Sample.Name%in%selected$x,]
meta$sample_id<-gsub("-","\\.",meta$Sample.Name)
ab_genus<-ab
ab_genus<-ab_genus[,meta$sample_id]

keep<-apply(ab_genus, 1, function(x){
  mean(x!=0)>0.1
})
ab_genus<-ab_genus[keep,]
meta$group<-ifelse(meta$diastolic_bp>=90|meta$Systolic_BP>=140,"hypertension","healthy") ##WHO
net <- netConstruct(data = t(ab_genus)[meta$group=="hypertension",],data2 = t(ab_genus)[meta$group=="healthy",], 
                    filtTax = "none",
                    #filtTaxPar = list(relFreq = 0),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta$group=="healthy")),
                    measure = "sparcc",
                    #measurePar = list(nlambda=10, 
                    #                  rep.num=10,
                    #                  Rmethod = "approx"),
                    normMethod = "clr", 
                    zeroMethod = "alrEM",
                    sparsMethod = "threshold",
                    thresh = 0,
                    dissFunc = "signed",
                    verbose = 3,
                    seed = 123456)
thresh=(sum(abs(net[["assoMat1"]])+abs(net[["assoMat2"]]))-2*nrow(net[["assoMat2"]]))/(2*nrow(net[["assoMat2"]])*(nrow(net[["assoMat2"]])-1))
net <- netConstruct(data = t(ab_genus)[meta$group=="hypertension",],data2 = t(ab_genus)[meta$group=="healthy",], 
                    filtTax = "none",
                    #filtTaxPar = list(relFreq = 0.1),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta$group=="healthy")),
                    measure = "sparcc",
                    #measurePar = list(nlambda=10, 
                    #                  rep.num=10,
                    #                  Rmethod = "approx"),
                    normMethod = "clr", 
                    zeroMethod = "alrEM",
                    sparsMethod = "threshold",
                    thresh = thresh,
                    dissFunc = "signed",
                    verbose = 3,
                    seed = 123456)

props <- netAnalyze(net, 
                    centrLCC = FALSE,
                    avDissIgnoreInf = TRUE,
                    sPathNorm = FALSE,
                    clustMethod = "cluster_fast_greedy",
                    hubPar = c("degree", "eigenvector"),
                    hubQuant = 0.9,
                    lnormFit = TRUE,
                    normDeg = FALSE,
                    normBetw = FALSE,
                    normClose = FALSE,
                    normEigen = FALSE)

comp <- netCompare(props, 
                   permTest = FALSE, 
                   verbose = FALSE,
                   seed = 123456)

summ_res<-summary(comp, 
                  groupNames = c("hypertension","healthy"),pAdjust=T,
                  #showCentr = c("degree", "between", "closeness"), 
                  showCentr="all",
                  numbNodes = 200)

diff <- diffnet(net,discordThresh=0.8,alpha = 0.05,lfdrThresh = 0.2,
                diffMethod = "permute", nPerm = 1000,cores = 14,
                storeCountsPerm = TRUE, 
                fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
                storeAssoPerm = TRUE,
                fileStoreAssoPerm = "assoPerm",
                adjust = "none",seed = 123)
saveRDS(diff,"data/specie_network_res1.rds")
diff<-readRDS("data/specie_network_res1.rds")

props_pears <- netAnalyze(net, 
                          clustMethod = "cluster_fast_greedy",
                          weightDeg = TRUE,
                          normDeg = FALSE,
                          gcmHeat = FALSE)


diffmat_sums <- rowSums(diff$diffAdjustMat)
diff_asso_names <- names(diffmat_sums[diffmat_sums > 0])


diffedge<-diff$assoMat1-diff$assoMat2
diffpadj<-diff[["pAdjustMat"]]
diffpadj[upper.tri(diffpadj)]<-1

diffpadj[diffpadj>0.05]<-1
diffedge<-reshape2::melt(diffedge)
colnames(diffedge)<-c("m1","m2","diff")
diffpadj<-reshape2::melt(diffpadj)
colnames(diffpadj)<-c("m1","m2","p")
diffpadj<-diffpadj[diffpadj$p<0.05&diffpadj$m1!=diffpadj$m2,]
diffedge<-diffedge[paste(diffedge$m1,diffedge$m2,sep = ",")%in%paste(diffpadj$m1,diffpadj$m2,sep = ","),]
difnetwork<-data.frame(Var1=diffedge$m1,Var2=diffedge$m2,cor=diffedge$diff)

nodes<-unique(c(difnetwork$Var1,difnetwork$Var2))


nodeatr<-data.frame(nodes = nodes,nodetype="difnetwork")
nodeatr<-nodeatr[nodeatr$nodes%in%unique(c(difnetwork$Var1,difnetwork$Var2)),]
importance<-read.csv("output/doubleml_res.csv")
importance<-gsub("^m_","",importance$X)
importance<-gsub("_m$","_e",importance)
flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv",header = T)
flux<-flux[flux$sample_id%in%selected$x,]
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
flux<-flux[flux$metabolite%in%importance,c(1,2,7)]
imp_species<-unique(flux$taxon)
spedif_nodes<-nodeatr
spedif_nodes$nodetype="Related"
spedif_nodes$nodetype[spedif_nodes$nodes%in%imp_species]="Importance"
dif_graph<-graph_from_data_frame(d=difnetwork,vertices = spedif_nodes,directed = F)
set.seed(444)
cluster<-cluster_louvain(dif_graph)
#cluster<-cluster_walktrap(dif_graph)
cluster<-membership(cluster)
table(cluster)

V(dif_graph)$cluster<-paste("cluster",cluster,sep = "_")
tmp<-V(dif_graph)$cluster
tmp[tmp=="cluster_3"]<-"cluster_3(FERM)"
V(dif_graph)$cluster<-factor(tmp,levels = unique(tmp))
g <- sample_islands(length(table(cluster)), max(table(cluster)), 0.5, length(table(cluster)))
g <- igraph::simplify(g)
set.seed(123)
bb <- layout_as_backbone(g, keep = 0.2)
pos_tmp<-bb$xy
pos<-data.frame()
nodelist<-V(dif_graph)$cluster%>%unique()
spedif_nodes$cluster<-V(dif_graph)$cluster
for(i in 1:length(table(cluster))){
  set.seed(123)
  x.pos<-sample((-7.5+i):(-6.5+i))
  tmp<-pos_tmp[((max(table(cluster))*(i-1))+1):(max(table(cluster))*i),]
  tmp<-tmp[sample(1:max(table(cluster)),sum(spedif_nodes$cluster==nodelist[i])),,drop=F]
  pos<-rbind(pos,tmp)
}
rownames(pos)<-names(V(dif_graph))[order(V(dif_graph)$cluster)]
pos<-pos[names(V(dif_graph)),]
difgraph<-ggraph(dif_graph,
                 layout = "manual",
                 x = pos[, 1],
                 y = pos[, 2]) +
  geom_edge_link0(aes(col = cor), width = 0.2) +
  geom_node_point(aes(fill = cluster), shape = 21, size = 3) +
  geom_mark_hull( 
    aes(x, y, group = cluster, fill = cluster),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("cluster_1"="orange3","cluster_2"="orchid4","cluster_3(FERM)"="red","cluster_4"="orangered4",
                                "cluster_5"="lightsteelblue4","cluster_6"="linen","cluster_7"="orange","cluster_8"="royalblue4",
                                "cluster_9"="skyblue4","cluster_10"="thistle","cluster_11"="turquoise4","cluster_12"="palegreen4","cluster_13"="mediumseagreen"),
                     labels=c(paste("guild",1:2,sep="_"),"guild_3(FERM)",paste("guild",4:13,sep="_"))) +
  scale_fill_manual(values = c("cluster_1"="orange3","cluster_2"="orchid4","cluster_3(FERM)"="red","cluster_4"="orangered4",
                               "cluster_5"="lightsteelblue4","cluster_6"="linen","cluster_7"="orange","cluster_8"="royalblue4",
                               "cluster_9"="skyblue4","cluster_10"="thistle","cluster_11"="turquoise4","cluster_12"="palegreen4","cluster_13"="mediumseagreen"),
                    labels=c(paste("guild",1:2,sep="_"),"guild_3(FERM)",paste("guild",4:13,sep="_"))) +
  scale_edge_colour_gradientn(name = "difcor",
                              colors = c("aquamarine3","grey90","mediumpurple"),
                              #limits = c(0, 0.28),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "none")+
  guides(fill = guide_legend(ncol = 1))+
  theme_graph(base_family = "sans")+theme(legend.box = 'horizontal',
                                          legend.box.just = 'top')+ggtitle("A.")+
  theme(plot.title = element_text(hjust = -0.02,vjust = 4.5),title = element_text(size = 16,face="bold"))


importance<-read.csv("output/doubleml_res.csv")
importance<-gsub("^m_","",unique(importance$X))
importance<-gsub("_m$","_e",importance)
flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv",header = T)
meta<-read.csv("data/discovery_cohort/16s/SraRunTable.csv",row.names = 1)
selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
meta<-meta[meta$Sample.Name%in%selected$x,]
flux<-flux[flux$sample_id%in%meta$Sample.Name,]
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
flux<-flux[flux$metabolite%in%importance,c(1,2,7)]
list<-sort(table(flux$taxon),decreasing = T)
#list<-list[-which(names(list)=="medium")]


species_cluster<-cluster
species_cluster<-as.data.frame(species_cluster)
species_cluster$spcies<-rownames(species_cluster)
species_cluster$x<-paste("cluster",species_cluster$x,sep = "_")
species_cluster<-split(as.data.table(species_cluster),by=c("x"))
species_cluster<-lapply(species_cluster,function(x){x[[2]]})
es<-fgsea(species_cluster, list,nperm=1000)


myRankedlist<-as.numeric(list)
names(myRankedlist)<-names(list)
cluster_3nodes<-intersect(names(cluster)[cluster==which.min(es$pval)],names(myRankedlist))
write.csv(cluster_3nodes,"output/cluster_3.csv",row.names = F)


meta<-read.table("data/discovery_cohort/16s/sampleda.txt",row.names = 1)
meta<-meta[meta$sample%in%selected$x,]
flux_h<-flux[flux$sample_id%in%meta$sample[meta$group=="healthy"],]
flux_d<-flux[flux$sample_id%in%meta$sample[meta$group=="hypertension"],]
flux_h$sample_id<-1
flux_d$sample_id<-1
aggreate_flux_h<-aggregate(sample_id~taxon+metabolite,data=flux_h,sum)
aggreate_flux_d<-aggregate(sample_id~taxon+metabolite,data=flux_d,sum)
aggreate_flux_h$metabolite<-paste(aggreate_flux_h$metabolite,"healthy",sep = "_")
aggreate_flux_d$metabolite<-paste(aggreate_flux_d$metabolite,"hypertension",sep = "_")
aggreate_flux<-rbind(aggreate_flux_h,aggreate_flux_d)

agg_mat<-as.data.frame(pivot_wider(aggreate_flux,names_from = taxon,values_from = sample_id))
agg_mat[is.na(agg_mat)]<-0
rownames(agg_mat)<-agg_mat$metabolite
agg_mat<-agg_mat[,-1]
agg_mat<-agg_mat[,order(colSums(agg_mat),decreasing = T)]
agg_mat<-agg_mat[order(rowSums(agg_mat)),]
rownames(agg_mat)<-gsub("_e","",rownames(agg_mat))


agg_mat1 <- reshape2::melt(as.matrix(agg_mat)) 
colnames(agg_mat1)[3]<-"Nsample"
custom_palette <- colorRampPalette(c("grey94", "#FEE08B", "#FDAE61", "orangered3"))
set.seed(444)
cluster<-cluster_louvain(dif_graph)
#cluster<-cluster_walktrap(dif_graph)
subg_node<-names(membership(cluster))[membership(cluster)==3]
subg_node<-intersect(subg_node,colnames(agg_mat))

selected_s<-subg_node
label_data <- agg_mat1 %>% 
  filter(Var2 %in% selected_s) %>% 
  distinct(Var2)
gheat<-ggplot(agg_mat1, aes(x = Var2, y = Var1)) +  
  geom_tile(aes(fill = Nsample), colour = "white") +  
  scale_fill_gradientn(name = "Detected\nsample count",colors = custom_palette(100)) +  
  theme_minimal() +  
  labs(y = "Important Features", x = "")+
  theme(axis.text.x = element_blank(),#axis.text.x = element_text(angle = 45,hjust = 1,size = 4),
        axis.title.y = element_text(face = "bold",size = 12),
        axis.text.y = element_text(size=9,face = "bold"),
        legend.key.size = unit(0.3, "cm")) +coord_cartesian(clip = "off") +
  theme(plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))+
  geom_text_repel(
    data = label_data,
    aes(x = Var2, y = 0.5, label = Var2),
    nudge_y = -3,  # 调整标签位置
    direction = "x",
    angle = 45,color="black",vjust=0.1,facet="bold",
    hjust = 0.6,size=3.5,box.padding = 0.1,parse = T,segment.alpha=0.2
  )+
  expand_limits(y = c(-6, 1.5))+
  #scale_y_discrete(breaks = selected_s,labels=selected_s) + 
  theme(plot.margin = ggplot2::margin(t=-2,b=1,r=3),
        title = element_text(size = 12,face="bold") )


aggreate_flux_h<-reshape2::melt(as.matrix(agg_mat[grep("healthy",rownames(agg_mat)),]))
colnames(aggreate_flux_h)<-c("metabolites","taxon","count")
aggreate_flux_h<-aggregate(count~taxon,data=aggreate_flux_h,sum)
aggreate_flux_d<-reshape2::melt(as.matrix(agg_mat[grep("hypertension",rownames(agg_mat)),]))
colnames(aggreate_flux_d)<-c("metabolites","taxon","count")
aggreate_flux_d<-aggregate(count~taxon,data=aggreate_flux_d,sum)

aggreate_flux_h$group<-"healthy"
aggreate_flux_d$group<-"hypertension"
species<-rbind(aggreate_flux_h,aggreate_flux_d)
colnames(species)[2]<-"counts"
species$taxon<-factor(species$taxon,levels = colnames(agg_mat))
bar_species<-ggplot()+
  geom_bar(data = species[species$group=="healthy",],aes(taxon,counts,fill=group),
           stat = "identity", color=NA,alpha=.6)+
  geom_bar(data = species[species$group=="hypertension",],aes(taxon,counts,fill=group),
           stat = "identity", color=NA,alpha=.6)+
  scale_fill_manual(values = c("#5b94c2","#e54445"))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x = element_blank(),
        # axis.text.x = element_text(face = "italic", color="black"),
        #axis.line = element_line(colour = "black"),
        plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
        axis.text.y = element_text(size=10, color="black")) +
  # scale_x_discrete(position = "top")+
  theme(axis.title.x = element_blank(),  
        title = element_text(size = 12,face = "bold"),
        legend.key.size = unit(0.5, "cm"),legend.text = element_text(size=10))+
  #ylim(c(0,900))+
  #geom_text(aes(label = freq),vjust = -0.5, size = 4.5, color = "black")+
  xlab("")+ylab("Sum of samples")

names(species_cluster)[3]<-"FERM"
term2gene <- stack(species_cluster)
colnames(term2gene) <- c("feature", "set")
term2gene<-term2gene[,c(2,1)]
keep<-names(table(term2gene$set))[table(term2gene$set)>=5]
term2gene$set<-as.character(term2gene$set)
term2gene<-term2gene[term2gene$set%in%keep,]
term2gene$set<-factor(term2gene$set,levels = unique(term2gene$set))
gsea_res <- GSEA(myRankedlist, 
                 TERM2GENE = term2gene,
                 minGSSize = 0,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "fdr",
                 seed = 123
)
gsea_df<-as.data.frame(gsea_res@result)
write.table(gsea_df,"output/gsea_df.txt",row.names = F,quote = F,sep = "\t")
source("script/gseaplot3.R")
gseap<-gseaplot3(gsea_res,geneSetID = "FERM",base_size=10,pvalue_table = T)

cluster3nodes<-read.csv("output/cluster_3.csv")
clu3net<-difnetwork[as.character(difnetwork$Var1)%in%cluster3nodes$x&as.character(difnetwork$Var2)%in%cluster3nodes$x,]
clus3no<-spedif_nodes[as.character(spedif_nodes$nodes)%in%as.character(cluster3nodes$x),]
c3_graph<-graph_from_data_frame(d=clu3net,vertices = clus3no)
V(c3_graph)$cluster<-factor(V(c3_graph)$cluster)
c3graph<-ggraph(c3_graph,
                layout = "manual",
                x = pos[names(V(c3_graph)), 1],
                y = pos[names(V(c3_graph)), 2]) +
  geom_edge_link0(aes(col = cor), width = 0.2) +
  geom_node_point(aes(fill = cluster), shape = 21, size = 3) +
  geom_mark_hull( # 注意这里哈
    aes(x, y, group = cluster, fill = cluster),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("red")) +
  scale_fill_manual(values = c("red")) +
  scale_edge_colour_gradientn(name = "difcor",
                              colors = c("grey90","skyblue1","mediumpurple"),
                              limits = c(0, 0.28),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "none")+
  guides(fill = guide_legend(ncol = 1))+
  #theme_graph(base_family = "sans")+theme(legend.position = "node")+
  scale_x_continuous(expand = expansion(add = 1))+
  scale_y_continuous(expand = expansion(add = 1))+
  theme(
    panel.background = element_rect(fill = "transparent"), # 设置绘图区背景为透明
    plot.background = element_rect(fill = "transparent", color = NA), # 设置整个图形背景为透明
    panel.grid.major = element_blank(), # 移除主要网格线
    panel.grid.minor = element_blank(), # 移除次要网格线
    legend.position = "none"
  )



gseap[[1]]<-gseap[[1]]+theme(axis.title.y = element_text(size=12,face = "bold"),
                             axis.text.y = element_text(size=10))

gsea<-wrap_plots(list(gseap[[1]],gseap[[2]],bar_species,gheat),ncol = 1,heights = c(0.38,0.05,0.25,1.2))
gsea<-ggdraw(gsea) +ggtitle("B.")+
  theme(plot.title = element_text(hjust = 0.01,vjust = 0),
        plot.margin = ggplot2::margin(t=15),
        title = element_text(size = 18,face="bold"))+
  draw_plot(c3graph,x = 0.3,y=0.8,width = 0.2,height = 0.12)



known<-read.csv("data/gutmdisorder.csv")
known<-unique(gsub(" .*","",known$Gut.Microbiota..ID.))
subg_node<-names(membership(cluster))[membership(cluster)==3]
subg_node<-intersect(subg_node,colnames(agg_mat))
dif_subgraph<-difnetwork[difnetwork$Var1%in%subg_node&difnetwork$Var2%in%subg_node,]
subg_nodeattr<-nodeatr[nodeatr$nodes%in%subg_node,]
subg_nodeattr$Nsample<-myRankedlist[match(subg_nodeattr$nodes,names(myRankedlist))]
subg_nodeattr$nodetype<-ifelse(as.character(subg_nodeattr$nodes)%in%known,"Documented","Not Documented")
TCG<-read.csv("data/TCG.csv")
TCG$TCG<-rep(c("C1B","C1A"),c(91,50))
TCG<-TCG[TCG$Genus!="",]
TCG<-TCG[!duplicated(TCG$Genus),]
TCG_A<-TCG[TCG$TCG=="C1A",]
TCG_B<-TCG[TCG$TCG=="C1B",]

genus_dif<-read.csv("output/genus_dif.csv",row.names = 1)
genus_dif<-genus_dif[!is.na(genus_dif$ef_lda),]
genus_h<-genus_dif[genus_dif$enrich_group=="healthy",]
genus_d<-genus_dif[genus_dif$enrich_group=="hypertension",]
rownames(genus_h)<-gsub("g__","",rownames(genus_h))
genus_h$sign<-case_when(genus_h$padj<=0.001~"***",
                        genus_h$padj<=0.01&genus_h$padj>0.001~"**",
                        genus_h$padj<=0.05&genus_h$padj>0.01~"*",
                        .default = "")
rownames(genus_d)<-gsub("g__","",rownames(genus_d))
genus_d$sign<-case_when(genus_d$padj<=0.001~"***",
                        genus_d$padj<=0.01&genus_d$padj>0.001~"**",
                        genus_d$padj<=0.05&genus_d$padj>0.01~"*",
                        .default = "")
sign_all<-rbind(genus_d,genus_h)
sign_all$nodes<-paste0(rownames(sign_all),sign_all$sign)
index<-which(subg_nodeattr$nodes%in%rownames(sign_all))
subg_nodeattr$nodes<-as.character(subg_nodeattr$nodes)
subg_nodeattr$nodes[index]<-sign_all$nodes[match(subg_nodeattr$nodes[index],rownames(sign_all))]
subg_nodeattr$change<-ifelse(subg_nodeattr$nodes%in%paste0(rownames(genus_h),genus_h$sign),"Down",ifelse(subg_nodeattr$nodes%in%paste0(rownames(genus_d),genus_d$sign),"Up","Not"))
subg_nodeattr<-subg_nodeattr[order(subg_nodeattr$change),]
subg_nodeattr$change<-factor(subg_nodeattr$change,levels = c("Down","Up","Not"))

dif_subgraph$Var1<-as.character(dif_subgraph$Var1)
dif_subgraph$Var2<-as.character(dif_subgraph$Var2)
index<-which(dif_subgraph$Var1%in%rownames(sign_all))
dif_subgraph$Var1[index]<-sign_all$nodes[match(dif_subgraph$Var1[index],rownames(sign_all))]
index<-which(dif_subgraph$Var2%in%rownames(sign_all))
dif_subgraph$Var2[index]<-sign_all$nodes[match(dif_subgraph$Var2[index],rownames(sign_all))]

subgraph<-graph_from_data_frame(d=dif_subgraph,vertices = subg_nodeattr)
V(subgraph)$nodes<-subg_nodeattr$nodes
V(subgraph)$change<-factor(V(subgraph)$change,levels = c("Down","Up","Not"))

legend_data <- data.frame(
  change = c("Down", "Up"),
  label = c("Down", "Up")
)
subgraph_p<-ggraph(subgraph, layout='linear', circular = TRUE) +
  geom_node_text(aes(x = 1.1 * x,
                     y = 1.1 * y,
                     label = nodes,color=change,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90),
                 hjust = "outward") +
  scale_color_manual(name="change",values=c("Up"="red","Down"="palegreen4","Not"="black"),
                     breaks = c("Up","Down"),
                     labels=c(expression(HT[Up]),expression(HT[Down])))+
  ggnewscale::new_scale_color()+
  geom_node_point(size=3,alpha=1,aes(color=nodetype))+
  geom_edge_arc(aes(colour = cor), width=0.1,
                alpha = 1.5) +
  scale_color_manual(name="node type",values=c("Documented"="red","Not Documented"="darkcyan"))+
  scale_edge_width_continuous(range = c(0,0.2)) +theme_void()+
  #theme_graph()+
  scale_edge_colour_gradientn(name = "Differential \ncorrelation",
                              colors = c("aquamarine3","grey90","mediumpurple"),
                              #limits = c(0, max(E(graph)$cor)),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "edge_colourbar")+
  expand_limits(x = c(-3.5, 2.5), y = c(-3, 3))+theme(plot.margin = ggplot2::margin(r=50))+
  ggtitle("C.")+labs(subtitle = "FERM")+
  theme(plot.title = element_text(hjust = 0.06,vjust = -6),
        plot.subtitle = element_text(hjust = 0.5),
        plot.margin = ggplot2::margin(l=0,r=5,t=5,b=5),
        title = element_text(size = 16,face="bold"))

data_prob<-data.frame(class=c("FERM","other guilds"),n=c(sum(V(subgraph)$change!="Not"),nrow(sign_all)-sum(V(subgraph)$change!="Not")))
data_prob$prop<-data_prob$n/sum(data_prob$n)
data_prob <- data_prob %>%
  arrange(desc(class)) %>% # 重排序
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
data_prob$per<-paste0(round(data_prob$prop,3)*100,"%")
mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")

per_p<-ggplot(data_prob, aes(x = 2, y = prop, fill = class)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = lab.ypos, label = per), color = "black",size=8)+
  scale_fill_manual(values = mycols) +
  theme_void()+labs(subtitle="Proportion of Differential Species")+
  ggtitle("D.")+
  theme(plot.title = element_text(hjust = -0.35,vjust = -2),
        plot.subtitle = element_text(hjust = 0.5,vjust=-5),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16,face="bold"))

spe_net_plot<-gridExtra::grid.arrange(difgraph,subgraph_p,ncol=1,heights=c(0.4,0.5))


ab_c3<-ab_genus
meta<-read.csv("data/discovery_cohort/16s/SraRunTable.csv")
selected<-read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
meta<-meta[meta$Sample.Name%in%selected$x,]
meta$sample_id<-gsub("-","\\.",meta$Sample.Name)
ab_c3<-ab_c3[names(V(dif_graph)),]
MEList <- WGCNA::moduleEigengenes(t(ab_c3), colors = V(dif_graph)$cluster)
MES0<-MEList$eigengenes
colnames(MES0)<-paste0("c",1:9)
meta<-cbind(meta,MES0)
meta$city<-as.factor(meta$city)
meta$sex<-as.factor(meta$sex)
meta$smoker<-as.factor(meta$smoker)
meta$Medicament_use<-as.factor(meta$Medicament_use)
stat_sys<-data.frame()
stat_dia<-data.frame()
for(i in 1:9){
  meta$tmp<-meta[[paste0("c",i)]]
  fit1<-lmerTest::lmer(Systolic_BP~tmp+AGE+Adiponectin+BMI+Body_Fat_Percentage+Glucose+Hb1Ac+
                         HDL+Insulin+Total_Cholesterol+Triglycerides+VLDL+Waist_circumference+
                         (1|city)+(1|sex)+(1|smoker)+(1|Medicament_use),data = meta)
  fit2<-lmerTest::lmer(diastolic_bp~tmp+AGE+Adiponectin+BMI+Body_Fat_Percentage+Glucose+Hb1Ac+
                         HDL+Insulin+Total_Cholesterol+Triglycerides+VLDL+Waist_circumference+
                         (1|city)+(1|sex)+(1|smoker)+(1|Medicament_use),data = meta)
  tmp1<-coef(summary(fit1))[2,]
  tmp2<-coef(summary(fit2))[2,]
  stat_sys<-rbind(stat_sys,tmp1)
  stat_dia<-rbind(stat_dia,tmp2)
}
colnames(stat_dia)<-colnames(stat_sys)<-names(tmp1)
rownames(stat_dia)<-rownames(stat_sys)<-c(paste("guild",1:2,sep="_"),"guild_3(FERM)",paste("guild",4:9,sep="_"))
stat_sys$group<-"Systolic_BP";stat_dia$group<-"Diastolic_BP"
stat_sys$guild<-stat_dia$guild<-rownames(stat_dia)
stat_all<-rbind(stat_dia,stat_sys)
stat_all$p_signif<-case_when(stat_all$`Pr(>|t|)`<0.001~"***",
                             stat_all$`Pr(>|t|)`<0.01~"**",
                             stat_all$`Pr(>|t|)`<0.05~"*",.default = "")
stat_all$`-Log10(Pvalue)`<- -log10(stat_all$`Pr(>|t|)`)
write.csv(stat_all,"output/stat_guild.csv",row.names = F,quote = F)
p_heat<-ggplot(stat_all,aes(guild,group))+
  geom_tile(fill="white",color="black")+
  geom_point(aes(size =`-Log10(Pvalue)`,fill=Estimate),shape=22)+
  geom_text(aes(label=p_signif),size=4,color="white",
            hjust=0.5,vjust=0.7)+
  scale_fill_gradient2(name="coefficient",low = "blue4",mid = "white",high = "red4")+
  scale_size_continuous(range = c(3,8))+
  ggtitle("A.")+
  theme(
        text = element_text(size = 16),axis.title.x = element_text(size = 14,face = "bold"),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold",size = 16),
        plot.margin = margin(t=30,b=0),legend.key.size = unit(4,"mm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 45,hjust=1),
        legend.position = "bottom",
        legend.direction = "vertical" ,
        plot.title.position = "plot",title = element_text(size = 16),
        plot.title = element_text(vjust = 5,hjust = 0.08,face = "bold"),
        )+
  labs(y="",x="")
p_heat
tmp1<-meta[,c("Systolic_BP","c3")]
tmp1$group<-"Systolic_BP"
tmp2<-meta[,c("diastolic_bp","c3")]
tmp2$group<-"Diastolic_BP"
colnames(tmp1)<-colnames(tmp2)<-c("value","c3","group")
lin_data<-rbind(tmp1,tmp2)
sim_text<-data.frame(x=c(0.14,0.14),y=c(170,200),
                     lab=c(paste0("Coefficient=",round(stat_all$Estimate[c(3,12)],2),"\n",
                                  "Pvalue=",round(stat_all$`Pr(>|t|)`[c(3,12)],3))),
                     group=c("Diastolic_BP","Systolic_BP"))
p_line<-ggplot(lin_data,aes(x = c3,y=value))+
  geom_point(aes(color=group))+
  ggpmisc::stat_poly_line(aes(color=group),formula = y ~ x,alpha=0.2)+
  geom_text(data = sim_text,mapping = aes(x=x,y=y,label = lab,color=group),size=5)+
  scale_color_manual(name="",values = c("#acd372","#fbb05b"))+
  theme_classic()+
  labs(x="FERM(eigenvector)",y="BP")+
  ggtitle("B.")+
  theme(
    axis.text = element_text(color = "black",size = 12),
    axis.title = element_text(size = 14,face = "bold"),
    plot.title.position = "plot",
    plot.title = element_text(vjust = -2,hjust = -0.015,face = "bold"),
    strip.text = element_text(face = "bold",size = 16),
    legend.key.size = unit(0.5,'cm'), title = element_text(size=16),
    legend.title = element_text(size = 14,face = "bold"),legend.background = element_blank(),
    plot.background = element_blank(),legend.position = c(0.8,0.15),
    legend.text = element_text(size = 12) 
  )

load("data/sampleMetadata.rda")
hypertension_meta<-sampleMetadata[which(sampleMetadata$disease=='hypertension'),]
healthy_meta<-hypertension_meta[which(hypertension_meta$age_category=='adult'),]
healthy_meta<-sampleMetadata[which(sampleMetadata$disease=='healthy'),]
healthy_meta<-healthy_meta[which(healthy_meta$age_category=='adult'),]
healthy_meta<-healthy_meta[which(healthy_meta$country=='CHN'|healthy_meta$country=='ITA'|healthy_meta$country=='AUT'),]
healthy_meta<-healthy_meta[which(healthy_meta$BMI>18.5&healthy_meta$BMI<23.9),]
healthy_meta<-healthy_meta[which(healthy_meta$age>55&healthy_meta$age<70),]
hypertension_meta<-rbind(hypertension_meta,healthy_meta)
write.csv(hypertension_meta,"output/sum1_table.csv",quote = F,row.names = F)
ab_hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance",rownames = "long")
abundance_hypertension<-ab_hypertension@assays@data@listData[["relative_abundance"]]

ab<-abundance_hypertension|>as.data.frame()
taxa<-strsplit(rownames(ab),"\\|")
genus<-unlist(lapply(taxa, function(x){x[6]}))
genus<-gsub("g__","",genus)
ab<-data.frame(genus=genus,ab)
ab_genus<-aggregate(.~genus,data=ab,sum)
rownames(ab_genus)<-ab_genus$genus
ab_genus<-ab_genus[,-1]
c3<-read.csv("output/cluster_3.csv")[,1]
c3<-intersect(c3,rownames(ab_genus))
pc<-WGCNA::moduleEigengenes(t(ab_genus[c3,]), colors = rep("FERM",16))
pc<-pc$eigengenes
pc$group<-hypertension_meta$disease
wilcox.test(MEFERM~group,data = pc)

p_validation<-ggplot(pc,aes(group,MEFERM,fill = group))+
  geom_boxplot()+
  scale_fill_manual(values = c("blue4","red4"))+
  coord_cartesian(ylim = c(min(pc$MEFERM),0.1))+
  stat_compare_means(comparisons = list(c("hypertension","healthy")),
                     label.y = 0.0,tip.length = 0)+
  theme_classic()+labs(x="",y="FERM(eigenvector)")+ggtitle("C.")+
  theme(
    text = element_text(size = 16),axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),plot.margin = margin(b=30,t=10,r=10,l=5),
    strip.text = element_text(face = "bold",size = 16),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    axis.text.x = element_blank(),legend.position = c(0.8,0.8),
    plot.title.position = "plot",title = element_text(size = 16),
    plot.title = element_text(vjust =-1,hjust = 0,face = "bold"),
  )

pdf("output/Supplementary2.pdf",width = 12,height = 4.5)
p_cor<-grid.arrange(p_heat,p_line,p_validation,ncol=3)
dev.off()

pdf("output/Fig3.pdf",width = 14,height = 10)
grid.arrange(spe_net_plot,gsea,ncol=2,widths=c(0.45,0.55))
dev.off()

pdf("output/Fig3.pdf",width = 14,height = 14)
grid.arrange(tmp,p_cor,nrow=2,heights=c(7,3))
dev.off()




##########Supplementary figures
sup<-list()
for(i in c(1:9)){
  select_cl<-paste0("c",i)
  tmp1<-meta[,c("Systolic_BP",select_cl)]
  tmp1$group<-"Systolic_BP"
  tmp2<-meta[,c("diastolic_bp",select_cl)]
  tmp2$group<-"Diastolic_BP"
  colnames(tmp1)<-colnames(tmp2)<-c("value","c3","group")
  lin_data<-rbind(tmp1,tmp2)
  sim_text<-data.frame(x=c(max(lin_data$c3)-0.25,max(lin_data$c3)-0.25),y=c(150,190),
                       lab=c(paste0("Coefficient=",round(stat_all$Estimate[c(i,i+9)],2),"\n",
                                    "Pvalue=",round(stat_all$`Pr(>|t|)`[c(i,i+9)],3))),
                       group=c("Diastolic_BP","Systolic_BP"))
  p_line<-ggplot(lin_data,aes(x = c3,y=value))+
    geom_point(aes(color=group))+
    ggpmisc::stat_poly_line(aes(color=group),formula = y ~ x,alpha=0.2)+
    geom_text(data = sim_text,mapping = aes(x=x,y=y,label = lab,color=group),size=5)+
    scale_color_manual(name="",values = c("#acd372","#fbb05b"))+
    theme_classic()+
    labs(x=paste0("guild_",i),y="BP")+
    theme(
      axis.text = element_text(color = "black",size = 12),
      axis.title = element_text(size = 14,face = "bold"),
      plot.title.position = "plot",
      plot.title = element_text(vjust = -2,hjust = -0.015,face = "bold"),
      strip.text = element_text(face = "bold",size = 16),
      legend.key.size = unit(0.5,'cm'), title = element_text(size=16),
      legend.title = element_text(size = 14,face = "bold"),
      plot.background = element_blank(),#legend.position = c(0.8,0.35),
      legend.text = element_text(size = 12) 
    )
  sup[[i]]<-p_line
}
sup<-sup[-c(3,12)]
sup[[1]]<-sup[[1]]+ggtitle("A.")
pdf("output/Supplementary2.pdf",width = 15,height = 10)
wrap_plots(sup)
dev.off()
