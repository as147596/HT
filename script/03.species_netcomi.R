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
showtext_auto()
ab<-read.csv("data/abundance.csv",row.names = 1,header = T)
samplename<-as.character(read.csv("data/abundance.csv",row.names = 1,header = F)[1,])
colnames(ab)<-samplename
taxa<-strsplit(rownames(ab),"\\|")
genus<-unlist(lapply(taxa, function(x){x[6]}))
genus<-gsub("g__","",genus)
ab<-data.frame(genus=genus,ab)
ab_genus<-aggregate(.~genus,data=ab,sum)
rownames(ab_genus)<-ab_genus$genus
ab_genus<-ab_genus[,-1]
meta<-read.csv("data/meta.csv")
meta$sample_id<-gsub("-","\\.",meta$sample_id)
ab_genus<-ab_genus[,meta$sample_id]

keep<-apply(ab_genus, 1, function(x){
  mean(x!=0)>0.05
})
ab_genus<-ab_genus[keep,]
net <- netConstruct(data = t(ab_genus)[meta$disease=="hypertension",],data2 = t(ab_genus)[meta$disease=="healthy",], 
                    filtTax = "none",
                    #filtTaxPar = list(relFreq = 0),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta$disease=="healthy")),
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
net <- netConstruct(data = t(ab_genus)[meta$disease=="hypertension",],data2 = t(ab_genus)[meta$disease=="healthy",], 
                    filtTax = "none",
                    #filtTaxPar = list(relFreq = 0.1),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta$disease=="healthy")),
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
                diffMethod = "permute", nPerm = 1000,cores = 15,
                storeCountsPerm = TRUE, 
                fileStoreCountsPerm = c("countsPerm1", "countsPerm2"),
                storeAssoPerm = TRUE,
                fileStoreAssoPerm = "assoPerm",
                adjust = "none",seed = 123456)
saveRDS(diff,"data/specie_network_res.rds")
diff<-readRDS("data/specie_network_res.rds")
plot(diff,
     cexNodes = 1, 
     cexLegend = 0.6,
     cexTitle = 2,
     mar = c(10,10,10,20),layout=NULL,
     legendGroupnames = c("hypertension", "healthy"),
     legendPos = c(1.1,0.8))
props_pears <- netAnalyze(net, 
                          clustMethod = "cluster_fast_greedy",
                          weightDeg = TRUE,
                          normDeg = FALSE,
                          gcmHeat = FALSE)


diffmat_sums <- rowSums(diff$diffAdjustMat)
diff_asso_names <- names(diffmat_sums[diffmat_sums > 0])

plot(props_pears, 
     nodeFilter = "names",
     nodeFilterPar = diff_asso_names,
     nodeColor = "gray",
     highlightHubs = FALSE,
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = FALSE, 
     nodeSize = "clr",
     edgeTranspHigh = 20,
     labelScale = FALSE,
     cexNodes = 1, 
     cexLabels = 0.8,
     cexTitle = 0.8,
     groupNames = c("No seasonal allergies", "Seasonal allergies"),
     hubBorderCol  = "gray40")






top=20
attr<-list(
  difdegree=rownames(summ_res$topProps$topDeg)[1:top],
  difBet=rownames(summ_res$topProps$topBetw)[1:top],
  difclose=rownames(summ_res$topProps$topClose)[1:top],
  difeig=rownames(summ_res$topProps$topEigen)[1:top]
)
top<-unique(unlist(attr))
topDeg<-summ_res$topProps$topDeg
colnames(topDeg)<-paste0("Degree_",colnames(topDeg))
topBetw<-summ_res$topProps$topBetw
colnames(topBetw)<-paste0("Betweenness_",colnames(topBetw))
topClose<-summ_res$topProps$topClose
colnames(topClose)<-paste0("Closeness_",colnames(topClose))
topEigen<-summ_res$topProps$topEigen
colnames(topEigen)<-paste0("Eigen_",colnames(topEigen))
summ_dif<-merge(topDeg,topBetw,by="row.names")
summ_dif<-merge(summ_dif,topClose,by.x="Row.names",by.y="row.names")
summ_dif<-merge(summ_dif,topEigen,by.x="Row.names",by.y="row.names")
summ_dif<-summ_dif[summ_dif$Row.names%in%top,]

degree<-summ_res[["topProps"]][["topDeg"]][top,1:2]
degree<-degree[order(degree[,1]-degree[,2]),]
degree<-reshape2::melt(as.matrix(degree))

bet<-summ_res[["topProps"]]$topBetw[top,1:2]
bet<-bet[order(bet[,1]-bet[,2]),]
bet<-reshape2::melt(as.matrix(bet))

clo<-summ_res[["topProps"]]$topClose[top,1:2]
clo<-clo[order(clo[,1]-clo[,2]),]
clo<-reshape2::melt(as.matrix(clo))

eig<-summ_res[["topProps"]]$topEigen[top,1:2]
eig<-eig[order(eig[,1]-eig[,2]),]
eig<-reshape2::melt(as.matrix(eig))

degree$group<-"Degree"
bet$group<-"Betweeness"
clo$group<-"Closeness"
eig$group<-"Eigenvector"
heat_data<-rbind(degree,bet,clo,eig)

topnodes<-rownames(summ_res[["topProps"]][["topDeg"]][top,1:2])
difnet<-(diff$assoMat1-diff$assoMat2)[topnodes,]

diffmat_sums <- colSums(difnet)
#diff_asso_names<-names(diffmat_sums[diffmat_sums > 0])
diff_asso_names <- unique(c(names(diffmat_sums[diffmat_sums > 0]),topnodes))


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


nodeatr<-data.frame(nodes = nodes,degree=summ_res[["topProps"]][["topDeg"]][nodes,3],
                    Between=summ_res[["topProps"]][["topBetw"]][nodes,3],nodetype="difnetwork")
nodeatr<-nodeatr[nodeatr$nodes%in%unique(c(difnetwork$Var1,difnetwork$Var2)),]
dif_graph<-graph_from_data_frame(d=difnetwork,vertices = nodeatr)
importance<-read.csv("data/xgb_rf_importance.csv")
importance<-gsub("^m_","",importance$x)
flux<-read.csv("data/discovery_cohort/exchanges.csv",header = T)
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
flux<-flux[flux$metabolite%in%importance,c(1,2,7)]
imp_species<-unique(flux$taxon)
spedif_nodes<-nodeatr
spedif_nodes$nodetype="Related"
spedif_nodes$nodetype[spedif_nodes$nodes%in%imp_species]="Importance"
dif_graph<-graph_from_data_frame(d=difnetwork,vertices = spedif_nodes)
cluster<-cluster_walktrap(dif_graph)
cluster<-membership(cluster)

spedif_nodes$cluster<-1
spedif_nodes$cluster[match(names(cluster),spedif_nodes$nodes)]<-paste("cluster",(unname(cluster)),sep = "_")
spedif_nodes$cluster<-factor(spedif_nodes$cluster,levels = paste("cluster",1:13,sep = "_"))
spedif_nodes<-spedif_nodes[order(spedif_nodes$cluster),]
dif_graph<-graph_from_data_frame(d=difnetwork,vertices = spedif_nodes)
tmp<-V(dif_graph)$cluster
tmp[tmp=="cluster_3"]<-"cluster_3(KEPR)"
V(dif_graph)$cluster<-factor(tmp,levels = unique(tmp))
g <- sample_islands(length(table(cluster)), max(table(cluster)), 0.5, 15)
g <- igraph::simplify(g)
set.seed(123)
bb <- layout_as_backbone(g, keep = 0.2)
pos_tmp<-bb$xy
pos<-data.frame()
nodelist<-spedif_nodes$cluster%>%unique()
for(i in 1:13){
  set.seed(123)
  x.pos<-sample((-7.5+i):(-6.5+i))
  tmp<-pos_tmp[((31*(i-1))+1):(31*i),]
  tmp<-tmp[sample(1:31,sum(spedif_nodes$cluster==nodelist[i])),,drop=F]
  pos<-rbind(pos,tmp)
}
difgraph<-ggraph(dif_graph,
       layout = "manual",
       x = pos[, 1],
       y = pos[, 2]) +
  geom_edge_link0(aes(col = cor), width = 0.2) +
  geom_node_point(aes(fill = cluster), shape = 21, size = 3) +
  geom_mark_hull( # 注意这里哈
    aes(x, y, group = cluster, fill = cluster),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("cluster_1"="orangered4","cluster_2"="orchid4","cluster_3(KEPR)"="red","cluster_4"="mediumpurple",
                                "cluster_5"="lightsteelblue4","cluster_6"="linen","cluster_7"="orange3","cluster_8"="royalblue4",
                                "cluster_9"="skyblue4","cluster_10"="thistle","cluster_11"="turquoise4","cluster_12"="palegreen4","cluster_13"="mediumseagreen"),
                     labels=paste("guild",1:13,sep="_")) +
  scale_fill_manual(values = c("cluster_1"="orangered4","cluster_2"="orchid4","cluster_3(KEPR)"="red","cluster_4"="mediumpurple",
                               "cluster_5"="lightsteelblue4","cluster_6"="linen","cluster_7"="orange3","cluster_8"="royalblue4",
                               "cluster_9"="skyblue4","cluster_10"="thistle","cluster_11"="turquoise4","cluster_12"="palegreen4","cluster_13"="mediumseagreen"),
                    labels=c(paste("guild",1:2,sep="_"),"guild_3(KEPR)",paste("guild",4:13,sep="_"))) +
  scale_edge_colour_gradientn(name = "difcor",
                              colors = c("aquamarine3","grey90","mediumpurple"),
                              #limits = c(0, 0.28),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "none")+
  guides(fill = guide_legend(ncol = 1))+
  theme_graph(base_family = "sans")+theme(legend.box = 'horizontal',
                                          legend.box.just = 'top')+ggtitle("A.")+
  theme(plot.title = element_text(hjust = -0.065,vjust = 4.5),title = element_text(size = 16,face="bold"))


importance<-read.csv("data/xgb_rf_importance.csv")
importance<-gsub("^m_","",importance$x)
importance<-gsub("_m$","_e",importance)
flux<-read.csv("data/discovery_cohort/exchanges.csv",header = T)
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
flux<-flux[flux$metabolite%in%importance,c(1,2,7)]
list<-sort(table(flux$taxon),decreasing = T)
#list<-list[-which(names(list)=="medium")]



cluster<-cluster_walktrap(dif_graph)  # 
#cluster<-cluster_leading_eigen(dif_graph)
table(cluster$membership)
plot(cluster,dif_graph,vertex.label=NA,
     edge.arrow.mode='-',    
     vertex.size = 5 )
species_cluster<-membership(cluster)
species_cluster<-as.data.frame(species_cluster)
species_cluster$spcies<-rownames(species_cluster)
species_cluster$x<-paste("cluster",species_cluster$x,sep = "_")
species_cluster<-split(as.data.table(species_cluster),by=c("x"))
species_cluster<-lapply(species_cluster,function(x){x[[2]]})
es<-fgsea(species_cluster, list,nperm=1000)
esplot<-plotEnrichment(species_cluster[[3]],list) + 
  labs(title=head(es[order(pval), ], 1)$pathway)

myRankedlist<-as.numeric(list)
names(myRankedlist)<-names(list)

cluster_3nodes<-names(cluster)[cluster==3]
write.csv(cluster_3nodes,"output/cluster_3.csv")

known<-read.csv("data/known.csv")
known<-unique(gsub(" .*","",known$Gut.Microbiota..ID.))
subg_node<-names(membership(cluster))[membership(cluster)==3]
dif_subgraph<-difnetwork[difnetwork$Var1%in%subg_node&difnetwork$Var2%in%subg_node,]
subg_nodeattr<-nodeatr[nodeatr$nodes%in%subg_node,]
subg_nodeattr$Nsample<-myRankedlist[match(subg_nodeattr$nodes,names(myRankedlist))]
subg_nodeattr$nodetype<-ifelse(subg_nodeattr$nodes%in%known,"Documented","Not Documented")
TCG<-read.csv("data/TCG.csv")
TCG$TCG<-rep(c("C1B","C1A"),c(91,50))
TCG<-TCG[TCG$Genus!="",]
TCG<-TCG[!duplicated(TCG$Genus),]
TCG_A<-TCG[TCG$TCG=="C1A",]
TCG_B<-TCG[TCG$TCG=="C1B",]
subg_nodeattr$TCG<-ifelse(subg_nodeattr$nodes%in%TCG$Genus,ifelse(subg_nodeattr$nodes%in%TCG_A$Genus,"C1A","C1B"),"Other")
subg_nodeattr<-subg_nodeattr[order(subg_nodeattr$TCG),]
subgraph<-graph_from_data_frame(d=dif_subgraph,vertices = subg_nodeattr)
V(subgraph)$nodes<-subg_nodeattr$nodes
subgraph_plot<-ggraph(subgraph, layout = "linear",circular = TRUE) +
  geom_edge_link(aes(colour = cor,
                     # width=difcor
  ),
  #  linetype = type), 
  alpha = 1) +
  
  scale_edge_colour_gradientn(name = "cor",
                              colors = c("#66C2A5", "#ABDDA4", "#E6F598", "#FEE08B", "#FDAE61", "#F46D43"),
                              #limits = c(0, max(E(graph)$cor)),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "edge_colourbar") +
  
  geom_node_point(aes(size = degree,
                      #size =  30,
                      #shape = nodetype,
                      colour = Between,
                      #fill = difclose,
                      stroke = 1.15)) +
  
  
  scale_colour_gradientn(name = "Betweenness",
                         colours = c("skyblue1","slateblue","red"),
                         # limits = c(min(V(tp_graph)$logpval), max(V(tp_graph)$logpval)),
                         #limits = c(min(V(graph)$Between), max(V(graph)$Between))
                         # oob = squish
  ) + coord_fixed() +
  #facet_nodes(~nodetype,nrow = 1)+
  geom_node_text(aes(filter = nodes%in%subg_node,
                     label = paste(nodes, sep = '')), 
                 colour="firebrick4", repel = T) +
  
  theme_void() +
  theme(plot.title = element_text(size = 5, face = "bold"))+
  theme(panel.spacing = unit(1,"lines"),
        #legend.position = 'bottom',
        legend.spacing.x = unit(2,'cm'),
        legend.key.size = unit(0.5,'cm'))


subgraph_p<-ggraph(subgraph, layout='linear', circular = TRUE) +
  geom_node_text(aes(x = 1.1 * x,
                     y = 1.1 * y,
                     label = nodes,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90),
                 hjust = "outward") +
  geom_node_point(size=3,alpha=1,aes(color=nodetype,shape=TCG))+
  geom_edge_arc(aes(colour = cor), width=0.1,
                alpha = 1.5) +
  scale_color_manual(name="node type",values=c("Documented"="red","Not Documented"="palegreen4"))+
  scale_edge_width_continuous(range = c(0,0.2)) +theme_void()+
  #theme_graph()+
  scale_edge_colour_gradientn(name = "Differential \ncorrelation",
                              colors = c("aquamarine3","grey90","mediumpurple"),
                              #limits = c(0, max(E(graph)$cor)),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "edge_colourbar")+
  expand_limits(x = c(-3.5, 2.5), y = c(-3, 3))+theme(plot.margin = ggplot2::margin(r=50))+
  ggtitle("C.")+labs(subtitle = "KEPR")+
  theme(plot.title = element_text(hjust = 0.07,vjust = -6),
        plot.subtitle = element_text(hjust = 0.5),
        title = element_text(size = 16,face="bold"))

spe_net_plot<-gridExtra::grid.arrange(difgraph,subgraph_p,ncol=1,heights=c(0.4,0.5))

flux<-read.csv("data/discovery_cohort/exchanges.csv",header = T)
flux$metabolite<-gsub("\\[e\\]","_e",flux$metabolite)
flux<-flux[flux$metabolite%in%importance,]
flux<-flux[flux$taxon!="medium",c(1,2,7)]
meta<-read.csv("data/meta.csv",row.names = 1)
flux_h<-flux[flux$sample_id%in%meta$sample_id[meta$disease=="healthy"],]
flux_d<-flux[flux$sample_id%in%meta$sample_id[meta$disease=="hypertension"],]
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
rownames(agg_mat)<-gsub("_e$","",rownames(agg_mat))

agg_mat1 <- reshape2::melt(as.matrix(agg_mat)) 
colnames(agg_mat1)[3]<-"Nsample"
custom_palette <- colorRampPalette(c("grey94", "#FEE08B", "#FDAE61", "orangered3"))
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
        axis.text.y = element_text(size=8),
        legend.key.size = unit(0.3, "cm")) +coord_cartesian(clip = "off") +
  theme(plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"))+
  geom_text_repel(
    data = label_data,
    aes(x = Var2, y = 0.5, label = Var2),
    nudge_y = -3,  # 调整标签位置
    direction = "x",
    angle = 45,color="black",vjust=0.1,
    hjust = 0.6,size=2.2,box.padding = 0.1,parse = T,segment.alpha=0.2
  )+
  expand_limits(y = c(-6, 1.5))+
  #scale_y_discrete(breaks = selected_s,labels=selected_s) + 
  theme(plot.margin = ggplot2::margin(t=-2,b=1,r=3),
        title = element_text(size = 12,face="bold") )


aggreate_flux_h<-aggregate(sample_id~taxon,data=aggreate_flux_h,sum)
aggreate_flux_d<-aggregate(sample_id~taxon,data=aggreate_flux_d,sum)
aggreate_flux_h$group<-"healthy"
aggreate_flux_d$group<-"hypertension"
species<-rbind(aggreate_flux_h,aggreate_flux_d)
colnames(species)[2]<-"counts"
species$taxon<-factor(species$taxon,levels = colnames(agg_mat))
bar_species<-ggplot(data = species,aes(taxon,counts,fill=group))+
  geom_bar(stat = "identity", color=NA,alpha=.6)+
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

names(species_cluster)[3]<-"KEPR"
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
                 #pAdjustMethod = "BH",
                 seed = 123
)
source("script/gseaplot3.R")
gseap<-gseaplot3(gsea_res,geneSetID = "KEPR",base_size=10,pvalue_table = T)

cluster3nodes<-read.csv("output/cluster_3.csv")
clu3net<-difnetwork[as.character(difnetwork$Var1)%in%cluster3nodes$x&as.character(difnetwork$Var2)%in%cluster3nodes$x,]
clus3no<-spedif_nodes[spedif_nodes$nodes%in%cluster3nodes$x,]
c3_graph<-graph_from_data_frame(d=clu3net,vertices = clus3no)
V(c3_graph)$cluster<-factor(V(c3_graph)$cluster)
c3graph<-ggraph(c3_graph,
                 layout = "manual",
                 x = pos[spedif_nodes$nodes%in%cluster3nodes$x, 1],
                 y = pos[spedif_nodes$nodes%in%cluster3nodes$x, 2]) +
  geom_edge_link0(aes(col = cor), width = 0.2) +
  geom_node_point(aes(fill = cluster), shape = 21, size = 3) +
  geom_mark_hull( # 注意这里哈
    aes(x, y, group = cluster, fill = cluster),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("cluster_3"="red")) +
  scale_fill_manual(values = c("cluster_3"="red")) +
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

gsea<-wrap_plots(list(gseap[[1]],gseap[[2]],bar_species,gheat),ncol = 1,heights = c(0.5,0.1,0.7,1.2))
gsea<-ggdraw(gsea) +ggtitle("B.")+
  theme(plot.title = element_text(hjust = 0.01,vjust = 0),
        plot.margin = ggplot2::margin(t=15),
        title = element_text(size = 18,face="bold"))+
  draw_plot(c3graph,x = 0.35,y=0.85,width = 0.14,height = 0.11)
pdf("output/Fig3.pdf",width = 14,height = 10)
grid.arrange(spe_net_plot,gsea,ncol=2,widths=c(0.5,0.55))
dev.off()

