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
library(ggtext)
library(aplot)
library(patchwork)
library(cowplot)
library(gridExtra)
library(ggforce)
library(limma)
library(ggrepel)
flux<-read.csv("data/discovery_cohort/exchanges.csv")
flux<-flux[,c(2,7,5)]
flux<-flux %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux<-flux%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
meta<-read.csv("data/meta.csv",row.names = 1)
flux<-flux[meta$sample_id,]
flux[is.na(flux)]<-0

keep<-apply(flux, 2, function(x){
  mean(x!=0)>0.2
})
flux<-flux[,keep]
net <- netConstruct(data = flux[meta$disease=="hypertension",],data2 = flux[meta$disease=="healthy",], 
                    filtTax = "none",
                    filtTaxPar = list(relFreq = 0),
                    filtSamp = "highestFreq",
                    filtSampPar = list(highestFreq = sum(meta$disease=="healthy")),
                    measure = "pearson",
                    # measurePar = list(nlambda=10, 
                    #                   rep.num=10,
                    #                   Rmethod = "approx"),
                    normMethod = "none", 
                    zeroMethod = "none",
                    sparsMethod = "threshold",
                    thresh = 0.5,
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
                  showCentr="all",showJacc=F,
                  numbNodes = 2000)
diff <- diffnet(net,discordThresh=0.8,alpha = 0.05,lfdrThresh = 0.05,
                diffMethod = "fisherTest", nPerm = 10,cores = 12,
                adjust = "lfdr")

props_pears <- netAnalyze(net, 
                          clustMethod = "cluster_fast_greedy",
                          weightDeg = TRUE,
                          normDeg = FALSE,
                          gcmHeat = FALSE)


top=10
difdegree=summ_res$topProps$topDeg
difBet=summ_res$topProps$topBetw
difclose=summ_res$topProps$topClose
difeig=summ_res$topProps$topEigen

attr<-list(
  difdegree=rownames(summ_res$topProps$topDeg)[order(difdegree$abs.diff.,decreasing = T)[1:top]],
  difBet=rownames(summ_res$topProps$topBetw)[order(difBet$abs.diff.,decreasing = T)[1:top]],
  difclose=rownames(summ_res$topProps$topClose)[order(difclose$abs.diff.,decreasing = T)[1:top]],
  difeig=rownames(summ_res$topProps$topEigen)[order(difeig$abs.diff.,decreasing = T)[1:top]]
)
top<-unique(unlist(attr))
attr<-list_to_data_frame(attr)
ggvenn(data = attr,show_percentage = F)+
  scale_fill_brewer(palette = "Set3")

deg<-difdegree[top,1:2]
deg<-deg[order(deg[,1]-deg[,2]),]
deg<-reshape2::melt(as.matrix(deg))
deg$group<-"Degree"

bet<-summ_res$topProps$topBetw[top,1:2]
bet<-bet[order(bet[,1]-bet[,2]),]
bet<-reshape2::melt(as.matrix(bet))
bet$group<-"Betweeness"

clo<-summ_res$topProps$topClose[top,1:2]
clo<-clo[order(clo[,1]-clo[,2]),]
clo<-reshape2::melt(as.matrix(clo))
clo$group<-"Closeness"

eig<-summ_res$topProps$topEigen[top,1:2]
eig<-eig[order(eig[,1]-eig[,2]),]
eig<-reshape2::melt(as.matrix(eig))
eig$group<-"Eigen"
attr_all<-rbind(deg,bet,clo,eig)
attr_all$Var2<-factor(attr_all$Var2,levels = c("healthy","hypertension"))


colnames(attr_all)[2]<-"Group"
barp<-ggplot() +
  geom_bar(data=attr_all[attr_all$Group=="healthy",], mapping = aes(y = Var1, x = log10(value+1),fill=Group),stat = "identity", alpha = 0.7) +
  geom_bar(data=attr_all[attr_all$Group=="hypertension",],mapping = aes(y = Var1, x = log10(value+1),fill=Group),stat = "identity", alpha = 0.7)+
  facet_grid(.~group, scales = "free") +
  #geom_jitter(aes(color = group),alpha = 0.5,
  #            position = position_jitterdodge(dodge.width = 1))+
  scale_fill_manual(name="group",values = c("hypertension"="#ea5b57", "healthy"="#5b94c2")) +
  #scale_color_manual(values = c("hypertension"="#ea5b57", "healthy"="#5b94c2")) +
  theme_bw() +
  labs(x="",y="")+
  theme(axis.text.x = element_text(size = 8,face = "bold"),
        axis.text.y = element_text(size = 12,face="bold"),
        legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        plot.margin = ggplot2::margin(b=0,t=10),
        strip.text = element_text(face = "bold",size = 12))

barp<-ggdraw(barp)+ggtitle("B.")+
  theme(plot.title = element_text(hjust = 0.01,vjust = -2),
        plot.title.position = "plot",
        title = element_text(size = 16,face="bold"))



importanceflux<-read.csv("data/xgb_rf_importance.csv")
importanceflux<-gsub("^m_","",importanceflux[,1])
difnetwork<-diff$assoMat1-diff$assoMat2
difnetwork[diff$pAdjustMat>=1e-6]<-0
rownames(difnetwork)<-gsub("\\[e\\]$","_e",rownames(difnetwork))
colnames(difnetwork)<-gsub("\\[e\\]$","_e",colnames(difnetwork))
difnetwork_imp<-difnetwork[importanceflux,]
difnetwork_imp<-difnetwork_imp[rowSums(difnetwork_imp)!=0,colSums(difnetwork_imp)!=0]
imp_nodes<-unique(c(rownames(difnetwork_imp),colnames(difnetwork_imp)))
difnode<-as.character(unique(attr_all$Var1[attr_all$Var1%in%imp_nodes]))
showlabel<-c(importanceflux)
imp_nodes<-data.frame(nodes=imp_nodes,nodetype=rep(c("Important metabolite","Neighbors"),c(9,length(imp_nodes)-9)))
#imp_nodes$nodetype[which(imp_nodes$nodes%in%difnode)]<-"Top10"
difnetwork_imp<-reshape2::melt(as.matrix(difnetwork_imp))
difnetwork_imp<-difnetwork_imp[difnetwork_imp$value!=0,]
#difnetwork_imp<-difnetwork_imp[abs(difnetwork_imp$value)>0.5,]
flux_difgraph<-graph_from_data_frame(difnetwork_imp,vertices = imp_nodes)
V(flux_difgraph)$nodes<-imp_nodes$nodes


fluxdifgraph<-ggraph(flux_difgraph, layout = "stress") + #layout_tbl_graph_igraph unrooted stress
  geom_edge_link(aes(colour = value,
                     # width=difcor
  ),
  #  linetype = type), 
  alpha = 1) +
  
  scale_edge_colour_gradientn(name = "Differential \ncorrelation",
                              colors = c("aquamarine4","grey93","mediumpurple"),
                              #limits = c(0, 0.7),
                              space = "Lab",
                              na.value = "grey50", 
                              guide = "edge_colourbar") +
  
  geom_node_point(aes(
    #size =  30,
    #shape = nodetype,
    colour = nodetype,
    #fill = difclose,
    stroke = 1.15),size=2,alpha=0.8) +
  
  scale_colour_manual(name = "node type",
                      values = c("Neighbors"="#0f8096","Important metabolite"="#f09d69")
  ) + #coord_fixed() +  # 增加图形范围
  geom_node_label(aes(filter = nodes%in%showlabel,
                     label = paste(nodes, sep = '')),max.overlaps = Inf,fill = scales::alpha("white", 0.6),
                 colour="black", repel = T,size=5) +
  theme_void() +ggtitle("C.")+
  guides(color = guide_legend(override.aes = list(size = 6)))+
  theme(plot.title = element_text(size = 16, face = "bold",vjust = 0,hjust = 0.01),
        #legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        plot.title.position = "plot",
        plot.margin = ggplot2::margin(t=0))
fluxdifgraph
pathway<-read.delim("data/ReactionDatabase.txt")
pathway<-pathway[pathway$Subsystem!="",]
subsys<-unique(pathway$Subsystem)
pathwaylist<-data.frame()
for(i in 1:length(subsys)){
  tmp<-pathway[pathway$Subsystem==subsys[i],]
  metabolite<-unlist(lapply(strsplit(tmp$Formula," -> | <=> | <=>|->"),function(x) x[1:2]))
  metabolite<-unlist(strsplit(metabolite," \\+ "))
  metabolite<-gsub("\\s*$","",unique(gsub("^[0-9.]+\\s*","",metabolite)))
  metabolite<-gsub("\\[.\\]","",metabolite)
  path_tmp<-data.frame(pathway=subsys[i],metabolite=unique(metabolite))
  pathwaylist<-rbind(pathwaylist,path_tmp)
}
implist<-gsub("_e$|_m$","",imp_nodes$nodes)
implist<-unique(gsub("^[0-9.]+","",implist))

flux<-read.csv("data/discovery_cohort/exchanges.csv")
background<-unique(gsub("^[0-9]+|_m$|\\[e\\]","",unique(flux$metabolite)))
source("script/enrichment_analysis.R")
node_degree<-summ_res[["topProps"]][["topDeg"]]
node_degree$diff<-node_degree$hypertension-node_degree$healthy
node_degree<-node_degree[-grep("_m$",rownames(node_degree)),]
rownames(node_degree)<-gsub("\\[e\\]","",rownames(node_degree))
up<-unique(gsub("\\[e\\]$|_m$","",rownames(node_degree)[node_degree$diff>0]))
down<-unique(gsub("\\[e\\]$|_m$","",rownames(node_degree)[node_degree$diff<0]))
imp_overpresentation_up<-enrichment_analysis(intersect(up,implist),background = background,pathway_data = pathwaylist)
imp_overpresentation_up<-imp_overpresentation_up[imp_overpresentation_up$P_Value<0.05,]
imp_overpresentation_up$change<-"UP"
imp_overpresentation_down<-enrichment_analysis(intersect(down,implist),background = background,pathway_data = pathwaylist)
imp_overpresentation_down<-imp_overpresentation_down[imp_overpresentation_down$P_Value<0.05,]
imp_overpresentation_down$change<-"DOWN"
imp_overpresentation_down1<-imp_overpresentation_down[-grep("Transport|Others",imp_overpresentation_down$pathway),]
imp_overpresentation_down<-imp_overpresentation_down1[1:10,]
#imp_overpresentation<-imp_overpresentation[imp_overpresentation$pathway!="Others",]
imp_overpresentation_top<-rbind(imp_overpresentation_up,imp_overpresentation_down)
imp_overpresentation_top$group<-"Important Metabolite Differential Subnetwork"

difdegree=summ_res$topProps$topDeg
difBet=summ_res$topProps$topBetw
difclose=summ_res$topProps$topClose
difeig=summ_res$topProps$topEigen

difdeg<-rownames(difdegree)[order(difdegree$abs.diff.,decreasing = T)][1:100]
difbet<-rownames(difBet)[order(difBet$abs.diff.,decreasing = T)][1:100]
difclo<-rownames(difclose)[order(difclose$abs.diff.,decreasing = T)][1:100]
difeig<-rownames(difeig)[order(difeig$abs.diff.,decreasing = T)][1:100]
dif_atr<-c(difdeg)#,difbet,difclo,difeig)
dif_atr<-unique(gsub("\\[e\\]|_m$","",dif_atr))
overpresentation_up<-enrichment_analysis(intersect(up,dif_atr),background = background,pathway_data = pathwaylist)
overpresentation_up<-overpresentation_up[overpresentation_up$P_Value<0.05,]
overpresentation_up$change<-"UP"
overpresentation_up$group<-"Top 100 Nodes with the Largest Degree Differences in the GMFiN"

overpresentation_down<-enrichment_analysis(intersect(down,dif_atr),background = background,pathway_data = pathwaylist)
overpresentation_down<-overpresentation_down[overpresentation_down$P_Value<0.05,]
overpresentation_down$change<-"DOWN"
overpresentation_down$group<-"Top 100 Nodes with the Largest Degree Differences in the GMFiN"
overpresentation_down<-overpresentation_down[-grep("Transport|Others",overpresentation_down$pathway),]
overpresentation_down<-overpresentation_down
overpresentation<-rbind(overpresentation_up,overpresentation_down)



all_overpres<-rbind(imp_overpresentation_top,overpresentation)
all_overpres$group<-factor(all_overpres$group,levels = c("Top 100 Nodes with the Largest Degree Differences in the GMFiN","Important Metabolite Differential Subnetwork"))
all_overpres$pvalue1<- -log10(all_overpres$P_Value)
all_overpres$pvalue1<-ifelse(all_overpres$change=="DOWN",-log10(all_overpres$pvalue1),log10(all_overpres$pvalue1))

flux_overpres<-all_overpres[all_overpres$group!="Top 100 Nodes with the Largest In-degree Differences in the GMxN",]
flux_enrichment<-ggplot(data = flux_overpres,aes(reorder(pathway, pvalue1), y = pvalue1,fill=change))+
  geom_bar(stat = "identity",position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(labels = function(y) ifelse(y == 0, 0, ifelse(y<0,10^(-y),10^(y))),breaks = c(-log10(10),-log10(5),0,log10(5),log10(10),log10(20),log(40)))+
  #facet_wrap(~group,ncol = 1,scales = "free",space="free") +
  facet_col(vars(group), scales = "free_y", space = "free")+
  coord_flip() +
  labs(title = "",
       x = "",
       y = "-log10(P_value)") +theme(legend.key.size = unit(0.5,"cm"))+
  scale_fill_manual(values = c(UP="firebrick",DOWN="dodgerblue4"))+
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
  geom_text(data = flux_overpres[flux_overpres$change == "UP",], aes(x = pathway, y = -0.01, label = pathway),
            hjust = 1, size = 4) +
  geom_text(data = flux_overpres[flux_overpres$change == "DOWN",], aes(x = pathway, y = 0.01, label = pathway),
            hjust = 0, size = 4)+#expand_limits(y=c(-2,3))+
  ggtitle("D.")+theme(
    plot.title = element_text(size = 16, face = "bold",vjust = -1,hjust = 0.01),
    plot.title.position = "plot",
  )


net_p<-grid::rasterGrob(png::readPNG("output/net_graph1.png"), interpolate = TRUE)
net_p<-ggdraw(net_p) +ggtitle("A.")+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = -0.02,vjust = 1),
        plot.margin = ggplot2::margin(r=30,l=20,b=2,t=10),
        title = element_text(size = 18,face="bold"))
tmp1<-grid.arrange(net_p,barp,ncol=1,heights=c(0.4,0.65))
tmp2<-grid.arrange(fluxdifgraph,flux_enrichment,ncol=1,heights=c(0.45,0.55))
dev.off()
pdf("output/Fig4.pdf",width = 15,height = 12)
grid.arrange(tmp1,tmp2,ncol=2,widths=c(0.55,0.45))
dev.off()



