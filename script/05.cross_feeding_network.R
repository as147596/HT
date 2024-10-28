library(ggplot2)
library(tidyverse)
library(tidyr)
library(reshape2)
library(igraph)
library(limma)
library(dplyr)
library(ggpubr)

flux<-read.csv("data/discovery_cohort/exchanges.csv")
flux<-flux[,c(2,1,7,5)]
colnames(flux)[2:3]<-c("source","target")
index<-which(flux$flux<0)
flux[index,2:3]<-flux[index,3:2]
meta<-read.csv("data/meta.csv",row.names = 1)
flux<-flux[flux$sample_id%in%meta$sample_id,]
flux$flux<-abs(flux$flux)
res<-list()
for(i in meta$sample_id){
  sample_net<-flux[flux$sample_id==i,]
  sample_net_g<-graph_from_data_frame(sample_net[,2:4],directed = T)
  eigen<- eigen_centrality(sample_net_g,directed=T)$vector
  bet<-betweenness(sample_net_g)
  deg_in<-igraph::degree(sample_net_g,mode="in")
  deg_out<-igraph::degree(sample_net_g,mode="out")
  close<-closeness(sample_net_g)
  res[[i]]<-data.frame(degree_in=deg_in,degree_out=deg_out,betweenness=bet,closeness=close,eigen_centrality=eigen)
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
  x[,1]
})
deg_inmat<-as.data.frame(t(as.data.frame(deg_inmat)))
colnames(deg_inmat)<-metab
deg_inmat$group<-meta$disease

deg_outmat<-lapply(res_metabolite,function(x){
  x[,2]
})
deg_outmat<-as.data.frame(t(as.data.frame(deg_outmat)))
colnames(deg_outmat)<-metab
deg_outmat$group<-meta$disease

bet_mat<-lapply(res_metabolite,function(x){
  x[,3]
})
bet_mat<-as.data.frame(t(as.data.frame(bet_mat)))
colnames(bet_mat)<-metab
bet_mat$group<-meta$disease

close_mat<-lapply(res_metabolite,function(x){
  x[,4]
})
close_mat<-as.data.frame(t(as.data.frame(close_mat)))
colnames(close_mat)<-metab
close_mat$group<-meta$disease

eig_mat<-lapply(res_metabolite,function(x){
  x[,5]
})
eig_mat<-as.data.frame(t(as.data.frame(eig_mat)))
colnames(eig_mat)<-metab
eig_mat$group<-meta$disease

dif_fun<-function(mat){
  dif_res<-data.frame()
  for(i in 1:(ncol(mat)-1)){
    tmp<-mat[,c(i,ncol(mat))]
    logfc<-log2(mean(tmp[tmp$group=="hypertension",1])/(mean(tmp[tmp$group=="healthy",1])+1e-6))
    p<-wilcox.test(tmp[,1]~tmp[,2])$p.val
    dif_res<-rbind(dif_res,data.frame(row.names = colnames(mat)[i],logfc=logfc,p=p))
  }
  dif_res$metabolite=rownames(dif_res)
  dif_res
}


dif_degin<-dif_fun(deg_inmat)
dif_degout<-dif_fun(deg_outmat)
dif_bet<-dif_fun(bet_mat)
dif_close<-dif_fun(close_mat)
dif_eig<-dif_fun(eig_mat)

top<-10

top<-unique(c(
         dif_bet$metabolite[order(dif_bet$p)[1:top]],
         dif_close$metabolite[order(dif_close$p)[1:top]],
         dif_degin$metabolite[order(dif_degin$p)[1:top]],
         dif_degout$metabolite[order(dif_degout$p)[1:top]],
         dif_eig$metabolite[order(dif_eig$p)[1:top]]))

degin_res<-dif_degin[top,]
degin_res$group<-"in_degree"
degout_res<-dif_degout[top,]
degout_res$group<-"out_degree"
bet_res<-dif_bet[top,]
bet_res$group<-"betweeness"
close_res<-dif_bet[top,]
close_res$group<-"closeness"
eig_res<-dif_eig[top,]
eig_res$group<-"eigen_centrality"
result<-rbind(bet_res,close_res,degin_res,degout_res,eig_res)

top<-c(top,"group")
degin_res<-deg_inmat[,top]
degin_res<-reshape2::melt(degin_res)
degin_res$type<-"In Degree"
degout_res<-deg_outmat[,top]
degout_res<-reshape2::melt(degout_res)
degout_res$type<-"Out Degree"
bet_res<-bet_mat[,top]
bet_res<-reshape2::melt(bet_res)
bet_res$type<-"Betweeness"
close_res<-close_mat[,top]
close_res<-reshape2::melt(close_res)
close_res$type<-"Closeness"
eig_res<-eig_mat[,top]
eig_res<-reshape2::melt(eig_res)
eig_res$type<-"Eigen"
result<-rbind(bet_res,close_res,degin_res,degout_res,eig_res)
symnum.args <- list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", ""))
result$variable<-factor(result$variable,levels = rev(top[-38]))

feeding_box<-ggplot(result, aes(y = log(value + 1), x = variable)) +
  geom_boxplot(aes(color = group), alpha = 1, outlier.colour = NA, position = position_dodge(1)) +
  #geom_jitter(aes(color = group),alpha = 0.5,
  #            position = position_jitterdodge(dodge.width = 1))+
  scale_fill_manual(values = c("hypertension"="#ea5b57", "healthy"="#5b94c2")) +
  scale_color_manual(values = c("hypertension"="#ea5b57", "healthy"="#5b94c2")) +
  theme_bw() +
  labs(x="")+
  facet_grid(.~type, scales = "free") +
  stat_compare_means(aes(group = group), method = "wilcox.test", 
                     symnum.args = symnum.args,size=3,
                     label = "p.signif",vjust = 0.5,hjust=0.8) +
  coord_flip() +
  expand_limits(x = unique(result$variable))+
  theme(text = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size = 12),
        legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        plot.margin = ggplot2::margin(b=0,t=10),
        strip.text = element_text(face = "bold",size = 12))
feeding_box
saveRDS(feeding_box,"output/feeding_box.rds")

top<-100

top<-unique(c(dif_degin$metabolite[order(dif_degin$p)[1:top]],
              dif_degout$metabolite[order(dif_degout$p)[1:top]],
              dif_bet$metabolite[order(dif_bet$p)[1:top]],
              dif_close$metabolite[order(dif_close$p)[1:top]],
              dif_eig$metabolite[order(dif_eig$p)[1:top]]))


metabolite_indegree<-dif_degin[which(dif_degin$p<0.05&abs(dif_degin$logfc)>0.5),]
metabolite_outdegree<-dif_degout[which(dif_degout$p<0.05&abs(dif_degout$logfc)>0.5),]
write.table(metabolite_indegree,"output/feeding_indegree.txt")
write.table(metabolite_outdegree,"output/feeding_outdegree.txt")


feeding<-readRDS("output/feeding_box.rds")
feeding<-feeding+ggtitle("B.")+theme(
  plot.title = element_text(size = 16, face = "bold",vjust = 1,hjust = 0.02),
  plot.title.position = "plot",
  legend.position = "top",
  legend.title = element_blank(),
  plot.margin = ggplot2::margin(r=50)
)


feeding<-readRDS("output/feeding_box.rds")
feeding<-feeding+ggtitle("C.")+
  theme(plot.title = element_text(hjust = 0.05,vjust = -2),
        plot.title.position = "plot",
        title = element_text(size = 16,face="bold"))
pathway<-read.delim("data/ReactionDatabase.txt")
pathway<-pathway[pathway$Subsystem!="",]
subsys<-unique(pathway$Subsystem)
flux<-read.csv("data/discovery_cohort/exchanges.csv")
background<-unique(gsub("^[0-9]+|_m$|\\[e\\]","",unique(flux$metabolite)))
source("script/enrichment_analysis.R")
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

feeding_indegree<-read.table("output/feeding_indegree.txt")
feeding_in_up<-rownames(feeding_indegree)[feeding_indegree$logfc>0]
feeding_in_up<-gsub("\\[e\\]$|_m","",feeding_in_up)
feeding_in_down<-gsub("\\[e\\]$|_m","",rownames(feeding_indegree)[feeding_indegree$logfc<0])

feeding_overpresentation_in_up<-enrichment_analysis(unique(feeding_in_up),background = background,pathway_data = pathwaylist2)
feeding_overpresentation_in_up<-feeding_overpresentation_in_up[feeding_overpresentation_in_up$P_Value<0.05,]
feeding_overpresentation_in_up$change<-"UP"
feeding_overpresentation_in_up$group<-"Top 100 Nodes with the Largest In-degree Differences in the GMxN"
feeding_overpresentation_in_down<-enrichment_analysis(unique(feeding_in_down),background = background,pathway_data = pathwaylist2)
feeding_overpresentation_in_down<-feeding_overpresentation_in_down[feeding_overpresentation_in_down$P_Value<0.05,]
feeding_overpresentation_in_down<-feeding_overpresentation_in_down[-grep("Transport|Others",feeding_overpresentation_in_down$pathway),]
feeding_overpresentation_in_down$change<-"DOWN"
feeding_overpresentation_in_down$group<-"Top 100 Nodes with the Largest In-degree Differences in the GMxN"
feeding_overpresentation_indegree<-rbind(feeding_overpresentation_in_up,feeding_overpresentation_in_down[1:10,])

all_overpres<-feeding_overpresentation_indegree
all_overpres$group<-factor(all_overpres$group,levels = c("Top 100 Nodes with the Largest In-degree Differences in the GMxN"))
all_overpres$pvalue1<- -log10(all_overpres$P_Value)
all_overpres$pvalue1<-ifelse(all_overpres$change=="DOWN",-log10(all_overpres$pvalue1),log10(all_overpres$pvalue1))

feeding_overpres=all_overpres[all_overpres$group=="Top 100 Nodes with the Largest In-degree Differences in the GMxN",]
feeding_enrich<-ggplot(data = feeding_overpres,aes(reorder(pathway, pvalue1), y = pvalue1,fill=change))+
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
  geom_text(data = feeding_overpres[feeding_overpres$change == "UP",], aes(x = pathway, y = -0.01, label = pathway),
            hjust = 1, size = 4) +
  geom_text(data = feeding_overpres[feeding_overpres$change == "DOWN",], aes(x = pathway, y = 0.01, label = pathway),
            hjust = 0, size = 4)+expand_limits(y=c(-1.5,2.5))+
  ggtitle("C.")+theme(
    plot.title = element_text(size = 16, face = "bold",vjust = -1,hjust = 0.01),
    plot.title.position = "plot",
  )



flux<-read.csv("data/discovery_cohort/exchanges.csv",header = T)
flux<-flux[,c(2,7,5)]
flux<-flux %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux<-flux%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
flux[is.na(flux)]<-0
flux<-flux[meta$sample_id,]
#flux<-as.data.frame(t(flux))
group<-meta$disease
flux$group<-group


dif_fun<-function(mat){
  dif_res<-data.frame()
  for(i in 1:(ncol(mat)-1)){
    tmp<-mat[,c(i,ncol(mat))]
    x1<-mean(quantile(tmp[tmp$group=="hypertension",1],probs = seq(0.1,0.9,0.05),na.rm = T))
    x2<-mean(quantile(tmp[tmp$group=="healthy",1],probs = seq(0.1,0.9,0.05),na.rm = T))
    if((x1>0&x2>0)|(x1<0&x2<0)){
      logfc<-log2(x1/x2)
      p<-wilcox.test(tmp[,1]~tmp[,2])$p.val
      if(x1>0&x2>0){
        dif_res<-rbind(dif_res,data.frame(row.names = colnames(mat)[i],logfc=logfc,p=p,direction="output"))
      }else{
        dif_res<-rbind(dif_res,data.frame(row.names = colnames(mat)[i],logfc=logfc,p=p,direction="input"))
      }
    }
    if(x1>0&x2<0){
      p<-wilcox.test(tmp[,1]~tmp[,2])$p.val
      dif_res<-rbind(dif_res,data.frame(row.names = colnames(mat)[i],logfc=5,p=p,direction="input2output"))
    }
    if(x1<0&x2>0){
      p<-wilcox.test(tmp[,1]~tmp[,2])$p.val
      dif_res<-rbind(dif_res,data.frame(row.names = colnames(mat)[i],logfc=5,p=p,direction="output2input"))
    }
    else{dif_res<-rbind(dif_res,data.frame(row.names = colnames(mat)[i],logfc=NA,p=NA,direction="NA"))}
  }
  dif_res$metabolite=rownames(dif_res)
  dif_res
}
tmp<-dif_fun(flux)
tmp<-tmp[!is.na(tmp$p),]

sig_out<-tmp[tmp$p<=0.05&tmp$direction%in%c("output","input2output"),]
sig_in<-tmp[tmp$p<=0.05&tmp$direction%in%c("input","output2input"),]
nosig<-tmp[which(tmp$p>0.05),]
out<-rbind(sig_out,nosig)
out$change=ifelse(out$p<=0.05&abs(out$logfc)>log2(1.3),ifelse(out$logfc>0,"Up","Down"),"Not")
out$change[out$direction=="input2output"&out$change!="Not"]<-"input_to_output"
out$change<-factor(out$change,levels = c("Down","Not","Up","input_to_output"))

inp<-rbind(sig_in,nosig)
inp$change=ifelse(inp$p<=0.05&abs(inp$logfc)>log2(1.3),ifelse(inp$logfc>0,"Up","Down"),"Not")
inp$change[inp$direction=="output2input"&inp$change!="Not"]<-"output_to_input"

out$class<-"Microbial output metabolite"
out$metabolite<-paste(out$metabolite,"")
rownames(out)<-out$metabolite
inp$class<-"Microbial input metabolite"

showlab<-out[abs(out$logfc)>log2(1.2)&out$p<0.05,]
showlab_out<-showlab[grep("MGlcn|^ac_m",showlab$metabolite),]
showlab<-inp[abs(inp$logfc)>log2(1.2)&inp$p<0.05,]
showlab<-showlab[order(inp$p),]
showlab_in<-showlab[grep("MGlcn|^ac_m",showlab$metabolite)[1:15],]
showlab<-rbind(showlab_out,showlab_in)
all_change<-rbind(inp,out)
metabolite_all<-ggplot(all_change, aes(x =logfc, y=-log10(p), colour=change)) +
  geom_point(alpha=0.5, size=2.5) +  #点的透明度和大小
  scale_color_manual(values=c(Down='#5b94c2',Not='gray',Up='#ea5b57',output_to_input="gold2",input_to_output="green4"),
                     labels=c("output_to_input"=expression(H["in"]~HT[out]),"input_to_output"=expression(HT["in"]~H[out]))) + #调整点的颜色
  #xlim(c(-11, 11)) +  ##调整x轴的取值范围，max(abs(BRCA_Match_DEG$logFC))最大值是多少，再进行取舍
  geom_vline(xintercept=c(-log2(1.2),log2(1.2)),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="logFC", y="-log10pvalue") +  #x、y轴标签
  theme_bw()+ylim(c(0,5))+xlim(c(-5,5))+
  facet_row(~class)+
  ggtitle("D.") + #标题
  geom_text_repel(data = subset(showlab, class == "Microbial input metabolite"), 
                  aes(label = metabolite),segment.alpha=0.6,
                  size = 3,box.padding = unit(1, "lines"),
                  point.padding = unit(0, "lines"),vjust=0.7,
                  segment.color = "black",
                  show.legend = FALSE, max.overlaps = 10000) +
  # 针对第二个 facet 使用 geom_text_repel
  geom_text_repel(data = subset(showlab, class == "Microbial output metabolite"), 
                  aes(label = metabolite),size = 3,box.padding = unit(1, "lines"),
                  point.padding = unit(1, "lines"),vjust=0.5,
                  segment.color = "black",
                  show.legend = FALSE, max.overlaps = 10000) +
  theme(axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(size = 14),
        panel.border = element_rect(linewidth = 1),
        panel.grid.major = element_line(linewidth = 1),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 12),
        legend.key.size = unit(1,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 12),
        plot.margin = ggplot2::margin(b=0,t=10),
        plot.title = element_text(size = 16, face = "bold",vjust = -1,hjust = 0.01),
        plot.title.position = "plot",
        plot.subtitle = element_text(size=12,hjust = 0.2)
  )
metabolite_all
c_f<-grid::rasterGrob(png::readPNG("output/corssfeeding.drawio.png"), interpolate = TRUE)
c_f<-ggdraw(c_f) +ggtitle("A.")+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = -0.02,vjust = 1),
        plot.margin = ggplot2::margin(r=30,l=20,b=20,t=10),
        title = element_text(size = 18,face="bold"))
c_f


tmp3<-grid.arrange(c_f,feeding,ncol=1,heights=c(0.3,0.68))
tmp4<-grid.arrange(feeding_enrich,metabolite_all,ncol=1,heights=c(0.45,0.45))
dev.off()

pdf("output/Fig5.pdf",height = 11,width = 17)
wrap_plots(tmp3,tmp4,ncol = 2,widths = c(0.5,0.45))
dev.off()
