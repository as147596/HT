library(tidyverse)
library(ggplot2)
library(reticulate)
library(cowplot)
library(patchwork)
library(microbiomeMarker)
library(phyloseq)
library(curatedMetagenomicData)
library(gridExtra)
library(balance)
# key_genus ----

flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv",header = T)
key_sp<-read.csv("output/cluster_3.csv")
flux_sel<-flux[flux$taxon%in%key_sp$x,]
flux_sel$metabolite<-gsub("\\[e\\]","",flux_sel$metabolite)
importance<-read.csv("output/doubleml_res.csv",row.names = 1)
importance$X<-gsub("_m$","",importance$X)
flux_sel<-flux_sel[flux_sel$metabolite%in%importance$X,]
flux_sel$sp_m<-paste(flux_sel$taxon,flux_sel$metabolite,sep = "_")
flux_sel$flux<-flux_sel$flux*flux_sel$abundance
flux_sel<-flux_sel[,c(2,9,5)]
flux_sel<-flux_sel %>%group_by(sample_id,sp_m)%>%
  summarise(flux=sum(flux))
flux_sel<-flux_sel%>% pivot_wider(names_from = sp_m,values_from = flux) %>% column_to_rownames("sample_id")
flux_sel[is.na(flux_sel)]<-0

write.csv(flux_sel,"data/key_sp_flux/key_sp_flux.csv",quote = F)

reticulate::use_condaenv("F:/software/anaconda/envs/doubleml/")



### step 2 ----

library(tidyverse)
library(ggplot2)
library(reticulate)
library(cowplot)
library(patchwork)
library(igraph)
library(ggraph)

ATE<-read.csv("output/keysp_doubleml_Systolic_.csv")
ATE$padj<-p.adjust(ATE$pvalue,method = "BH")
sign<-sign(ATE[,c("ATE","CI2.5","CI97.5","theta.lower","theta.upper")])
keep<-apply(sign, 1, function(x){
  length(unique(x))==1
})
sig<-ATE[ATE$padj<0.05&ATE$RV>0.05,]
sig$color<-ifelse(sig$ATE>0,"R","P")
sig$metabolite<-factor(sig$X,levels = sig$X[order(sig$ATE)])
sig_SBP<-sig

ATE<-read.csv("output/keysp_doubleml_diastolic.csv")
ATE$padj<-p.adjust(ATE$pvalue,method = "BH")
sign<-sign(ATE[,c("ATE","CI2.5","CI97.5","theta.lower","theta.upper")])
keep<-apply(sign, 1, function(x){
  length(unique(x))==1
})
sig<-ATE[ATE$padj<0.05&ATE$RV>0.05,]
sig$color<-ifelse(sig$ATE>0,"R","P")
sig$metabolite<-factor(sig$X,levels = sig$X[order(sig$ATE)])
sig_DBP<-sig

sig_DBP$bp<-"DBP"
sig_SBP$bp<-"SBP"
x<-sub("^[^_]*_", "", sig_SBP$metabolite)
importance<-read.csv("output/doubleml_res.csv")
imp_up<-gsub("_m$","",importance$X)[importance$ATE>0]|>unique()
imp_down<-gsub("_m$","",importance$X)[importance$ATE<0]|>unique()
sig_SBP<-sig_SBP[x%in%gsub("_m$","",importance$X)[importance$class=="SBP"],]
sig<-rbind(sig_DBP,sig_SBP)

x<-sub("_",";",sig$metabolite)
sig$sp<-gsub(";.*","",x)
sig$metabolite<-gsub(".*;","",x)

sig<-sig[(sig$metabolite%in%imp_up&sig$ATE>0)|
           (sig$metabolite%in%imp_down&sig$ATE<0),]


selected=read.csv("data/res_grow_wd_pFBA/selected_sample.csv")
meta<-read.csv("data/discovery_cohort/16s/SraRunTable.csv",row.names = 1)
meta<-meta[meta$Sample.Name%in%selected$x,]
meta$disease<-ifelse(meta$diastolic_bp>=90|meta$Systolic_BP>=140,"hypertension","healthy")
dir<-read.csv("data/res_grow_wd_pFBA/exchanges.csv")
dir<-dir[dir$sample_id%in%selected$x,]
dir<-dir[dir$taxon%in%sig$sp,]
dir$metabolite<-gsub("\\[e\\]","",dir$metabolite)
dir<-dir[dir$metabolite%in%sig$metabolite,]
dir$id<-paste(dir$taxon,gsub("\\[e\\]","",dir$metabolite),sep=";")
direction<-c()
for(i in 1:length(x)){
  tmp<-dir[dir$id==x[i],]
  d<-tmp[tmp$sample_id%in%meta$Sample.Name[meta$disease=="hypertension"],]
  h<-tmp[tmp$sample_id%in%meta$Sample.Name[meta$disease=="healthy"],]
  d<-table(d$direction)[which.max(table(d$direction))]|>names()
  h<-table(h$direction)[which.max(table(h$direction))]|>names()
  
}

edge<-data.frame(sp=sig$sp,metabolite=sig$metabolite,
                 ATE=sig$ATE,padj=sig$padj,BP=sig$bp,RV=sig$RV)
#nodes<-data.frame(nodes=unique(c(sig$sp,sig$metabolite)),
#                  type=rep(c("Genus","Metabolites"),
#                           c(length(unique(sig$sp)),length(unique(sig$metabolite)))))
#
#g<-graph_from_data_frame(edge,vertices = nodes)
edge$signif<-case_when(edge$padj<0.001~"***",
                       edge$padj>=0.001&edge$padj<0.01~"**",
                       edge$padj>=0.01&edge$padj<0.05~"*",
                       .default = ""
)
edge$metabolite<-factor(edge$metabolite,levels = c(imp_down,imp_up))

p1<-ggplot(edge,aes(sp,metabolite))+
  geom_point(aes(color=ATE,size = RV))+
  facet_grid(BP~.,scale="free_y",space="free_y")+
  scale_color_gradient2(low = "#4575B4",mid = "white",high = "#D73027")+
  geom_text(aes(label = signif),hjust=0.5,vjust=0.65)+
  theme_test()+labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.text.y = element_text(size = 10))+ggtitle("A.")+
  theme(axis.text = element_text(color = "black",size = 10),
        plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = 0,face = "bold"),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 14),
        legend.key.size = unit(0.8,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = ggplot2::margin(b=0,t=10))


# val key_genus 2 ----

flux<-read.csv("data/res_grow_wd_pFBA/exchanges.csv",header = T)
load("data/sampleMetadata.rda")
hypertension_meta<-sampleMetadata[which(sampleMetadata$disease=='hypertension'),]
healthy_meta<-hypertension_meta[which(hypertension_meta$age_category=='adult'),]
healthy_meta<-sampleMetadata[which(sampleMetadata$disease=='healthy'),]
healthy_meta<-healthy_meta[which(healthy_meta$age_category=='adult'),]
healthy_meta<-healthy_meta[which(healthy_meta$country=='CHN'|healthy_meta$country=='ITA'|healthy_meta$country=='AUT'),]
healthy_meta<-healthy_meta[which(healthy_meta$BMI>18.5&healthy_meta$BMI<23.9),]
healthy_meta<-healthy_meta[which(healthy_meta$age>55&healthy_meta$age<70),]
hypertension_meta<-rbind(hypertension_meta,healthy_meta)

ab_hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance",rownames = "long")

abund_hypertension<-ab_hypertension@assays@data@listData[["relative_abundance"]]|>as.data.frame()
genus<-sapply(rownames(abund_hypertension),function(x){
  gsub("^g__","",strsplit(x,"\\|")[[1]][6])
})|>unname()
abund_hypertension$genus<-genus
genus_df<-aggregate(.~genus,abund_hypertension,sum)
rownames(genus_df)<-genus_df[,1]
genus_df<-genus_df[,-1]
hypertension_meta$sample_id<-gsub("-","_",hypertension_meta$sample_id)
hypertension_meta$sample_id[hypertension_meta$sample_id=="SID21.TO.V.new"]="SID21_TO_V_new"
inter_sample<-intersect(hypertension_meta$sample_id,flux$sample_id)
flux<-flux[flux$sample_id%in%inter_sample,]
genus<-flux$taxon|>unique()
genus_df<-genus_df[genus[-length(genus)],]
keep<-colSums(genus_df)>80
genus_df<-genus_df[,keep]
colnames(genus_df)<-gsub("-","_",colnames(genus_df))
colnames(genus_df)[colnames(genus_df)=="SID21.TO.V.new"]="SID21_TO_V_new"
hypertension_meta<-hypertension_meta[hypertension_meta$sample_id%in%colnames(genus_df),]
flux<-flux[flux$sample_id%in%hypertension_meta$sample_id,]

key_sp<-read.csv("output/cluster_3.csv")
flux_sel<-flux[flux$taxon%in%key_sp$x,]
flux_sel$metabolite<-gsub("\\[e\\]","",flux_sel$metabolite)
importance<-read.csv("output/doubleml_res.csv",row.names = 1)
importance$X<-gsub("_m$","",importance$X)
flux_sel<-flux_sel[flux_sel$metabolite%in%importance$X,]
flux_sel$sp_m<-paste(flux_sel$taxon,flux_sel$metabolite,sep = ";")
flux_sel$flux<-flux_sel$flux*flux_sel$abundance
flux_sel<-flux_sel[,c(2,9,5)]
flux_sel<-flux_sel %>%group_by(sample_id,sp_m)%>%
  summarise(flux=sum(flux))
flux_sel<-flux_sel%>% pivot_wider(names_from = sp_m,values_from = flux) %>% column_to_rownames("sample_id")
flux_sel[is.na(flux_sel)]<-0
flux_sel<-flux_sel[hypertension_meta$sample_id,]
dif_res<-data.frame()
for(i in 1:ncol(flux_sel)){
  p<-wilcox.test(as.numeric(flux_sel[,i][hypertension_meta$disease=="hypertension"]),
                 as.numeric(flux_sel[,i][hypertension_meta$disease!="hypertension"]))$p.value
  tmp<-data.frame(feature=colnames(flux_sel)[i],
                  logfc=(mean(flux_sel[,i][hypertension_meta$disease=="hypertension"])-
                           mean(flux_sel[,i][hypertension_meta$disease!="hypertension"])),
                  pvalue=p)
  dif_res<-rbind(dif_res,tmp)
}
sig<-dif_res[dif_res$pvalue<0.05,]
sig$metabolite<-gsub(".*;","",sig$feature)
sig<-sig[(sig$metabolite%in%imp_up&sig$logfc>0)|
           (sig$metabolite%in%imp_down&sig$logfc<0),]
sig$sp<-gsub(";.*","",sig$feature)
sig$signif<-case_when(sig$pvalue<0.001~"***",
                      sig$pvalue>=0.001&sig$pvalue<0.01~"**",
                      sig$pvalue>=0.01&sig$pvalue<0.05~"*",
                      .default = ""
)
sig$metabolite<-factor(sig$metabolite,levels = c(imp_down,imp_up))
sig$logFC<-ifelse(sig$logfc>0,"logFC>0","logFC<0")
p2=ggplot(sig,aes(sp,metabolite))+
  geom_point(aes(color=logFC,size = -log10(pvalue)))+
  scale_color_manual(values = c("#4575B4", "#D73027"))+
  geom_text(aes(label = signif),hjust=0.5,vjust=0.65)+
  theme_test()+labs(x="",y="")+ggtitle("B.")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 10),
        axis.text.y = element_text(size = 10))+
  theme(axis.text = element_text(color = "black",size = 10),
        plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = 0,face = "bold"),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 14),
        legend.key.size = unit(0.8,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = ggplot2::margin(b=0,t=10))

inter_metab<-list(discovery=edge$metabolite,
                  validation=sig$metabolite)
p3<-ggvenn::ggvenn(inter_metab)+
  labs(title = "C.")+coord_cartesian(clip = "off")+
  theme(plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = 3,face = "bold"))

inter_sp<-list(discovery=edge$sp,
                  validation=sig$sp)
p4<-ggvenn::ggvenn(inter_sp)+coord_cartesian(clip = "off")+
  theme(plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = 0,face = "bold"))




# balance ----

vol_prepare<-function(vol,key_genus){
  key_genus_ind<-sapply(vol$z,function(x){
    any(strsplit(x,"_and_")[[1]]%in%key_genus$x)
  })
  vol$sign<-ifelse(vol$pvalue>0.05,"Not Significant",
                   ifelse(vol$pvalue<=0.05&key_genus_ind,"Imbalanced genus within\nthe FERM guild","Other imbalance"))
  vol
}
otu<-read.table("data/discovery_cohort/16s/microbio_selected.otus",row.names = 1,header = T)
meta_test<-read.csv("data/discovery_cohort/16s/SraRunTable.csv",header = T)
rownames(meta_test)<-meta_test$Sample.Name
otu<-otu[rownames(meta_test),]
meta_test$group<-ifelse(meta_test$diastolic_bp>=90|meta_test$Systolic_BP>=140,"hypertension","healthy")
taxa<-read.table("data/discovery_cohort/16s/microbio_selected.taxonomy",row.names = 1,header = T)
taxatmp<-strsplit(taxa$Taxonomy,";")
taxatmp<-as.data.frame(taxatmp)
taxa<-t(taxatmp)
rownames(taxa)<-colnames(otu)

otu<-otu|>t()|>as.data.frame()
otu$genus<-taxa[,6]

genus<-aggregate(.~genus,otu,sum)
rownames(genus)<-gsub("g__","",genus$genus)
genus<-genus[,-1]
genus<-prop.table(as.matrix(genus),2)|>as.data.frame()


## Discovery balance ----

key_genus<-read.csv("output/cluster_3.csv")
sbp<-sbp.fromADBA(t(genus),meta_test$group)
sbp<-sbp.subset(sbp)
z<-balance.fromSBP(t(genus),sbp)
pvalue<-c()
logfc<-c()
for(i in 1:ncol(z)){
  p=wilcox.test(z[,i]~meta_test$group)$p.value
  pvalue<-c(pvalue,p)
  logfci<-mean(z[meta_test$group=="hypertension",i])-
    mean(z[meta_test$group!="hypertension",i])
  logfc<-c(logfc,logfci)
}
sum(pvalue<0.05,na.rm = T)
vol<-data.frame(z=colnames(sbp),pvalue=pvalue,logfc=logfc)
vol1<-vol_prepare(vol,key_genus = key_genus)
sum(pvalue<0.05)
z<-z[,pvalue<0.05]|>t()|>as.data.frame()
colnames(z)<-meta_test$Sample.Name
z$pvalue<-pvalue[pvalue<0.05]
z$logfc<-logfc[pvalue<0.05]
write.csv(z,"output/imbalance_dis.csv",quote = F)
tmp<-colnames(sbp)[pvalue<0.05]
tmp<-sapply(tmp,function(x){
  strsplit(x,"_and_")[[1]]
})
tmp<-unlist(tmp)|>unique()
prop<-data.frame(group=c("FERM guild","Other guilds"),
                 value=c(sum(tmp%in%key_genus$x)/19,
                         sum(!tmp%in%key_genus$x)/(430-19)))
prop$value<-prop$value*100
prop$label<-paste0(round(prop$value,2),"%")
prop1<-prop
prop1$cohort<-"Discovery"

## Validation balance ----

### data prepare ----

load("data/sampleMetadata.rda")
hypertension_meta<-sampleMetadata[which(sampleMetadata$disease=='hypertension'),]
healthy_meta<-hypertension_meta[which(hypertension_meta$age_category=='adult'),]
healthy_meta<-sampleMetadata[which(sampleMetadata$disease=='healthy'),]
healthy_meta<-healthy_meta[which(healthy_meta$age_category=='adult'),]
healthy_meta<-healthy_meta[which(healthy_meta$country=='CHN'|healthy_meta$country=='ITA'|healthy_meta$country=='AUT'),]
healthy_meta<-healthy_meta[which(healthy_meta$BMI>18.5&healthy_meta$BMI<23.9),]
healthy_meta<-healthy_meta[which(healthy_meta$age>55&healthy_meta$age<70),]
hypertension_meta<-rbind(hypertension_meta,healthy_meta)

ab_hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance",rownames = "long")

abund_hypertension<-ab_hypertension@assays@data@listData[["relative_abundance"]]|>as.data.frame()
genus<-sapply(rownames(abund_hypertension),function(x){
  gsub("^g__","",strsplit(x,"\\|")[[1]][6])
})|>unname()
abund_hypertension$genus<-genus
genus_df<-aggregate(.~genus,abund_hypertension,sum)
rownames(genus_df)<-genus_df[,1]
genus_df<-genus_df[,-1]
genus_df<-genus_df[,hypertension_meta$sample_id]

sbp<-sbp.fromADBA(t(genus_df),hypertension_meta$disease)
sbp<-sbp.subset(sbp)
z1<-balance.fromSBP(t(genus_df),sbp)
rownames(z1)<-hypertension_meta$sample_id
pvalue<-c()
logfc<-c()
for(i in 1:ncol(z1)){
  p=wilcox.test(z1[,i]~hypertension_meta$disease)$p.value
  pvalue<-c(pvalue,p)
  logfci<-mean(z1[hypertension_meta$disease=="hypertension",i])-
    mean(z1[hypertension_meta$disease!="hypertension",i])
  logfc<-c(logfc,logfci)
}
vol<-data.frame(z=colnames(sbp),pvalue=pvalue,logfc=logfc)
vol2<-vol_prepare(vol,key_genus = key_genus)
sum(pvalue<0.05)
z1<-t(z1[,pvalue<0.05])|>as.data.frame()
z1$pvalue<-pvalue[pvalue<0.05]
z1$logfc<-logfc[pvalue<0.05]
write.csv(z1,"output/imbalance_val.csv",quote = F)
tmp<-colnames(sbp)[pvalue<0.05]
tmp<-sapply(tmp,function(x){
  strsplit(x,"_and_")[[1]]
})
tmp<-unlist(tmp)|>unique()
prop2<-data.frame(group=c("FERM guild","Other guilds"),
                  value=c(sum(tmp%in%key_genus$x)/sum(key_genus$x%in%rownames(genus_df)),
                          sum(!tmp%in%key_genus$x)/(nrow(genus_df)-sum(key_genus$x%in%rownames(genus_df)))))
prop2$value<-prop2$value*100
prop2$label<-paste0(round(prop2$value,2),"%")
prop2$cohort<-"Validation"


## merge and plot ----

vol1$cohort<-"Discovery"
vol2$cohort<-"Validation"
vol<-rbind(vol1,vol2)
p5<-ggplot(vol,aes(logfc,-log10(pvalue),color=sign))+
  geom_point(size=2)+
  facet_grid(.~cohort,scale="free_x")+
  scale_color_manual(name="",values = c("red","grey","mediumorchid"))+
  theme_test()+ggtitle("D.")+
  labs(x="logFC",y="-log10(Pvalue)")+
  theme(axis.text = element_text(color = "black",size = 10),
        plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = 0,face = "bold"),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 14),
        legend.key.size = unit(0.8,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = ggplot2::margin(b=3,t=0,l=2))

prop<-rbind(prop1,prop2)
p6<-ggplot(prop,aes(value,group))+
  geom_col(fill="#aa6190",alpha=0.8)+
  geom_text(aes(x=value+1,label = label))+
  theme_classic()+xlim(c(0,90))+
  labs(y="",x="Percentage of genus with balance imbalance")+
  facet_grid(cohort~.)+ggtitle("E.")+
  theme(axis.text = element_text(color = "black",size = 10),
        plot.title.position = "plot",
        plot.title = element_text(size=14,vjust = 0,face = "bold"),
        strip.background = element_rect(linewidth = 1),
        strip.text = element_text(face = "bold",size = 12),
        legend.key.size = unit(0.8,'cm'), 
        legend.title = element_text(size = 14,face = "bold"),
        legend.text = element_text(size = 10),
        plot.margin = ggplot2::margin(b=3,t=0,r=5))


p1_4<-grid.arrange(p1,plot_grid(p2,(p3+p4), nrow=2,rel_heights = c(0.7,0.3)),nrow=1)
p5_6<-grid.arrange(p5,p6,nrow=1,widths=c(0.6,0.4))

dev.off()
pdf("output/fig5.pdf",height = 8,width = 10)
grid.arrange(p1_4,p5_6,nrow=2,heights=c(0.65,0.3))
dev.off()
