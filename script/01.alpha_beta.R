library(curatedMetagenomicData)
library(ggplot2)
library(MicrobiotaProcess)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(tidytree)
library(ggstar)
library(forcats)
library(tidyverse)
library(gridExtra)


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
  returnSamples("relative_abundance",rownames = "short")
abundance_hypertension<-ab_hypertension@assays@data@listData[["relative_abundance"]]
write.csv(abundance_hypertension,"data/abundance_short.csv")

ab_hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance",rownames = "long")
saveRDS(ab_hypertension,"data/rel_ab.rds")
abundance_hypertension<-ab_hypertension@assays@data@listData[["relative_abundance"]]
write.csv(abundance_hypertension,"data/abundance.csv")


ab_hypertension<-readRDS("data/rel_ab.rds")
mpse<-as.MPSE(ab_hypertension)
mpse@assays@data@listData[["Abundance"]]<-mpse@assays@data@listData[["Abundance"]]*1e5
mpse %<>% mp_rrarefy()

##### alpha index and visualization#####
mpse %<>% 
  mp_cal_alpha(.abundance=RareAbundance)

alpha_p<-mpse %>% 
  mp_plot_alpha(
    .group=disease, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou),map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)
  ) +ggtitle("A.")+
  theme(plot.title = element_text(hjust = 0.02,vjust = -0.5),
        plot.title.position = "plot",
        title = element_text(size = 16,face="bold"))+
  scale_fill_manual(values=c("#4f6db7","orangered3"), guide="none") +
  scale_color_manual(values=c("#4f6db7","orangered3"), guide="none")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13,face = "bold"))


#####beta#####
mpse %<>% 
  mp_decostand(.abundance=Abundance)
mpse %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")


mpse %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")

mpse %<>%
  mp_adonis(.abundance=hellinger, .formula=~disease, distmethod="bray", permutations=9999, action="add")
adonis<-mpse %>% mp_extract_internal_attr(name=adonis)

adonis <-paste0("adonis R2: ",round(adonis$R2,2),"; P-value: ", adonis$`Pr(>F)`)

beta_p<-mpse %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = disease, 
    .color = disease, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = F # don't display the legend of stat_ellipse 
  ) +labs(subtitle = adonis)+
  scale_alpha(breaks=c(1,2,3,4))+
  scale_fill_manual(
    values = c("#4f6db7","orangered3"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16))
  ) +
  scale_color_manual(
    values=c("#4f6db7","orangered3"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16))
  ) +
  scale_size_continuous(
    limits =c(50, 250),breaks = c(50,100,150,200,250),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16))
  )+
  scale_alpha_continuous(limits=c(1,4),breaks=c(1,2,3,4),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16)))+
  theme(axis.text.x = element_text(size = 14,face = "bold"),
          legend.spacing.y = unit(0.03, "cm"),
          #legend.position = "none",
          axis.text.y = element_text(size = 14,face = "bold"))+
  theme(axis.title.x = element_text(margin = margin(t = -40)),
        axis.title.y = element_text(margin = margin(r = -40)),
        plot.margin = ggplot2::margin(b=30,l=30))+
  guides(
    fill = guide_legend(order = 3, keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 15)),    # 调整fill图例顺序
    color = guide_legend(order = 2, keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 15)),   # 调整color图例顺序
    size = guide_legend(order = 1, keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 15))     # 调整size图例顺序
  )+
  ggtitle("C.")+
  theme(plot.title = element_text(hjust = -0.02,vjust = -0.5),
        plot.title.position = "plot",legend.position = "none",
        plot.subtitle = element_text(hjust=0.25,size=12),
        plot.margin = ggplot2::margin(r=40,l=30,t=10,b=25),
        title = element_text(size = 16,face="bold"))


####test data
otu<-read.table("data/test_cohort/test_16s/microbio_selected.otus",row.names = 1,header = T)
meta_test<-read.csv("data/test_cohort//test_16s/SraRunTable.csv",header = T)
meta_healthy<-meta_test[meta_test$AGE>35&
                          meta_test$BMI>18.5&meta_test$BMI<25&
                          meta_test$Cardiometabolic_status=="Healthy"&
                          (meta_test$diastolic_bp<80&meta_test$Systolic_BP<130)&
                          meta_test$Medicament_use=="No"
                        #meta_test$Total_Cholesterol<200
                        ,]
meta_hy<-meta_test[meta_test$AGE>35&
                     #meta_test$BMI>18.5&meta_test$BMI<25&
                     (meta_test$diastolic_bp>90&meta_test$Systolic_BP>140)&
                     #meta_test$Total_Cholesterol>200&
                     meta_test$Medicament_use=="No"
                   ,]

meta_test<-rbind(meta_healthy,meta_hy)
rownames(meta_test)<-meta_test$Sample.Name
meta_test$group<-"median"
meta_test$group[meta_test$diastolic_bp<90&meta_test$Systolic_BP<130]<-"healthy"
meta_test$group[meta_test$diastolic_bp>90&meta_test$Systolic_BP>140]<-"hypertension"
meta_test<-meta_test[meta_test$group!="median",]
write.csv(meta_test,"data/test_cohort/test_16s/meta_test_select.csv")
taxa<-read.table("data/test_cohort/test_16s/microbio_selected.taxonomy",row.names = 1,header = T)
taxatmp<-strsplit(taxa$Taxonomy,";")
taxatmp<-as.data.frame(taxatmp)
taxa<-t(taxatmp)
rownames(taxa)<-colnames(otu)
colnames(taxa)<-c("Kingdom", "Phylum","Class","Order","Family" ,"Genus","Species")
taxa[,1]<-gsub("k__","",taxa[,1])
taxa[,2]<-gsub("p__","",taxa[,2])
taxa[,3]<-gsub("c__","",taxa[,3])
taxa[,4]<-gsub("o__","",taxa[,4])
taxa[,5]<-gsub("f__","",taxa[,5])
taxa[,6]<-gsub("g__","",taxa[,6])
taxa[,7]<-gsub("s__","",taxa[,7])
taxa[taxa=="unclassified"]<-NA

test_mpse<-mp_import_dada2(seqtab = otu,taxatab = taxa,sampleda = meta_test)
test_mpse %<>% mp_rrarefy()

test_mpse %<>% 
  mp_cal_alpha(.abundance=RareAbundance)

test_alpha_p<-test_mpse %>% 
  mp_plot_alpha(
    .group=group, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)
  ) +
  scale_fill_manual(values=c("#4f6db7","orangered3"), guide="none") +
  scale_color_manual(values=c("#4f6db7","orangered3"), guide="none")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13,face = "bold"))+
  ggtitle("B.")+
  theme(plot.title = element_text(hjust = 0,vjust = -0.5),
        plot.title.position = "plot",
        title = element_text(size = 16,face="bold"))

#####test_beta
test_mpse %<>% 
  mp_decostand(.abundance=Abundance)
test_mpse %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")


test_mpse %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
test_mpse %<>%
  mp_adonis(.abundance=hellinger, .formula=~group, distmethod="bray", permutations=9999, action="add")
adonis<-test_mpse %>% mp_extract_internal_attr(name=adonis)

adonis <-paste0("adonis R2: ",round(adonis$R2,2),"; P-value: ", adonis$`Pr(>F)`)

test_beta_p<-test_mpse %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = group, 
    .color = group, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +labs(subtitle = adonis)+
  scale_alpha(breaks=c(1,2,3,4))+
  scale_fill_manual(
    values = c("#4f6db7","orangered3"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16))
  ) +
  scale_color_manual(
    values=c("#4f6db7","orangered3"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16))
  ) +
  scale_size_continuous(
    limits = c(50, 250),breaks = c(50,100,150,200,250),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16))
  )+
  scale_alpha_continuous(limits=c(1,4),breaks=c(1,2,3,4),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=16)))+
  theme(axis.text.x = element_text(size = 14,face = "bold"),
          legend.spacing.y = unit(0.03, "cm"),
          #legend.position = "none",
          axis.text.y = element_text(size = 16,face = "bold"))+
  theme(axis.title.x = element_text(margin = margin(t = -40)),
        axis.title.y = element_text(margin = margin(r = -40)),
        plot.margin = ggplot2::margin(b=30,l=30))+
  guides(
    fill = guide_legend(order = 3, keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 15)),    # 调整fill图例顺序
    color = guide_legend(order = 2, keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 15)),   # 调整color图例顺序
    size = guide_legend(order = 1, keywidth = 0.6, keyheight = 0.6, label.theme = element_text(size = 15))     # 调整size图例顺序
  )+
  ggtitle("D.")+
  theme(plot.title = element_text(hjust = 0.07,vjust = -0.5),
        plot.subtitle = element_text(hjust=0.2,size=13),
        plot.title.position = "plot",
        plot.margin = ggplot2::margin(r=30,l=30,t=10,b=25),
        title = element_text(size = 16,face="bold"))


p_all1<-grid.arrange(alpha_p,test_alpha_p,ncol=2)
p_all2<-grid.arrange(beta_p,test_beta_p,ncol=2,widths=c(0.4,0.5))
dev.off()
pdf("output/Fig1.pdf",width = 16,height = 10)
grid.arrange(
  p_all1, 
  grid::nullGrob(),  # 插入一个空白的 grob
  p_all2, 
  heights = unit(c(2, 0.5, 4), c("null", "cm", "null"))  # 调整两图间的空白高度为 1 cm
)
dev.off()


rel_ab<-mpse@assays@data@listData[["RelRareAbundanceBySample"]]

library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
flux<-read.csv("data/discovery_cohort/exchanges.csv")
flux<-flux[,c(2,7,5)]
flux<-flux %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux<-flux%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
meta<-read.csv("data/meta.csv",row.names = 1)
flux<-flux[meta$sample_id,]
flux[is.na(flux)]<-0
flux<-flux[meta$sample_id,]
flux[flux!=0]<-1

ab<-read.csv("data/abundance_short.csv",row.names = 1)
ab_sample<-as.character(read.csv("data/abundance_short.csv",row.names = 1,header = F)[1,])
colnames(ab)<-ab_sample
ab<-ab[,meta$sample_id]
ab[ab!=0]<-1
data<-cbind(t(ab),flux)



flux_test<-read.csv("data/test_cohort/test_16s/res_grow_g/exchanges.csv")
flux_test<-flux_test[,c(2,7,5)]
flux_test<-flux_test %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux_test<-flux_test%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
meta_test<-read.csv("data/test_cohort/test_16s/meta_test_select.csv",header = T,row.names = 1)
flux_test<-flux_test[meta_test$Sample.Name,]
flux_test[is.na(flux_test)]<-0
flux_test[flux_test!=0]<-1

test_ab<-read.csv("data/test_cohort/test_16s/species_16s.csv",row.names = 1,header = T)
test_ab<-as.data.frame(t(test_ab))
test_ab<-test_ab[meta_test$Sample.Name,]
test_ab[test_ab!=0]<-1
test_ab<-test_ab[,colSums(test_ab)>0]
data_test<-cbind(test_ab,flux_test)


flux<-as.data.frame(t(flux))
flux_test<-as.data.frame(t(flux_test))
flux_all<-merge(flux,flux_test,by="row.names",all=T)
rownames(flux_all)<-flux_all$Row.names
flux_all<-flux_all[,-1]
flux_all[is.na(flux_all)]<-0

ab<-as.data.frame(ab)
ab_test<-as.data.frame(t(test_ab))
ab_all<-merge(ab,ab_test,by="row.names",all=T)
ab_all[is.na(ab_all)]<-0
rownames(ab_all)<-ab_all$Row.names
ab_all<-ab_all[,-1]

meta<-data.frame(row.names = c(meta$sample_id,meta_test$Sample.Name),sample=c(meta$sample_id,meta_test$Sample.Name),group=c(meta$disease,meta_test$group))
meta_all<-rbind(meta[meta$group=="hypertension",],meta[meta$group=="healthy",])
ab_all<-ab_all[,meta$sample]
flux_all<-flux_all[,meta$sample]
data_all<-t(rbind(flux_all,ab_all))

row_ha = rowAnnotation(disease = meta_all$group,col = list(disease = c("healthy"="#4f6db7","hypertension"="orangered4")),
                       annotation_legend_param = list(
                         disease = list(
                           title = "disease",
                           title_gp = gpar(fontsize = 14),
                           labels_gp = gpar(fontsize = 14))
                       ),
                       annotation_name_side ="top"
)

col_ha<-HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c("#00A087FF","#3C5488FF")),
                                           labels = c("Metabolite flux", "Species"), 
                                           labels_gp = gpar(col = "grey95", fontsize = 18)))
row_ann<-rowAnnotation(foo = anno_block(gp = gpar(fill = c("#4f6db7","orangered4")),
                                        labels = c("Healthy", "Hypertension"), 
                                        labels_gp = gpar(col = "grey95", fontsize = 14)))

p_heat_all<-Heatmap(data_all, name = "Existence", cluster_rows = F,cluster_columns = F,
                     column_split = factor(rep(c(" ", "  "), c(nrow(flux_all),nrow(ab_all))),levels = c(" ","  ")),
                     row_split = factor(rep(c(" ", "  "), c(sum(meta_all$group=="hypertension"),sum(meta_all$group=="healthy"))),levels = c("  "," ")),
                     left_annotation = row_ann,
                     top_annotation = col_ha,
                     show_heatmap_legend = T,
                     col=c("grey94","lightsteelblue3"),
                     heatmap_legend_param = list(
                       title_gp = gpar(fontsize = 14),
                       labels_gp = gpar(fontsize = 14)),
                     column_title_gp = gpar(fontsize=18),
                     show_row_names = F,show_column_names = F)
saveRDS(p_heat_all,"output/sparse.rds")
