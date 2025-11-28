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
library(microbiomeMarker)

otu<-read.table("data/test_cohort/test_16s/microbio_selected.otus",row.names = 1,header = T)
meta_test<-read.csv("data/test_cohort//test_16s/SraRunTable.csv",header = T)
#meta_healthy<-meta_test[meta_test$AGE>35&
#                          meta_test$BMI>18.5&meta_test$BMI<25&
#                          meta_test$Cardiometabolic_status=="Healthy"&
#                          (meta_test$diastolic_bp<80&meta_test$Systolic_BP<130)&
#                          meta_test$Medicament_use=="No"
#                        #meta_test$Total_Cholesterol<200
#                        ,]
#meta_hy<-meta_test[meta_test$AGE>35&
#                     #meta_test$BMI>18.5&meta_test$BMI<25&
#                     (meta_test$diastolic_bp>90&meta_test$Systolic_BP>140)&
#                     #meta_test$Total_Cholesterol>200&
#                     meta_test$Medicament_use=="No"
#                   ,]
#
#meta_test<-rbind(meta_healthy,meta_hy)
rownames(meta_test)<-meta_test$Sample.Name
meta_test$group<-ifelse(meta_test$diastolic_bp>=90|meta_test$Systolic_BP>=140,"hypertension","healthy")

meta_test<-meta_test[meta_test$group!="median",]
#write.csv(meta_test,"data/test_cohort/test_16s/meta_test_select.csv")
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
alpha_df<-test_mpse@colData@listData|>as.data.frame()
alpha_df<-alpha_df[,56:61]|>as.data.frame()
rownames(alpha_df)<-test_mpse$Sample.Name
write.csv(alpha_df,"output/alpha_df.csv",)

test_alpha_p<-test_mpse %>% 
  mp_plot_alpha(
    .group=group, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05)
  ) +labs(y="Alpha Index")+
  scale_fill_manual(values=c("#4f6db7","orangered3"), guide="none") +
  scale_color_manual(values=c("#4f6db7","orangered3"), guide="none")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13,face = "bold"))+
  ggtitle("A.")+
  theme(plot.title = element_text(hjust = 0,vjust = -0.5),
        plot.title.position = "plot",
        title = element_text(size = 16,face="bold"))

#####test_beta
#pvalue<-c()
#otu<-otu[meta_test$Sample.Name,]
#for(i in 1:nrow(otu)){
#  p<-wilcox.test(otu[meta_test$group=="healthy",i],otu[meta_test$group=="hypertension",i])$p.value
#  pvalue<-c(pvalue,p)
#}
#
test_mpse1<-test_mpse#[which(pvalue<0.01),]
test_mpse1 %<>% 
  mp_decostand(.abundance=Abundance)
test_mpse1 %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")


test_mpse1 %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
test_mpse1 %<>%
  mp_adonis(.abundance=hellinger, .formula=~group, distmethod="bray", permutations=9999, action="add")
adonis<-test_mpse1 %>% mp_extract_internal_attr(name=adonis)

adonis <-paste0("adonis R2: ",round(adonis$R2,2),"; P-value: ", adonis$`Pr(>F)`)

test_beta_p<-test_mpse1 %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = group, 
    .color = group, 
    #.size = Chao1, 
    .alpha = ACE,
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
  theme(axis.text.x = element_text(size = 12,face = "bold"),
        legend.spacing.y = unit(0.03, "cm"),
        #legend.position = "bottom",
        axis.text.y = element_text(size = 12,face = "bold"))+
  theme(axis.title.x = element_text(margin = margin(t = -40)),
        axis.title.y = element_text(margin = margin(r = -40)))+
  ggtitle("B.")+
  theme(plot.title = element_text(hjust = -0.08,vjust = 3),
        plot.subtitle = element_text(hjust=0.2,size=13),
        plot.title.position = "plot",#legend.direction = "vertical",
        plot.margin = ggplot2::margin(r=10,l=40,t=20,b=30),
        title = element_text(size = 16,face="bold"),
        legend.background = element_blank())


p_all1<-grid.arrange(test_alpha_p,test_beta_p,ncol=1,heights=c(0.3,0.7))

dev.off()
#pdf("output/Fig1.pdf",width = 16,height = 10)
#grid.arrange(
#  p_all1, 
#  grid::nullGrob(),  # 插入一个空白的 grob
#  p_all2, 
#  heights = unit(c(2, 0.5, 4), c("null", "cm", "null"))  # 调整两图间的空白高度为 1 cm
#)
#dev.off()

phy<-as.phyloseq(test_mpse)
mm_lefse <- run_lefse(
  phy,
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2,
  taxa_rank = "Genus"
)
write.csv(data.frame(mm_lefse@marker_table),"output/genus_dif.csv",row.names = F,quote = F)

mm_lefse@marker_table$lda<-ifelse(mm_lefse@marker_table$enrich_group=="healthy",
                                  -mm_lefse@marker_table$ef_lda,
                                  mm_lefse@marker_table$ef_lda)
marker_table<-mm_lefse@marker_table
lefse_bar<-ggplot(marker_table,aes(x=lda,y=reorder(feature,lda),fill=enrich_group))+
  geom_col()+
  theme_test()+
  ggtitle("C.")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14),
        axis.title = element_text(size=14,face = "bold"),
        axis.text = element_text(size=10),
        plot.margin = margin(t=5,b=15,l=0,r=20),
        legend.position = "none",
        plot.title = element_text(hjust = 0,vjust = 0),
        title = element_text(size = 16,face="bold"))+
  geom_text(data = marker_table[marker_table$enrich_group == "healthy",], aes(y = feature, x = 0.1, label = feature),
            hjust = 0, size = 4) +
  geom_text(data = marker_table[marker_table$enrich_group == "hypertension",], aes(y = feature, x = -0.1, label = feature),
            hjust = 1, size = 4)+expand_limits(x=c(-6,5))+
  labs(x="LDA score",y="")+
  scale_fill_manual(values=c(hypertension="firebrick",healthy="dodgerblue4"))




pdf("output/fig1.pdf",width = 14,height = 8)
grid.arrange(p_all1,lefse_bar,ncol=2,widths=c(0.55,0.5))
dev.off()








mm_lefse <- run_lefse(
  phy,
  wilcoxon_cutoff = 1,
  group = "group",
  kw_cutoff = 1,
  multigrp_strat = TRUE,
  lda_cutoff = 0,
  taxa_rank = "all"
)
markers<-mm_lefse@marker_table|>as.matrix()|>as.data.frame()
markers$ef_lda<-as.numeric(markers$ef_lda)
markers$pvalue<-as.numeric(markers$pvalue)
markers$padj<-as.numeric(markers$padj)
write.csv(markers,"output/all_dif_species.csv",row.names = F,quote = F)
tax<-phy@tax_table|>as.data.frame()
markers$Phylum<-tax$Phylum[match(markers$feature,tax$Species)]
markers$sign<-ifelse(markers$padj<0.05&abs(markers$ef_lda)>=2,markers$enrich_group,"Non")
markers$ef_lda[markers$enrich_group=="healthy"]<- -markers$ef_lda[markers$enrich_group=="healthy"]

count_table<-as.data.frame(table(markers$sign[markers$sign!="Non"]))|>column_to_rownames("Var1")|>t()
rownames(count_table)<-""
tb_count<-gridExtra::tableGrob(count_table,theme = ttheme_default(base_size=10))

ggplot(markers,aes(ef_lda,y = -log10(padj),color=sign))+
  geom_point(size=2)+
  annotation_custom(tb_count,xmin = -2,xmax = 2,ymin = 3,ymax = 4)+
  scale_color_manual(name="",values=c(healthy="blue",Non="grey",hypertension="red"))+
  theme_cowplot()+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = c(-2,2),linetype="dashed")+
  scale_x_continuous(breaks = c(-4,-2,0,2,4,6),labels = c(4,2,0,2,4,6))
ggsave("output/summplement_volcano.pdf",height = 4,width = 6)
