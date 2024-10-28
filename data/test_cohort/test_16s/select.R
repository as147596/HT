otu<-read.table("microbio_selected.otus",row.names = 1,header = T)
otu<-as.data.frame(t(otu))
anno<-read.table("microbio_selected.taxonomy",header = T)
anno<-anno[grep(";s__",anno$Taxonomy),]
otu<-otu[anno$OTU,]
otu<-data.frame(taxo=anno$Taxonomy,otu)
otu<-aggregate(.~taxo,data=otu,sum)

otu$taxo<-gsub(".*s__","",otu$taxo)
otu$taxo<-gsub(";$","",otu$taxo)
otu<-otu[!duplicated(otu$taxo),]
rownames(otu)<-otu$taxo
otu<-otu[,-1]
otu<-prop.table(as.matrix(otu),2)*100
write.csv(otu,"species_16s.csv")


otu<-read.table("microbio_selected.otus",row.names = 1,header = T)
otu<-as.data.frame(t(otu))
anno<-read.table("microbio_selected.taxonomy",header = T)
anno$Taxonomy<-gsub(";s__.*","",anno$Taxonomy)
anno$Taxonomy<-gsub(";unclassified;$","",anno$Taxonomy)

anno<-anno[grep(";g__",anno$Taxonomy),]
otu<-otu[anno$OTU,]
otu<-data.frame(taxo=anno$Taxonomy,otu)
otu<-aggregate(.~taxo,data=otu,sum)

otu$taxo<-gsub(".*g__","",otu$taxo)
#otu$taxo<-gsub(";$","",otu$taxo)
otu<-otu[!duplicated(otu$taxo),]
rownames(otu)<-otu$taxo
otu<-otu[,-1]
otu1<-prop.table(as.matrix(otu),2)*100
write.csv(otu1,"genus_16s.csv")
