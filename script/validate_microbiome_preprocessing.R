library(curatedMetagenomicData)
hypertension_meta<-sampleMetadata[which(sampleMetadata$disease=='hypertension'),]
healthy_meta<-hypertension_meta[which(hypertension_meta$age_category=='adult'),]
hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance")
abundance_hypertension<-hypertension@assays@data@listData[["relative_abundance"]]
temp<-strsplit(rownames(abundance_hypertension),'__')

for(i in 1:nrow(abundance_hypertension)){
  rownames(abundance_hypertension)[i]<-temp[[i]][8]
}
temp<-strsplit(rownames(abundance_hypertension),'_')
for(i in 1:nrow(abundance_hypertension)){
  rownames(abundance_hypertension)[i]<-temp[[i]][length(temp[[i]])]
}
setwd('D:/microbiome/micom/agora103_species/4ce8f8cf-98ea-4f45-bd80-86bf36efb506/data')
filenames<-list.files()
for(i in 1:length(filenames)){
  filenames[i]<-substr(filenames[i],1,nchar(filenames[i])-5)
}

abundance_hypertension<-abundance_hypertension[intersect(rownames(abundance_hypertension),filenames),]
abundance_hypertension<-abundance_hypertension[,as.numeric(which(colSums(abundance_hypertension)>60))]##163

data_hypertension<-matrix(0,nrow(abundance_hypertension)*ncol(abundance_hypertension),4)
data_hypertension<-as.data.frame(data_hypertension)
colnames(data_hypertension)<-c('id','sample_id','abundance','species')

id<-rep(rownames(abundance_hypertension),ncol(abundance_hypertension))
sample_id<-c()
for(i in 1:ncol(abundance_hypertension)){
  sample_id<-c(sample_id,rep(colnames(abundance_hypertension)[i],nrow(abundance_hypertension)))
}
data_hypertension$id<-id
data_hypertension$sample_id<-sample_id
data_hypertension$abundance<-abundance_hypertension %>% as.matrix() %>% as.numeric()
data_hypertension$species<-id



data_hypertension$id<-rep(rownames(abundance_hypertension),ncol(abundance_hypertension))
data_hypertension$sample_id<-sampleids
data_hypertension$abundance<-as.numeric(abundance_hypertension)
data_hypertension$species<-rep(rownames(abundance_hypertension),ncol(abundance_hypertension))

write.csv(data_hypertension,'output/data_hypertension.csv')

##########
healthy_meta<-sampleMetadata[which(sampleMetadata$disease=='healthy'),]
healthy_meta<-healthy_meta[which(healthy_meta$age_category=='adult'),]
healthy_meta<-healthy_meta[which(healthy_meta$country=='CHN'|healthy_meta$country=='ITA'|healthy_meta$country=='AUT'),]
healthy_meta<-healthy_meta[which(healthy_meta$BMI>18.5&healthy_meta$BMI<23.9),]
healthy_meta<-healthy_meta[which(healthy_meta$age>55&healthy_meta$age<70),]

healthy<-healthy_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance")
abundance_healthy<-healthy@assays@data@listData[["relative_abundance"]]
temp<-strsplit(rownames(abundance_healthy),'__')

for(i in 1:nrow(abundance_healthy)){
  rownames(abundance_healthy)[i]<-temp[[i]][8]
}
temp<-strsplit(rownames(abundance_healthy),'_')
for(i in 1:nrow(abundance_healthy)){
  rownames(abundance_healthy)[i]<-temp[[i]][length(temp[[i]])]
}
setwd('D:/microbiome/micom/agora103_species/4ce8f8cf-98ea-4f45-bd80-86bf36efb506/data')
filenames<-list.files()
for(i in 1:length(filenames)){
  filenames[i]<-substr(filenames[i],1,nchar(filenames[i])-5)
}

abundance_healthy<-abundance_healthy[intersect(rownames(abundance_healthy),filenames),]
abundance_healthy<-abundance_healthy[,as.numeric(which(colSums(abundance_healthy)>60))]##106

data_healthy<-matrix(0,nrow(abundance_healthy)*ncol(abundance_healthy),4)
data_healthy<-as.data.frame(data_healthy)
colnames(data_healthy)<-c('id','sample_id','abundance','species')

id<-rep(rownames(abundance_healthy),ncol(abundance_healthy))
sample_id<-c()
for(i in 1:ncol(abundance_healthy)){
  sample_id<-c(sample_id,rep(colnames(abundance_healthy)[i],nrow(abundance_healthy)))
}
data_healthy$id<-id
data_healthy$sample_id<-sample_id
data_healthy$abundance<-abundance_healthy %>% as.matrix() %>% as.numeric()
data_healthy$species<-id

data<-rbind(data_healthy,data_hypertension)
write.csv(data,'output/data_all.csv')
###
setwd('D:/microbiome/micom/agora103_species/4ce8f8cf-98ea-4f45-bd80-86bf36efb506/data')
filenames<-list.files()
manifest<-matrix(0,length(filenames),3)%>%as.data.frame()
colnames(manifest)<-c('file','species','summary_rank')
manifest$file<-filenames
for(i in 1:length(filenames)){
  filenames[i]<-substr(filenames[i],1,nchar(filenames[i])-5)
}
manifest$species<-filenames
manifest$summary_rank<-'species'
write.csv(manifest,'output/manifest.csv')

