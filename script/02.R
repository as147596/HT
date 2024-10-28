library(tidyverse)
library(tidyr)
library(dplyr)
library(randomForest)
library(pROC)
library(caret)
library(ggplot2)
library(vegan)
library(SHAPforxgboost)
library(xgboost)
library(tibble)
library(shapviz)
library(Rtsne)
library(glmnet)
library(ggvenn)
library(gridExtra)
library(cowplot)
library(ComplexHeatmap)


flux<-read.csv("data/discovery_cohort/exchanges.csv",header = T)
flux<-flux[,c(2,7,5)]

test_flux<-read.csv("data/test_cohort/test_16s/res_grow_g/exchanges.csv",header = T)
test_flux<-test_flux[,c(2,7,5)]

flux_all<-rbind(flux,test_flux)

flux_all<-flux_all %>%group_by(sample_id,metabolite)%>%
  summarise(flux=sum(flux))
flux_all<-flux_all%>% pivot_wider(names_from = metabolite,values_from = flux) %>% column_to_rownames("sample_id")
flux_all[is.na(flux_all)]<-0

meta_test<-read.csv("data/test_cohort/test_16s/SraRunTable.csv",header = T)

meta_healthy<-meta_test[meta_test$AGE>35&
                          meta_test$BMI>18.5&meta_test$BMI<25&
                          meta_test$Cardiometabolic_status=="Healthy"&
                          (meta_test$diastolic_bp<80&meta_test$Systolic_BP<130)&
                          meta_test$Medicament_use=="No"
                        #meta_test$Total_Cholesterol<200
                        ,]
meta_hy<-meta_test[meta_test$AGE>35&
                     #   meta_test$BMI>18.5&meta_test$BMI<25&
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

meta<-read.csv("data/meta.csv")

colnames(flux_all)<-paste("m",colnames(flux_all),sep = "_")
colnames(flux_all)<-gsub("\\[e]$","_e",colnames(flux_all))

flux<-flux_all[meta$sample_id,]
test_flux<-flux_all[meta_test$Sample.Name,]

#test_flux<-test_flux[,-grep("_m$",colnames(test_flux))]
set.seed(123)
index<-sample(1:nrow(flux),nrow(flux)*0.8,replace = F)
train<-flux[index,]
test<-flux[-index,]
train$label<-as.factor(meta$disease[index])

test$label<-as.factor(meta$disease[-index])


fit <- glmnet(as.matrix(train[,-ncol(train)]), train[,ncol(train)], family = "binomial", nlambda = 300)
# 在测试集评估模型
res=assess.glmnet(fit, newx = as.matrix(test[,-ncol(test)]), newy = test[,ncol(test)])
pred=predict(fit, newx = as.matrix(train[, -ncol(train)]), type = "response",s=fit$lambda[which.max(res$auc)])

cfit <- cv.glmnet(as.matrix(train[,-ncol(train)]), train[,ncol(train)], family = "binomial",keep = TRUE, nlambda = 300)
cf_res=assess.glmnet(cfit, newx = as.matrix(test[,-ncol(test)]), newy = test[,ncol(test)],s = "lambda.min")

rocs <- roc.glmnet(cfit$fit.preval, newy = train[,ncol(train)])
best <- cfit$index["min",] # 提取AUC最大的lambda值
plot(rocs[[best]], type = "l") # 画出AUC最大的ROC曲线
invisible(sapply(rocs, lines, col="grey")) # 把所有的ROC都画出来
lines(rocs[[best]], lwd = 2,col = "red4") # 把AUC最大的标红
text(0.8,0.2,round(cf_res$auc,4))
cof<-as.data.frame(as.matrix(coef(fit)))


set.seed(123)
rf_model <- randomForest(as.factor(label) ~ ., data=train, importance=TRUE,ntree=1000)

pr_rf<-predict(rf_model, test[,-ncol(test)],type = 'prob')
rf_roc_curve<-roc(test$label,pr_rf[,1],smooth=T)
rf_roc_curve
pr_rf2<-predict(rf_model, test_flux,type = 'prob')
rf_roc_curve2<-roc(as.factor(meta_test$group),pr_rf2[,1],smooth=T)
rf_roc_curve2
rf_importance<-importance(rf_model)

rf_importance<-importance(rf_model)
varImpPlot(rf_model,n.var = 15,main = "Top 15 flux")

top<-rownames(rf_importance)[order(rf_importance[,3],decreasing = T)][1:10]
top_flux<-rbind(train[,top],test[,top])
meta_top_flux<-c(train$label,test$label)


model_xgboost = xgboost(
  data = as.matrix(train[,-ncol(train)]),#训练集的自变量矩阵
  label = as.numeric(train$label)-1,
  max_depth = 3, 
  eta = 0.01, 
  nthread = 12, 
  nrounds = 7000,
  objective = "binary:logistic")

pred_xgb1<-predict(model_xgboost, as.matrix(test[,-ncol(test)]))
xgb_roc_curve1<-roc(test$label,pred_xgb1,smooth=T)
xgb_roc_curve1
pred_xgb2<-predict(model_xgboost, as.matrix(test_flux))
xgb_roc_curve2<-roc(meta_test$group,pred_xgb2,smooth=T)
xgb_roc_curve2

dtrain <- xgb.DMatrix(data = as.matrix(train[,-ncol(train)]),label=as.numeric(train$label)-1) #构造数据
paras <- list(eta=0.01,
              max_depth=3,
              nthread=12,
              objective="binary:logistic",
              eval_metric="auc")
xgb <- xgb.train(dtrain,params=paras,nrounds = 7000)
xgb_pred <- predict(xgb,dtrain)
xg_importance <- xgb.importance(model = xgb)
xgb.ggplot.importance(importance_matrix = xg_importance,top_n = 20)

dtest<-xgb.DMatrix(data = as.matrix(test[,-ncol(test)]),label=as.numeric(test$label)-1)
xgb_pred <- predict(xgb,dtest)
xgb_roc_curve<-roc(test$label,xgb_pred,smooth=T)
dtest<-xgb.DMatrix(data = as.matrix(test_flux),label=as.numeric(as.factor(meta_test$group))-1)
xgb_pred<-predict(xgb,dtest)
roc(as.factor(meta_test$group),xgb_pred,smooth=T)
performance<-list(`RF Validation AUC:0.87`=rf_roc_curve,`XGB Validation AUC:0.93`=xgb_roc_curve1,`RF Testing AUC:0.84`=rf_roc_curve2,`XGB Testing AUC:0.80`=xgb_roc_curve2)
p_roc<-ggroc(performance,legacy.axes=T,alpha=1,size=2)+
  theme_bw(base_size = 12)+
  scale_colour_manual(values = c("steelblue4","slateblue4","salmon4","hotpink4"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+#去掉背景网格线
  theme(legend.position = c(0.7,0.39),legend.background=element_rect(fill = "transparent"))+#图例的坐标位置
  #scale_color_tron()+
  theme(legend.text = element_text(size = 10),#图例文字
        legend.key.height = unit(0.5, "cm"),
        axis.text = element_text(size = 10))+#刻度文本
  theme(legend.title = element_blank())
rf_importance[,3]<-rf_importance[,3]/max(abs(rf_importance[,3]))
xg_importance$Importance<-xg_importance$Importance/max(abs(xg_importance$Importance))

imp_rf<-rf_importance[order(rf_importance[,3],decreasing = T)[1:20],]
imp_rf<-data.frame(Feature=rownames(imp_rf),Importance=imp_rf[,3])
imp_rf$Feature<-gsub("^m_","",imp_rf$Feature)
imp_rf$Feature<-factor(imp_rf$Feature,levels = rev(imp_rf$Feature))
imp_xg<-xg_importance[order(xg_importance$Importance,decreasing = T)[1:20],]
imp_xg$Feature<-gsub("^m_","",imp_xg$Feature)
imp_xg$Feature<-paste0(" ",imp_xg$Feature)
imp_xg$Feature<-factor(imp_xg$Feature,levels = rev(imp_xg$Feature))
imp_all<-rbind(imp_rf,imp_xg[,c(1,5)])
imp_all$class<-rep(c("RandomForest(RF)","Xgboost(XGB)"),each=20)

p_importan=ggplot(imp_all, aes(Importance, Feature))+
  geom_segment(aes(x=0, xend=Importance, y=Feature, yend=Feature), color="steelblue4", cex=1.5) +
  geom_point(aes(color=Importance),cex=3.5)+
  facet_wrap(.~class,scales = "free_y")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.3))+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size = 14),
    axis.text.x=element_text(color="black",size=12),
        axis.text.y=element_text(color="black",size=12),
        strip.text = element_text(size = 12,face="bold"))+
  scale_colour_gradient(low="#008280FF",high="#BB0021FF")+
  labs( y="",x="Scaled Importance")


#roc_imp<-grid.arrange(p_roc,p_importan,ncol=2,widths=c(4,7))
rf_importance_20<-rownames(rf_importance)[order(rf_importance[,3],decreasing = T)[1:20]]
xg_importance_20<-xg_importance$Feature[order(xg_importance$Importance,decreasing = T)[1:20]]
x<-list(RF=rf_importance_20,XGB=xg_importance_20)
x<-list_to_data_frame(x)
ven<-ggvenn(data = x,show_percentage = F,fill_color = c("steelblue4","slateblue4"),
            fill_alpha = 0.3,text_size = 5,set_name_size = 5)+expand_limits(y=2)
ven



intersection<-intersect(rf_importance_20,xg_importance_20)
write.csv(intersection,"data/xgb_rf_importance.csv",row.names = F,quote = F)
imp_rf<-rf_importance[intersection,]
imp_xg<-xg_importance[match(intersection,xg_importance$Feature),]
wight<-prop.table(c(as.numeric(rf_roc_curve2$auc),as.numeric(xgb_roc_curve2$auc)))
importance=imp_rf[,3]*wight[1]+imp_xg$Importance*wight[2]
importance_all<-data.frame(flux=rownames(imp_rf),importance=importance)
importance_all$flux<-str_sub(importance_all$flux,3,nchar(importance_all$flux))
importance_all$group<-ifelse(gsub(".*_","",importance_all$flux)=="m","medium","extracellular space")
importance_all$flux<-factor(importance_all$flux,levels = importance_all$flux[order(importance_all$importance)])
importance_all$group<-factor(importance_all$group,levels = c("medium","extracellular space"))
bar<-ggplot(data = importance_all,aes(importance,flux))+
  geom_bar(stat = "identity",fill="royalblue4", color="black",alpha=.4) +
  facet_grid(group~.,scales = "free_y",space = "free_y")+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
    axis.text.y = element_text(size=12, color="black"),
    axis.text.x = element_text(size=12, color="black")) +
 # scale_x_discrete(position = "top")+
  theme(axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14,face="bold"))+ylab("Metabolites Flux")+xlab("Importance")

bar<-bar+labs(title = "D.")+theme(title = element_text(size = 15,face = "bold"),plot.title.position = "plot",
                                  plot.title = element_text(hjust = 0,vjust = -3))
bar_ven<-ggdraw() +
  draw_plot(bar) +
  draw_plot(ven, x = 0.53, y = 0.2, width = 0.35, height = 0.35)


sparse<-readRDS("output/sparse.rds")
sparse<-grid::grid.grabExpr(draw(sparse,heatmap_legend_side = "right"))
sparse<-ggdraw(sparse)+ggtitle("A.") +
  theme(plot.title = element_text(hjust = 0.03,vjust = -4),
        title = element_text(size = 16,face="bold")
        )
p_roc1<-p_roc+theme(plot.margin = ggplot2::margin(t=15,r=10),plot.title.position = "plot",
                    plot.title = element_text(hjust = 0,vjust = 1),
                    legend.background = element_rect(fill = "transparent"),
                    title = element_text(size = 15,face="bold"))+labs(title = "B.")
#barp<-grid::rasterGrob(png::readPNG("output/bar_imp.png"), interpolate = TRUE)
p_importan<-p_importan+labs(title = "C.")+theme(title = element_text(size = 15,face = "bold"),
                                                plot.title.position = "plot", ###按照图形左对齐，默认是按照图形面板左对齐
                                                plot.title = element_text(hjust = 0.01,vjust = -4))
importance_mat<-flux[,names(importance)]
pvalues<-c()
for(i in 1:ncol(importance_mat)){
  p<-wilcox.test(as.numeric(importance_mat[,i])~meta$disease)
  pvalues<-c(pvalues,p$p.val)
}
pvalues<-ifelse(pvalues<0.0001,"***","")
importance_mat[importance_mat>1]<-1
importance_mat[importance_mat< -1]<- -1
importance_mat<-scale(importance_mat)
row_ann<-rowAnnotation(foo = anno_block(gp = gpar(fill = c("#4f6db7","orangered4")),
                                        labels = c("Healthy", "Hypertension"), 
                                        labels_gp = gpar(col = "grey95", fontsize = 13)))
colnames(importance_mat)<-gsub("^m_","",colnames(importance_mat))
colnames(importance_mat)<-paste0(colnames(importance_mat),pvalues)
importance_heatmap<-Heatmap(importance_mat, name = "scaled\nflux", cluster_rows = F,cluster_columns = T,
                    row_split = factor(rep(c(" ", "  "), c(sum(meta$disease=="hypertension"),sum(meta$disease=="healthy"))),levels = c("  "," ")),
                    left_annotation = row_ann,column_names_rot = 45,
                    show_heatmap_legend = T,
                    #col=c("grey94","lightsteelblue3"),
                    heatmap_legend_param = list(
                      title_gp = gpar(fontsize = 14),
                      labels_gp = gpar(fontsize = 14)),
                    column_title_gp = gpar(fontsize=18),
                    show_row_names = F,show_column_names = T)
importance_heatmap<-grid::grid.grabExpr(draw(importance_heatmap))
importance_heatmap<-ggdraw(importance_heatmap)+ggtitle("E.") +
  theme(plot.title = element_text(hjust = 0.03,vjust = -4),
        title = element_text(size = 16,face="bold"),axis.title.x = element_text(angle=45)
  )
tmp1<-grid.arrange(sparse, p_importan, ncol = 1,heights=c(0.5,0.85))
tmp2<-grid.arrange(p_roc1, bar_ven,importance_heatmap, ncol = 1,heights=c(0.5,0.8,0.75))
dev.off()
pdf("output/Fig2.pdf",width = 12,height = 9.5)
grid.arrange(tmp1,grid::nullGrob(),tmp2,ncol=3,widths=c(6,0.5,4))
dev.off()
