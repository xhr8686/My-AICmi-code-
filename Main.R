# TCGA-PAAD mRNA, miRNA and clinical data were downloaded by gdc-client and from TCGA website. After pre-processing the downloaded data, we obtained exprmiRNA.Rdata as the miRNA data.

exprmiRNA<-read.csv('exprmiRNA.csv')
rownames(exprmiRNA)<-exprmiRNA[,1]
exprmiRNA<-exprmiRNA[,-1]
tumor<-exprmiRNA[,substr(colnames(exprmiRNA),14,14)=='0']
normal<-exprmiRNA[,substr(colnames(exprmiRNA),14,14)=='1']
dim(tumor)
dim(normal)
expr<-cbind(normal,tumor)
library(limma)
library(edgeR)
group = c(rep('normal',4),rep('tumor',179))
DGElist <- DGEList( counts = expr, group = group )
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(DGElist$counts)
contrast.matrix <- makeContrasts(paste0(unique(group),collapse = '-'),levels =design)
contrast.matrix[1,1] = -1	
contrast.matrix[2,1] = 1
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
fit <- lmFit(v, design)	
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)	
tempOutput <- topTable(fit2, coef =1, n = Inf)
nrDEG = na.omit(tempOutput)
save(nrDEG,file='nrDEG.Rdata')
nrDEGsig<-nrDEG[abs(nrDEG[,1])>1&nrDEG[,4]<0.05,]
write.csv(nrDEGsig,file='nrDEGsig.csv')

#We obtained 69 significantly differentially expressed miRNAs between pancreatic tumor tissues and adjacent normal tissues. Then we draw a volcano plot based on the nrDEG file. 

dataset<-nrDEG
cut_off_pvalue = 0.05
cut_off_logFC = 1
dataset$change = ifelse(dataset$P.Value < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, 
                     ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                     'Stable')
ggplot(
  dataset, 
  aes(x = logFC, 
      y = -log10(P.Value), 
      colour=change)) +
      geom_point(alpha=0.4, size=3.5) +
      scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+

  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +

  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+

  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
)

# We used 69 DEmiRNAs and the clinical data to build a LASSO-based model. 
library('DESeq2')
library('survival')
library('survminer')
library('dplyr')
library('glmnet')
library('ggplot2')
library('GGally')
library('rms')
library('survivalROC')
library('plotROC')
survival<-read.csv('survival.csv')
exprSet<-expr
group_list<-as.factor(c(rep('tumor',ncol(tumor)),rep('normal',ncol(normal))))
colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                                                      colData = colData,#Rows of colData correspond to columns of countData
                                                      design = ~ group_list)
dds2 <- DESeq(dds)

dds_DE<-dds2
rld <- varianceStabilizingTransformation(dds_DE, blind = T)
expr_norm <- assay(rld)
DESeq_norm_vst_for_survival <- expr_norm[rownames(nrDEGsig),]
data<- DESeq_norm_vst_for_survival[,substr(colnames(DESeq_norm_vst_for_survival),14,14)=='0']
data1<-data
colnames(data1)<-substr(colnames(data1),1,12)
data1<-t(data1)

duplicated_rows <- duplicated(rownames(data1)) | duplicated(rownames(data1), fromLast = TRUE)
data1_avg <- aggregate(data1, by = list(rownames(data1)), FUN = mean)
rownames(data1_avg) <- data1_avg$Group.1
data1_avg <- data1_avg[, -1]  
data1_avg <-data1_avg[rownames(data1),]
data1<-data1_avg

result<-survival
result1<-result
result1<-result1[!duplicated(result1[,4]),]
rownames(result1)<-result1[,4]
rownames(result1)<-gsub('-','.',rownames(result1))
result1<-result1[,-1]
result1<-result1[,-3]
result1[,2]<-result1[,2]+1
com<-intersect(rownames(data1),rownames(result1))
result1<-result1[com,]
data1<-data1[com,]
colnames(result1)[1]<-'time'
colnames(result1)[2]<-'status'

x<-as.matrix(data1)
y=as.double(result1$time)
status<-as.double(result1$status)
surv=Surv(y,status)
cv.fit<-cv.glmnet(x,surv,family="cox",)
plot(cv.fit)
fit<-glmnet(x,surv,family="cox",alpha=0.3)
coefficients<-coef(fit,s=cv.fit$lambda.min)
coefficients
coef <-coefficients
nz  <- coef[, 1] != 0
Dataset2<-data.frame(miRNA   = rownames(coef)[nz],coefficient = coef[nz, 1])
write.csv(Dataset2,file='Dataset2.csv')
newcombine<-cbind(result1,data1)
newcombine[,2]<-newcombine[,2]-1

# The 69 DEmiRNAs and clinical data were input for an XGBoost training. 
library(xgboost)
library(survival)
library(survminer)
#set.seed(123)
colnames(newcombine)[1]<-'time'
colnames(newcombine)[2]<-'status'
trainData <-newcombine
train.x <- as.matrix(trainData[, colnames(newcombine)[3:ncol(newcombine)]])
train.y <- ifelse(trainData$status == 1, trainData$time, -trainData$time)
trainMat <- xgb.DMatrix(data = train.x, label = train.y)

param <- list(objective = "survival:cox",
              booster = "gbtree",
              eval_metric = "cox-nloglik",
              eta = 0.03,
              max_depth = 3,
              subsample = 1,
              colsample_bytree = 1,
              gamma = 0.5)

#set.seed(1)

xgb.fit <- xgb.train(params = param, data = trainMat, nrounds = 1000, watchlist = list(val2 = trainMat),
    early_stopping_rounds = 50)
riskScore <- predict(xgb.fit, newdata = train.x)
hist(riskScore)

groups <- ifelse(riskScore>median(riskScore),"high","low")
f <- survfit(Surv(trainData$time, trainData$status) ~ groups)
ggsurvplot(f,
           data = trainData,
           surv.median.line = "hv", 
           #legend.title = "Risk Group",
           #legend.labs = c("Low Risk", "High Risk"),
           pval = TRUE, 
           ggtheme = theme_bw()
)

impMatrix <- xgb.importance(feature_names = dimnames(train.x)[[2]], model = xgb.fit)
impMatrix
xgb.plot.importance(impMatrix, main = "Gain by Feature")

# The 69 DEmiRNAs and clinical data were input for a GBM training.
library("gbm")
train<-newcombine
set.seed(201)
gbm_model <- gbm(Surv(time, status) ~ .,
             distribution = "coxph",
             data = train,
             n.trees = 1000,
             shrinkage = 0.1,
             interaction.depth = 5,
             n.minobsinnode = 5,
             cv.folds = 10
)
summary(gbm_model)

best.iter <- gbm.perf(gbm_model, plot.it = TRUE, method = "cv")
pred_train <- predict(gbm_model, train, n.trees = best.iter)

groups <- ifelse(pred_train >0,"high","low")
f <- survfit(Surv(train$time, train$status) ~ groups)

ggsurvplot(f,
           data = train,
           surv.median.line = "hv", 
           #legend.title = "Risk Group",
           #legend.labs = c("Low Risk", "High Risk"),
           pval = TRUE, 
           ggtheme = theme_bw()
)

#K-means Clustering
mydata<-newcombine[,c('hsa-mir-203a','hsa-mir-4787')]

library(factoextra)
library(cluster)
library(ggplot2)
library(cluster)

mydata_scaled <- scale(mydata)
head(mydata_scaled)
wss <- sapply(1:10, function(k){kmeans(mydata_scaled, k, nstart=25)$tot.withinss})

fviz_nbclust(mydata_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2) +
  labs(title = "K means")
sil_width <- sapply(2:8, function(k){
  model <- kmeans(mydata_scaled, centers = k, nstart = 25)
  ss <- silhouette(model$cluster, dist(mydata_scaled))
  mean(ss[, 3])
})

plot(2:8, sil_width, type = "b", xlab = "Number of clusters", 
     ylab = "Average Silhouette Width", main = "K means")

set.seed(123)
kmeans_result <- kmeans(mydata_scaled, centers = 3)
cluster_labels <- kmeans_result$cluster  

library(pheatmap)
data_with_cluster <- cbind(mydata_scaled, cluster = factor(cluster_labels))
data_with_cluster<-as.data.frame(data_with_cluster)
ordered_data <- data_with_cluster[order(data_with_cluster$cluster), ]
a<-data.frame(Cluster = ordered_data$cluster)
rownames(a)<-rownames(ordered_data[, -ncol(ordered_data)])

pheatmap(
    ordered_data[, -ncol(ordered_data)], 
    cluster_rows = FALSE,           
    cluster_cols = TRUE,             
    annotation_row = a, 
    show_rownames = FALSE,        
    main = "K-means Clustering Heatmap"
)




