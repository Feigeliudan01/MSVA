################### R code for MSVA: A Novel Approach to Removing Unwanted Variation###################
##################### in Metabolomics Data Based on Surrogate Variable Analysis #######################

########################### PART1: function used in this codes #######################
### Group RSDs into 7 groups
RSD_group<-function(data){
  if (class(data)!="numeric"){
    stop("data must be numeric")
  }
  group<-vector(mode="numeric",length=length(data))
  for (i in (1:length(data))){
    if (0<data[i] & data[i]<=0.15){
      group[i]<-1
    } else if (data[i]<=0.30){
      group[i]<-2
    } else if (data[i]<=0.45){
      group[i]<-3
    } else if (data[i]<=0.60){
      group[i]<-4
    } else if (data[i]<=0.75){
      group[i]<-5
    } else if (data[i]<=0.9){
      group[i]<-6
    } else (group[i]<-7)
  }
  return(group)
}

## group RSDs into 3 groups
RSD_group3<-function(data){
  if (class(data)!="numeric"){
    stop("data must be numeric")
  }
  group<-vector(mode="numeric",length=length(data))
  for (i in (1:length(data))){
    if (0<data[i] & data[i]<=0.15){
      group[i]<-3
    } else if (data[i]<=0.30){
      group[i]<-2
    } else  (group[i]<-1)
  }
  return(group)
}

### group colors
color_group<-function(data,ratio){
  if (class(data)!="numeric"){
    stop("data must be numeric")
  }
  group<-vector(mode="character",length=length(data))
  for (i in (1:length(data))){
    if (0<data[i] & data[i]<=ratio){
      group[i]<-"red"
    } else (group[i]<-"grey")
  }
  return(group)
}





###################### PART 2: Removal of unwanted intra-batch variation #################
### Import data 
setwd("/home/stat-dengkui/intrabatch")
data_T3_bad<-read.csv("intra_batch_data.csv",header=T,row.names = 1,stringsAsFactors = F)
data_ZF_order<-read.csv("sample.information.csv",header=T,row.names = 1,stringsAsFactors = F)
id<-read.csv("id.csv",header=T,stringsAsFactors = F)
data_T3_bad<-data_T3_bad[,-c(1:11)]
data_T3_bad_tran<-t(data_T3_bad)

count_zero1<-apply(data_T3_bad_tran,2,function(x) any(x==0))
data_T3_bad_tran<-data_T3_bad_tran[,!count_zero1]

data_T3_bad<-data_T3_bad_tran
data_T3_bad<-data_T3_bad

#### Match sample information
# Change the row name of the data
location<-str_locate_all(rownames(data_T3_bad),"_")
locate_lower<-sapply(location,function(x) x[1,1])
locate_upper<-sapply(location,function(x) x[2,1])

rownames(data_T3_bad)<-str_sub(rownames(data_T3_bad),locate_lower+1,locate_upper-1)
data_T3_bad<-as.data.frame(data_T3_bad,stringsAsFactors = F)
# Change the row names of sample information
data_ZF_order_QC<-data_ZF_order[substr(rownames(data_ZF_order),1,2)=="QC",]
data_ZF_order_sample<-data_ZF_order[substr(rownames(data_ZF_order),1,2)!="QC",]

rownames(data_ZF_order_sample)<-str_sub(rownames(data_ZF_order_sample),str_locate(rownames(data_ZF_order_sample),"_")[,1]+1)
data_ZF_order<-rbind(data_ZF_order_QC,data_ZF_order_sample)

## merge with sample information
data_ZF_order$name<-rownames(data_ZF_order)
data_T3_bad$name<-rownames(data_T3_bad)
data_bad_merge<-merge(data_ZF_order,data_T3_bad,by="name")

data_bad<-data_bad_merge[,-c(1:5)]

data_bad<-data_bad[-138,]
data_bad_merge<-data_bad_merge[-138,]

#### PCA score plot for original data
label<-vector(mode="numeric",length=nrow(data_bad))
label[data_bad_merge$group==0]<-1
label[data_bad_merge$group==1]<-2
label[data_bad_merge$group=="QC"]<-3


bitmap(file="PCA_original.jpeg",type="jpeg",res=800)
pca(data_bad,id=label,plot=1)
dev.off()

#################### MSVA method####################
library(sva)
set.seed(1234)

mod = model.matrix(~as.factor(label), data=data_bad)
mod0 = model.matrix(~1,data=data_bad)
svobj_bad<-sva(t(data_bad),mod=mod,mod0=mod0)


data_bad_model<-cbind(data_bad,svobj_bad$sv,label=as.matrix(label))
colnames(data_bad_model)[6113:6123]<-c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","label")

#### MSVA for QC samples and subject samples
data_bad_model_QC<-data_bad_model[label==3,]
data_bad_model_sample<-data_bad_model[label!=3,]

### MSVA for QC samples
data_bad_sva_QC<-matrix(NA,17,6112)
for (i in (1:6112)){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_bad_model)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_bad_model_QC))
  data_bad_sva_QC[,i]<-data_bad_model_QC[,i]-predict(model_lm)+mean(data_bad_model_QC[,i],na.rm=T)
}

### MSVA for subject samples
data_bad_sva_sample<-matrix(NA,136,6112)
for (i in (1:6112)){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_bad_model)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_bad_model_QC))
  data_bad_sva_sample[,i]<-data_bad_model_sample[,i]-predict(model_lm,newdata=data_bad_model_sample[,c(i,6113:6123)])+mean(data_bad_model_sample[,i],na.rm=T)
}

data_bad_sva<-rbind(data_bad_sva_sample,data_bad_sva_QC)

### PCA score plot 
label_bad<-c(label[label!=3],label[label==3])
bitmap(file="PCA_msva.jpeg",type="jpeg",res=800)
pca(data_bad_sva,id=label_bad,plot=1)
dev.off()


### Pearson coefficient of QC samples
data_bad_QC<-data_bad[data_bad_merge$group=="QC",]
data_correlation_bad_sva<-correlation(data_bad_QC,data_bad_sva_QC)
mean(data_correlation_bad_sva[,3])
mean(data_correlation_bad_sva[,4])


### RSD of QCs
RSD_original_bad<-apply(data_bad_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
RSD_bad_sva<-apply(data_bad_sva_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))


RSD_original_bad_0.15<-sum(RSD_original_bad<0.15)/length(RSD_original_bad)
RSD_bad_sva_0.15<-sum(RSD_bad_sva<0.15)/length(RSD_bad_sva)

RSD_original_bad_0.30<-sum(RSD_original_bad<0.30)/length(RSD_original_bad)
RSD_sva_bad_0.30<-sum(RSD_bad_sva<0.30)/length(RSD_bad_sva)


##################### quantile normalization ####################
library(preprocessCore)
data_bad_quantile<-normalize.quantiles(t(data_bad),copy=F)

#PCA plot
bitmap(file="PCA_quantile.jpeg",res=800,type="jpeg")
pca(t(data_bad_quantile),id=label,plot=1)
dev.off()

#correlation of QC
data_bad_quantile_QC<-t(data_bad_quantile)[label==3,]

data_correlation_bad_quantile<-correlation(data_bad_QC,data_bad_quantile_QC)
mean(data_correlation_bad_quantile[,3])
mean(data_correlation_bad_quantile[,4])

# RSDs of QCs
RSD_bad_quantile<-apply(data_bad_quantile_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
RSD_bad_quantile_0.15<-sum(RSD_bad_quantile<0.15)/length(RSD_bad_quantile)
RSD_bad_quantile_0.30<-sum(RSD_bad_quantile<0.30)/length(RSD_bad_quantile)

############# Regression Calibration Method  #############
# Seperate the original data into QC and subject samples
data_bad_QC<-data_bad[label==3,]
data_bad_merge_QC<-data_bad_merge[label==3,]
data_bad_QC_R<-cbind(data_bad_QC,order=data_bad_merge_QC$injection.order)

data_bad_sample<-data_bad[label!=3,]
data_bad_merge_sample<-data_bad_merge[label!=3,]
data_bad_sample_R<-cbind(data_bad_sample,order=data_bad_merge_sample$injection.order)

###for QC samples ####
data_bad_regression_QC<-matrix(NA,17,6112)
for (i in (1:(ncol(data_bad_QC_R)-1))){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_bad_QC_R)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"order",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_bad_QC_R))
  data_bad_regression_QC[,i]<-data_bad_QC_R[,i]-predict(model_lm)+mean(data_bad_QC_R[,i],na.rm=T)
  
}

###For subject samples  #####
data_bad_regression_sample<-matrix(NA,136,6112)
for (i in (1:(ncol(data_bad_sample_R)-1))){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_bad_sample_R)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"order",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_bad_sample_R))
  data_bad_regression_sample[,i]<-data_bad_sample_R[,i]-predict(model_lm)+mean(data_bad_sample_R[,i],na.rm=T)
  
}


### PCA score plot 
data_bad_regression<-rbind(data_bad_regression_sample,data_bad_regression_QC)
label_R<-c(label[label!=3],label[label==3])
bitmap(file="PCA_regression.jpeg",type="jpeg",res=800)
pca(data_bad_regression,id=label_R,plot=1)
dev.off()
### pearson coefficients of QCs
data_bad_regression_correlation<-correlation(data_bad_QC_R[,-ncol(data_bad_QC_R)],data_bad_regression_QC)
mean(data_bad_regression_correlation[,3])
mean(data_bad_regression_correlation[,4])

## RSDs of QCs
RSD_bad_regression<-apply(data_bad_regression_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))

RSD_bad_regression_0.15<-sum(RSD_bad_regression<0.15)/length(RSD_bad_regression)

RSD_bad_regression_0.30<-sum(RSD_bad_regression<0.30)/length(RSD_bad_regression)



###################### median scaling #############
data_bad_median<-apply(data_bad,1,function(x) x/median(x,na.rm=T))

# PCA score plot
bitmap(file="PCA_median.jpeg",res=800,type="jpeg")
pca(t(data_bad_median),id=label,plot=1)
dev.off()

# Pearson correlation coefficients of QCs
data_bad_median_QC<-t(data_bad_median)[label==3,]

data_bad_correlation_median<-correlation(data_bad_QC,data_bad_median_QC)
mean(data_bad_correlation_median[,3])
mean(data_bad_correlation_median[,4])

# RSDs of QCs
RSD_bad_median<-apply(data_bad_median_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
RSD_bad_median_0.15<-sum(RSD_bad_median<0.15)/length(RSD_bad_median)
RSD_bad_median_0.30<-sum(RSD_bad_median<0.30)/length(RSD_bad_median)


######### Median Fold Change ############
## For QC samples
data_bad_QC<-as.matrix(data_bad_QC)
median_var_QC<-apply(data_bad_QC,2,function(x) x/median(x,na.rm=T))
median_sample_QC<-apply(median_var_QC,1,median,na.rm=T)

data_bad_FC_QC<-matrix(0,17,6112)

for (i in (1:nrow(data_bad_QC))){
  cat(paste("########",i,"###########\n"))
  data_bad_FC_QC[i,]<-data_bad_QC[i,]/median_sample_QC[i]
}

## For subject samples
data_bad_sample<-as.matrix(data_bad_sample)
median_var_sample<-apply(data_bad_sample,2,function(x) x/median(x,na.rm=T))
median_sample_sample<-apply(median_var_sample,1,median,na.rm=T)

data_bad_FC_sample<-matrix(0,136,6112)

for (i in (1:nrow(data_bad_sample))){
  cat(paste("########",i,"###########\n"))
  data_bad_FC_sample[i,]<-data_bad_sample[i,]/median_sample_sample[i]
}

### PCA score plot
data_bad_FC<-rbind(data_bad_FC_sample,data_bad_FC_QC)
label_fc<-c(label[label!=3],label[label==3])

bitmap(file="PCA_FC.jpeg",type="jpeg",res=800)
pca(data_bad_FC,id=label_fc,plot=1)
dev.off()

# Pearson correlation coefficients of QCs
data_bad_correlation_FC<-correlation(data_bad_QC,data_bad_FC_QC)
mean(data_bad_correlation_FC[,3])
mean(data_bad_correlation_FC[,4])

## RSDs of QCs
RSD_bad_fc<-apply(data_bad_FC_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))

RSD_bad_FC_0.15<-sum(RSD_bad_fc<0.15)/length(RSD_bad_fc)
RSD_bad_FC_0.30<-sum(RSD_bad_fc<0.30)/length(RSD_bad_fc)



###################### Figures ##########################
############## RSD ##########
###### RSD barplot

##original
# group RSDs
#Original data
library(ggplot2)
data<-RSD_original_bad
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x1<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,4000),breaks=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Raw")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

## MSVA data
data<-RSD_bad_sva
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")

x2<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,4000),breaks=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="MSVA")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

##quantile
data<-RSD_bad_quantile
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")

x3<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,4000),breaks=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Quantile")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

#### Regression
data<-RSD_bad_regression
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x4<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,4000),breaks=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Regression")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))


#### Median scaling
data<-RSD_bad_median
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x5<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,4000),breaks=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Median scaling")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

###median fold change

data<-RSD_bad_fc
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")

x6<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,4000),breaks=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Median fold change")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

#### put them together
bitmap(file="RSD_barplot_zong_label.jpeg",type="jpeg",res=800,height=7,width=12)
library(grid)
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))
print(x1,vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
print(x2,vp=viewport(layout.pos.row = 1,layout.pos.col = 2))
print(x5,vp=viewport(layout.pos.row = 1,layout.pos.col = 3))

print(x6,vp=viewport(layout.pos.row = 2,layout.pos.col = 1))
print(x3,vp=viewport(layout.pos.row = 2,layout.pos.col = 2))
print(x4,vp=viewport(layout.pos.row = 2,layout.pos.col = 3))

dev.off()
############################################
##### RSD percentchart ##########
data<-RSD_original_bad
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Raw",nrow(RSD_data)))
RSD_data1<-RSD_data
rownames(RSD_data1)<-names(RSD_original_bad)

data<-RSD_bad_sva
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("MSVA",nrow(RSD_data)))
RSD_data2<-RSD_data
rownames(RSD_data2)<-names(RSD_bad_sva)




data<-RSD_bad_quantile
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Quantile",nrow(RSD_data)))
RSD_data3<-RSD_data
rownames(RSD_data3)<-names(RSD_bad_quantile)


data<-RSD_bad_regression
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Regression",nrow(RSD_data)))
RSD_data4<-RSD_data
rownames(RSD_data4)<-names(RSD_bad_regression)

data<-RSD_bad_median
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Median scaling",nrow(RSD_data)))
RSD_data5<-RSD_data
rownames(RSD_data5)<-names(RSD_bad_median)

data<-RSD_bad_fc
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Median fold change",nrow(RSD_data)))
RSD_data6<-RSD_data
rownames(RSD_data6)<-names(RSD_bad_fc)

data_integral<-rbind(RSD_data1,RSD_data2,RSD_data3,RSD_data4,RSD_data5,RSD_data6)

bitmap(file="RSD percent barPlot_label.jpeg",type="jpeg",res=800,height=8,width=9)
ggplot(data=data_integral,aes(x=group,fill=factor(x)))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D"),name="RSD",guide=guide_legend(reverse=T),label=c(">30%","15-30%","0-15%"))+
  theme(axis.text=element_text(size=rel(1.2),face="bold"),axis.title=element_text(size=rel(1.5),face="bold"))+
  theme(legend.background   = element_rect(fill="gray90"),legend.title=element_text(size=rel(1.2),face="bold"),legend.text =element_text(size=rel(1.2),face="bold"))+
  
  scale_x_discrete(limits=c("Raw","MSVA","Median scaling","Median fold change","Quantile","Regression"),labels=c("Raw","MSVA","Median\n scaling","Median\n fold change","Quantile","Regression"))+
  scale_y_continuous(breaks=c(0,0.25,0.50,0.75,1),labels=c(0,25,50,75,100),name="% of Peaks")+
  labs(x="")


dev.off()


####PCA score plot

###Original data
library(ggfortify)
data_auto_PCA_original<-cbind(group=label,data_bad)
data_auto_PCA_original$group<-as.factor(data_auto_PCA_original$group)
PC_original<-prcomp(scale(data_bad),retx=T)

bitmap(file="PCA.original.jpeg",type="jpeg",res=1000)
autoplot(PC_original,data=data_auto_PCA_original,colour="group")+scale_colour_discrete(name="",label=c("CE","CRC","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(2.5,"cm"),legend.key.height = unit(1,"cm"))+
  labs(title="Raw")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )
dev.off()

s1<-autoplot(PC_original,data=data_auto_PCA_original,colour="group")+scale_colour_discrete(name="",label=c("CE","CRC","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Raw")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )

### MSVA
library(ggfortify)
data_auto_PCA_sva<-cbind(group=label_bad,data_bad_sva)
data_auto_PCA_sva<-as.data.frame(data_auto_PCA_sva)
data_auto_PCA_sva$group<-as.factor(data_auto_PCA_sva$group)
PC_sva<-prcomp(scale(data_bad_sva),retx=T)

s2<-autoplot(PC_sva,data=data_auto_PCA_sva,colour="group")+scale_colour_discrete(name="",label=c("CE","CRC","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="MSVA")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )

###Quantile 
library(ggfortify)
data_auto_PCA_quantile<-cbind(group=label,t(data_bad_quantile))
data_auto_PCA_quantile<-as.data.frame(data_auto_PCA_quantile)
data_auto_PCA_quantile$group<-as.factor(data_auto_PCA_quantile$group)
PC_quantile<-prcomp(scale(t(data_bad_quantile)),retx=T)


s3<-autoplot(PC_quantile,data=data_auto_PCA_quantile,colour="group")+scale_colour_discrete(name="",label=c("CE","CRC","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Quantile")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )


###Regression

library(ggfortify)
data_auto_PCA_regression<-cbind(group=label_R,data_bad_regression)
data_auto_PCA_regression<-as.data.frame(data_auto_PCA_regression)
data_auto_PCA_regression$group<-as.factor(data_auto_PCA_regression$group)
PC_regression<-prcomp(scale(data_bad_regression),retx=T)


s4<-autoplot(PC_regression,data=data_auto_PCA_regression,colour="group")+scale_colour_discrete(name="",label=c("CE","CRC","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Regression")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )


###Median scaling

library(ggfortify)
data_auto_PCA_median<-cbind(group=label,t(data_bad_median))
data_auto_PCA_median<-as.data.frame(data_auto_PCA_median)
data_auto_PCA_median$group<-as.factor(data_auto_PCA_median$group)
PC_median<-prcomp(scale(t(data_bad_median)),retx=T)


s5<-autoplot(PC_median,data=data_auto_PCA_median,colour="group")+scale_colour_discrete(name="",label=c("CE","CRC","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Median scaling")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.6) )

###Median fold change

library(ggfortify)
data_auto_PCA_FC<-cbind(group=label_fc,data_bad_FC)
data_auto_PCA_FC<-as.data.frame(data_auto_PCA_FC)
data_auto_PCA_FC$group<-as.factor(data_auto_PCA_FC$group)
PC_FC<-prcomp(scale(data_bad_FC),retx=T)


s6<-autoplot(PC_FC,data=data_auto_PCA_FC,colour="group")+scale_colour_discrete(name="",label=c("CE","CRC","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Median fold change")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )

########## put them together

bitmap(file="PCA_zong.jpeg",type="jpeg",res=800,height=7,width=12)
library(grid)
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=3)))
print(s1,vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
print(s2,vp=viewport(layout.pos.row = 1,layout.pos.col = 2))
print(s5,vp=viewport(layout.pos.row = 1,layout.pos.col = 3))

print(s6,vp=viewport(layout.pos.row = 2,layout.pos.col = 1))
print(s3,vp=viewport(layout.pos.row = 2,layout.pos.col = 2))
print(s4,vp=viewport(layout.pos.row = 2,layout.pos.col = 3))

dev.off()

####################### Pearson correlation coefficients
################# heatmaps
###original
data_bad_merge_QC<-data_bad_merge[str_sub(data_bad_merge$name,1,2)=="QC",]
data_bad_merge_QC$name[nchar(data_bad_merge_QC$name)==3]<-paste("QC",0,str_sub(data_bad_merge_QC$name[nchar(data_bad_merge_QC$name)==3],3),sep="")
data_bad_QC<-data_bad_QC[order(data_bad_merge_QC$name),]

rownames(data_bad_QC)<-paste("QC",1:17,sep="")
cor_original<-cor(t(data_bad_QC))
bitmap(file="heatmap_original.jpeg",type="jpeg",res=1000)

library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
autoplot(cor_original,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.7,1),name="r value")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:17,sep=""))+
  scale_y_discrete(limits=paste("QC",1:17,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  labs(x="",y="")+
  guides(fill=F)
dev.off()




### MSVA 
data_bad_sva_QC<-data_bad_sva_QC[order(data_bad_merge_QC$name),]
rownames(data_bad_sva_QC)<-paste("QC",1:17,sep="")
cor_sva<-cor(t(data_bad_sva_QC))
bitmap(file="heatmap_msva.jpeg",type="jpeg",res=1000)
autoplot(cor_sva,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.7,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:17,sep=""))+
  scale_y_discrete(limits=paste("QC",1:17,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)
dev.off()


######Quantile
data_bad_quantile_QC<-data_bad_quantile_QC[order(data_bad_merge_QC$name),]

rownames(data_bad_quantile_QC)<-paste("QC",1:17,sep="")
cor_quantile<-cor(t(data_bad_quantile_QC))
bitmap(file="heatmap_quantile.jpeg",type="jpeg",res=1000)
autoplot(cor_quantile,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.7,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:17,sep=""))+
  scale_y_discrete(limits=paste("QC",1:17,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)

dev.off()


#### Regression
data_bad_regression_QC<-data_bad_regression_QC[order(data_bad_merge_QC$name),]


rownames(data_bad_regression_QC)<-paste("QC",1:17,sep="")
cor_regression<-cor(t(data_bad_regression_QC))
bitmap(file="heatmap_regression1.jpeg",type="jpeg",res=1000)
autoplot(cor_regression,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.7,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:17,sep=""))+
  scale_y_discrete(limits=paste("QC",1:17,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)

dev.off()


#### Median
data_bad_median_QC<-data_bad_median_QC[order(data_bad_merge_QC$name),]

rownames(data_bad_median_QC)<-paste("QC",1:17,sep="")
cor_median<-cor(t(data_bad_median_QC))
bitmap(file="heatmap_median.jpeg",type="jpeg",res=1000)
autoplot(cor_median,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.7,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:17,sep=""))+
  scale_y_discrete(limits=paste("QC",1:17,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)

dev.off()

#### Median fold change
data_bad_FC_QC<-data_bad_FC_QC[order(data_bad_merge_QC$name),]

rownames(data_bad_FC_QC)<-paste("QC",1:17,sep="")
cor_FC<-cor(t(data_bad_FC_QC))
bitmap(file="heatmap_FC.jpeg",type="jpeg",res=1000)
autoplot(cor_FC,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.7,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:17,sep=""))+
  scale_y_discrete(limits=paste("QC",1:17,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)

dev.off()


######## The Pearson correlation coefficient for selected 6 QCs
##### original
data_bad_QC_tran<-t(data_bad_QC[c(2,7,9,10,11,12),])
original_cor<-cor(data_bad_QC_tran)

bitmap(file="scatterPlot_Original.jpeg",type="jpeg",res=1000)
pairs(~QC2+QC7+QC9+QC10+QC11+QC12,data=log(data_bad_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Original.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

#### MSVA
data_bad_QC_tran<-t(data_bad_sva_QC[c(2,7,9,10,11,12),])
original_cor<-cor(data_bad_QC_tran)

bitmap(file="scatterPlot_MSVA.jpeg",type="jpeg",res=1000)
pairs(~QC2+QC7+QC9+QC10+QC11+QC12,data=log(data_bad_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_MSVA.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

####Quantile
data_bad_QC_tran<-t(data_bad_quantile_QC[c(2,7,9,10,11,12),])
original_cor<-cor(data_bad_QC_tran)

bitmap(file="scatterPlot_Quantile.jpeg",type="jpeg",res=1000)
pairs(~QC2+QC7+QC9+QC10+QC11+QC12,data=log(data_bad_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Quantile.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()


###Regression
data_bad_QC_tran<-t(data_bad_regression_QC[c(2,7,9,10,11,12),])
original_cor<-cor(data_bad_QC_tran)

bitmap(file="scatterPlot_Regression.jpeg",type="jpeg",res=1000)
pairs(~QC2+QC7+QC9+QC10+QC11+QC12,data=log(data_bad_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Regression.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

###median
data_bad_QC_tran<-t(data_bad_median_QC[c(2,7,9,10,11,12),])
original_cor<-cor(data_bad_QC_tran)

bitmap(file="scatterPlot_Median.jpeg",type="jpeg",res=1000)
pairs(~QC2+QC7+QC9+QC10+QC11+QC12,data=log(data_bad_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Median.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2)
dev.off()

###Median fold change
data_bad_QC_tran<-t(data_bad_FC_QC[c(2,7,9,10,11,12),])
original_cor<-cor(data_bad_QC_tran)

bitmap(file="scatterPlot_FC.jpeg",type="jpeg",res=1000)
pairs(~QC2+QC7+QC9+QC10+QC11+QC12,data=log(data_bad_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_FC.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

########################## PART3: Removal of unwanted inter-batch variation #################
#### Import data
setwd("/home/stat-dengkui/interbatch")
data_xinxueguan<-read.csv("inter_batch_data.csv",header=T,stringsAsFactors = F,row.names = 1)
sample_infor<-read.csv("sample.information.csv",header=T,stringsAsFactors = F)
data_xinxueguan_tran<-t(data_xinxueguan)
data_xinxueguan_tran<-as.data.frame(data_xinxueguan_tran)
data_xinxueguan_tran$sample<-rownames(data_xinxueguan_tran)
data_xinxueguan_merge<-base::merge(data_xinxueguan_tran,sample_infor,by="sample")
rownames(data_xinxueguan_merge)<-data_xinxueguan_merge[,1]
data_xinxueguan_merge<-data_xinxueguan_merge[,-1]
data_xinxueguan_new<-data_xinxueguan_merge[,-which(colnames(data_xinxueguan_merge)=="order")]
data_xinxueguan_new<-data_xinxueguan_new[order(data_xinxueguan_merge$order),]

batch<-vector(length=162,mode="numeric")
batch[substr(rownames(data_xinxueguan_new),1,2)=="qc"]<-c(rep(1,25),rep(2,15))
batch[substr(rownames(data_xinxueguan_new),1,2)!="qc"]<-c(rep(1,97),rep(2,25))

group<-vector(length=162,mode="numeric")
group[substr(rownames(data_xinxueguan_new),1,2)=="qc"]<-3
group[substr(rownames(data_xinxueguan_new),1,7)=="disease"]<-2
group[substr(rownames(data_xinxueguan_new),1,7)=="control"]<-1

data_xinxueguan_QC<-data_xinxueguan_new[group==3,]
data_xinxueguan_sample<-data_xinxueguan_new[group!=3,]


library(sva)
set.seed(1234)

mod = model.matrix(~as.factor(group)+as.factor(batch), data=data_xinxueguan_new)
mod0 = model.matrix(~as.factor(batch),data=data_xinxueguan_new)
svobj<-sva(t(data_xinxueguan_new),mod=mod,mod0=mod0)

data_xinxueguan_model<-cbind(data_xinxueguan_new,svobj$sv,label=as.matrix(group),batch=as.matrix(batch))
colnames(data_xinxueguan_model)[2317:2319]<-c("PC1","label","batch")

#### MSVA for QC samples and subject samples
data_xinxueguan_model_QC<-data_xinxueguan_model[group==3,]
data_xinxueguan_model_sample<-data_xinxueguan_model[group!=3,]


### MSVA of QCs
data_xinxueguan_sva_QC<-matrix(NA,40,2316)
for (i in (1:2316)){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_xinxueguan_model)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"PC1+batch",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_xinxueguan_model_QC))
  data_xinxueguan_sva_QC[,i]<-data_xinxueguan_model_QC[,i]-predict(model_lm)+mean(data_xinxueguan_model_QC[,i],na.rm=T)
}

### MSVA of subject samples
data_xinxueguan_sva_sample<-matrix(NA,122,2316)
for (i in (1:2316)){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_xinxueguan_model)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"PC1+batch",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_xinxueguan_model_QC))
  data_xinxueguan_sva_sample[,i]<-data_xinxueguan_model_sample[,i]-predict(model_lm,newdata=data_xinxueguan_model_sample[,c(i,2317:2319)])+mean(data_xinxueguan_model_sample[,i],na.rm=T)
}


### PCA score plot 
data_xinxueguan_sva<-rbind(data_xinxueguan_sva_sample,data_xinxueguan_sva_QC)
label<-c(group[group!=3],group[group==3])
bitmap(file="PCA_msva.jpeg",type="jpeg",res=800)
pca(data=data_xinxueguan_sva,id=label,plot=1)
dev.off()

bitmap(file="PCA_original.jpeg",type="jpeg",res=800)
pca(data_xinxueguan_new,id=group,plot=1)
dev.off()
### Pearson correlation coefficient for QC sampls
data_xinxueguan_QC<-data_xinxueguan_new[group==3,]
data_correlation<-correlation(data_xinxueguan_QC,data_xinxueguan_sva_QC)
mean(data_correlation[,3])
mean(data_correlation[,4])

## RSDs of QCs
RSD_original<-apply(data_xinxueguan_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
RSD_sva<-apply(data_xinxueguan_sva_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))

RSD_original_0.15<-sum(RSD_original<0.15)/length(RSD_original)
RSD_sva_0.15<-sum(RSD_sva<0.15)/length(RSD_sva)

RSD_original_0.30<-sum(RSD_original<0.30)/length(RSD_original)
RSD_sva_0.30<-sum(RSD_sva<0.30)/length(RSD_sva)



##################### quantile normalization ####################
library(preprocessCore)
data_xinxueguan_quantile<-normalize.quantiles(t(data_xinxueguan_new),copy=F)

# PCA score plot
bitmap(file="PCA_quantile.jpeg",res=800,type="jpeg")
pca(t(data_xinxueguan_quantile),id=group,plot=1)
dev.off()

# Pearson correlation coefficient of QCs
data_xinxueguan_quantile_QC<-t(data_xinxueguan_quantile)[group==3,]

data_correlation_quantile<-correlation(data_xinxueguan_QC,data_xinxueguan_quantile_QC)
mean(data_correlation_quantile[,3])
mean(data_correlation_quantile[,4])

# RSDs of QCs
RSD_quantile<-apply(data_xinxueguan_quantile_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
RSD_quantile_0.15<-sum(RSD_quantile<0.15)/length(RSD_quantile)
RSD_quantile_0.30<-sum(RSD_quantile<0.30)/length(RSD_quantile)


############# Regression Calibration Method  #############
sample_infor<-read.csv("sample.information.csv",header=T)
data_xinxueguan_tran<-t(data_xinxueguan)
data_xinxueguan_tran<-as.data.frame(data_xinxueguan_tran)
data_xinxueguan_tran$sample<-rownames(data_xinxueguan_tran)
data_xinxueguan_merge<-base::merge(data_xinxueguan_tran,sample_infor,by="sample")
rownames(data_xinxueguan_merge)<-data_xinxueguan_merge[,1]
data_xinxueguan_merge<-data_xinxueguan_merge[,-1]
data_xinxueguan_new_R<-data_xinxueguan_merge[,-which(colnames(data_xinxueguan_merge)=="order")]


group_R<-vector(length=162,mode="numeric")
group_R[substr(rownames(data_xinxueguan_new_R),1,2)=="qc"]<-3
group_R[substr(rownames(data_xinxueguan_new_R),1,7)=="disease"]<-2
group_R[substr(rownames(data_xinxueguan_new_R),1,7)=="control"]<-1

## separated into QCs and subject samples
data_xinxueguan_QC_R<-data_xinxueguan_merge[group_R==3,]
data_xinxueguan_sample_R<-data_xinxueguan_merge[group_R!=3,]

group_QC_R<-group_R[group_R==3]
group_sample_R<-group_R[group_R!=3]

### for QC ####
data_xinxueguan_regression_QC<-matrix(NA,40,2316)
for (i in (1:(ncol(data_xinxueguan_QC_R)-1))){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_xinxueguan_QC_R)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"order",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_xinxueguan_QC_R))
  data_xinxueguan_regression_QC[,i]<-data_xinxueguan_QC_R[,i]-predict(model_lm)+mean(data_xinxueguan_QC_R[,i],na.rm=T)
  
}

### For subject samples #####
data_xinxueguan_regression_sample<-matrix(NA,122,2316)
for (i in (1:(ncol(data_xinxueguan_sample_R)-1))){
  cat(paste("#######",i,"###########\n"))
  y<-colnames(data_xinxueguan_sample_R)[i]
  formula<-paste(y,"~",sep="")
  formula1<-paste(formula,"order",sep="")
  model_lm<-lm(formula1,data=as.data.frame(data_xinxueguan_sample_R))
  data_xinxueguan_regression_sample[,i]<-data_xinxueguan_sample_R[,i]-predict(model_lm)+mean(data_xinxueguan_sample_R[,i],na.rm=T)
  
}


### PCA score plot 
data_xinxueguan_regression<-rbind(data_xinxueguan_regression_sample,data_xinxueguan_regression_QC)
label_R<-c(group_sample_R,group_QC_R)
bitmap(file="PCA_regression.jpeg",type="jpeg",res=800)
pca(data_xinxueguan_regression,id=label_R,plot=1)
dev.off()
### Pearson correlation coefficients of QCs
data_regression_correlation<-correlation(data_xinxueguan_QC_R[,-ncol(data_xinxueguan_QC_R)],data_xinxueguan_regression_QC)
mean(data_regression_correlation[,3])
mean(data_regression_correlation[,4])

## RSDs of QCs
RSD_regression<-apply(data_xinxueguan_regression_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))

RSD_regression_0.15<-sum(RSD_regression<0.15)/length(RSD_regression)

RSD_regression_0.30<-sum(RSD_regression<0.30)/length(RSD_regression)



###################### median scaling #############
data_xinxueguan_median<-apply(data_xinxueguan_new,1,function(x) x/median(x,na.rm=T))

# PCA score plot
bitmap(file="PCA_median.jpeg",res=800,type="jpeg")
pca(t(data_xinxueguan_median),id=group,plot=1)
dev.off()

# Pearson correlation coefficients of QCs
data_xinxueguan_median_QC<-t(data_xinxueguan_median)[group==3,]

data_correlation_median<-correlation(data_xinxueguan_QC,data_xinxueguan_median_QC)
mean(data_correlation_median[,3])
mean(data_correlation_median[,4])

# RSDs of QCs
RSD_median<-apply(data_xinxueguan_median_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
RSD_median_0.15<-sum(RSD_median<0.15)/length(RSD_median)
RSD_median_0.30<-sum(RSD_median<0.30)/length(RSD_median)


######### Median Fold Change ############
## for QC
data_xinxueguan_QC<-as.matrix(data_xinxueguan_QC)
median_var_QC<-apply(data_xinxueguan_QC,2,function(x) x/median(x,na.rm=T))
median_sample_QC<-apply(median_var_QC,1,median,na.rm=T)

data_xinxueguan_FC_QC<-matrix(0,40,2316)

for (i in (1:nrow(data_xinxueguan_QC))){
  cat(paste("########",i,"###########\n"))
  data_xinxueguan_FC_QC[i,]<-data_xinxueguan_QC[i,]/median_sample_QC[i]
}

## for sample
data_xinxueguan_sample<-as.matrix(data_xinxueguan_sample)
median_var_sample<-apply(data_xinxueguan_sample,2,function(x) x/median(x,na.rm=T))
median_sample_sample<-apply(median_var_sample,1,median,na.rm=T)

data_xinxueguan_FC_sample<-matrix(0,122,2316)

for (i in (1:nrow(data_xinxueguan_sample))){
  cat(paste("########",i,"###########\n"))
  data_xinxueguan_FC_sample[i,]<-data_xinxueguan_sample[i,]/median_sample_sample[i]
}

### PCA score plot
data_xinxueguan_FC<-rbind(data_xinxueguan_FC_sample,data_xinxueguan_FC_QC)
label<-c(group[group!=3],group[group==3])

bitmap(file="PCA_FC.jpeg",type="jpeg",res=800)
pca(data_xinxueguan_FC,id=label,plot=1)
dev.off()

### Pearson correlation coefficient of QCs
data_correlation_FC<-correlation(data_xinxueguan_QC,data_xinxueguan_FC_QC)
mean(data_correlation_FC[,3])
mean(data_correlation_FC[,4])

## RSDs of QCs
RSD_fc<-apply(data_xinxueguan_FC_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))

RSD_FC_0.15<-sum(RSD_fc<0.15)/length(RSD_fc)
RSD_FC_0.30<-sum(RSD_fc<0.30)/length(RSD_fc)


############## Combat ##################
library(sva)
data_xinxueguan_eb<-ComBat(dat=t(data_xinxueguan_new),batch=batch)

bitmap(file="PCA_Combat.jpeg",type="jpeg",res=800)
pca(t(data_xinxueguan_eb),id=group,plot=1)
dev.off()

data_xinxueguan_eb_QC<-t(data_xinxueguan_eb)[group==3,]

### Pearson correlation coefficients of QCs
data_correlation_eb<-correlation(data_xinxueguan_QC,data_xinxueguan_eb_QC)
mean(data_correlation_eb[,3])
mean(data_correlation_eb[,4])

## RSDs of QCs
RSD_eb<-apply(data_xinxueguan_eb_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))

RSD_eb_0.15<-sum(RSD_eb<0.15)/length(RSD_eb)
RSD_eb_0.30<-sum(RSD_eb<0.30)/length(RSD_eb)




###################### Figures ##########################
############## RSDs ##########
###### RSDs barplot
library(ggplot2)
##original
# group RSDs
# Original data
library(ggplot2)
data<-RSD_original
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")

x1<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,1000),breaks=c(0,200,400,600,800,1000),labels=c(0,200,400,600,800,1000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Raw")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

## MSVA data
data<-RSD_sva
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x2<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,1000),breaks=c(0,200,400,600,800,1000),labels=c(0,200,400,600,800,1000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="MSVA")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

## Combat
data<-RSD_eb
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x3<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,1000),breaks=c(0,200,400,600,800,1000),labels=c(0,200,400,600,800,1000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Combat")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

##quantile
data<-RSD_quantile
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")

x4<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,1000),breaks=c(0,200,400,600,800,1000),labels=c(0,200,400,600,800,1000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Quantile")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

#### Regression
data<-RSD_regression
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x5<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,1000),breaks=c(0,200,400,600,800,1000),labels=c(0,200,400,600,800,1000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Regression")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))


#### Median scaling
data<-RSD_median
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x6<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,1000),breaks=c(0,200,400,600,800,1000),labels=c(0,200,400,600,800,1000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Median scaling")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

### Median fold change

data<-RSD_fc
RSD_data<-RSD_group(data)
RSD_data<-as.data.frame(RSD_data)
RSD_data<-cbind(RSD_data,color_group(data,0.30))
colnames(RSD_data)<-c("x","col")


x7<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
  scale_y_continuous(limits=c(0,1000),breaks=c(0,200,400,600,800,1000),labels=c(0,200,400,600,800,1000))+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title="Median fold change")+
  theme(axis.text.y=element_text(face="bold"),axis.text.x=element_text(angle=45,face="bold"),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))

#### Put them together


bitmap(file="RSD_barplot_zong.jpeg",type="jpeg",res=800,height=7,width=12)
library(grid)
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=4)))
print(x1,vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
print(x2,vp=viewport(layout.pos.row = 1,layout.pos.col = 2))
print(x6,vp=viewport(layout.pos.row = 1,layout.pos.col = 3))

print(x7,vp=viewport(layout.pos.row = 1,layout.pos.col = 4))
print(x4,vp=viewport(layout.pos.row = 2,layout.pos.col = 1))
print(x5,vp=viewport(layout.pos.row = 2,layout.pos.col = 2))
print(x3,vp=viewport(layout.pos.row = 2,layout.pos.col = 3))

dev.off()
############################################
##### RSD percentchart ##########
data<-RSD_original
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Raw",nrow(RSD_data)))
RSD_data1<-RSD_data
rownames(RSD_data1)<-names(RSD_original)

data<-RSD_sva
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("MSVA",nrow(RSD_data)))
RSD_data2<-RSD_data
rownames(RSD_data2)<-names(RSD_sva)



data<-RSD_eb
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Combat",nrow(RSD_data)))
RSD_data3<-RSD_data
rownames(RSD_data3)<-names(RSD_eb)


data<-RSD_quantile
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Quantile",nrow(RSD_data)))
RSD_data4<-RSD_data
rownames(RSD_data4)<-names(RSD_quantile)


data<-RSD_regression
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Regression",nrow(RSD_data)))
RSD_data5<-RSD_data
rownames(RSD_data5)<-names(RSD_regression)

data<-RSD_median
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Median scaling",nrow(RSD_data)))
RSD_data6<-RSD_data
rownames(RSD_data6)<-names(RSD_median)

data<-RSD_fc
RSD_data<-RSD_group3(data)
RSD_data<-as.data.frame(RSD_data)
colnames(RSD_data)<-c("x")
RSD_data$group<-c(rep("Median fold change",nrow(RSD_data)))
RSD_data7<-RSD_data
rownames(RSD_data7)<-names(RSD_fc)

data_integral<-rbind(RSD_data1,RSD_data2,RSD_data3,RSD_data4,RSD_data5,RSD_data6,RSD_data7)

bitmap(file="RSD percent barPlot.jpeg",type="jpeg",res=800,height=8,width=9)
ggplot(data=data_integral,aes(x=group,fill=factor(x)))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#619CFF","#00BA38","#F8766D"),name="RSD",guide=guide_legend(reverse=T),label=c(">30%","15-30%","0-15%"))+
  theme(axis.text=element_text(size=rel(1.2),face="bold"),axis.title=element_text(size=rel(1.5),face="bold"))+
  theme(legend.background   = element_rect(fill="gray90"),legend.title=element_text(size=rel(1.2),face="bold"),legend.text =element_text(size=rel(1.2),face="bold"))+
  labs(x="")+
  scale_x_discrete(limits=c("Raw","MSVA","Median scaling","Median fold change","Quantile","Regression","Combat"),labels=c("Raw","MSVA","Median\n scaling","Median\n fold change","Quantile","Regression","Combat"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,25,50,75,100),name="% of Peaks")+
  guides(fill=F)


dev.off()


#### PCA score plot

### Original data
library(ggfortify)
data_auto_PCA_original<-cbind(group=group,data_xinxueguan_new)
data_auto_PCA_original$group<-as.factor(data_auto_PCA_original$group)
PC_original<-prcomp(scale(data_xinxueguan_new),retx=T)


s1<-autoplot(PC_original,data=data_auto_PCA_original,colour="group")+scale_colour_discrete(name="",label=c("AS","UA","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Raw")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )

### MSVA
library(ggfortify)
data_auto_PCA_sva<-cbind(group=label,data_xinxueguan_sva)
data_auto_PCA_sva<-as.data.frame(data_auto_PCA_sva)
data_auto_PCA_sva$group<-as.factor(data_auto_PCA_sva$group)
PC_sva<-prcomp(scale(data_xinxueguan_sva),retx=T)


s2<-autoplot(PC_sva,data=data_auto_PCA_sva,colour="group")+scale_colour_discrete(name="",label=c("AS","UA","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="MSVA")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )

###Quantile 
library(ggfortify)
data_auto_PCA_quantile<-cbind(group=group,t(data_xinxueguan_quantile))
data_auto_PCA_quantile<-as.data.frame(data_auto_PCA_quantile)
data_auto_PCA_quantile$group<-as.factor(data_auto_PCA_quantile$group)
PC_quantile<-prcomp(scale(t(data_xinxueguan_quantile)),retx=T)


s3<-autoplot(PC_quantile,data=data_auto_PCA_quantile,colour="group")+scale_colour_discrete(name="",label=c("AS","UA","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Quantile")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )


###Regression

library(ggfortify)
data_auto_PCA_regression<-cbind(group=label_R,data_xinxueguan_regression)
data_auto_PCA_regression<-as.data.frame(data_auto_PCA_regression)
data_auto_PCA_regression$group<-as.factor(data_auto_PCA_regression$group)
PC_regression<-prcomp(scale(data_xinxueguan_regression),retx=T)

s4<-autoplot(PC_regression,data=data_auto_PCA_regression,colour="group")+scale_colour_discrete(name="",label=c("AS","UA","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Regression")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )


###Median scaling

library(ggfortify)
data_auto_PCA_median<-cbind(group=group,t(data_xinxueguan_median))
data_auto_PCA_median<-as.data.frame(data_auto_PCA_median)
data_auto_PCA_median$group<-as.factor(data_auto_PCA_median$group)
PC_median<-prcomp(scale(t(data_xinxueguan_median)),retx=T)

s5<-autoplot(PC_median,data=data_auto_PCA_median,colour="group")+scale_colour_discrete(name="",label=c("AS","UA","QC"))+
  theme(legend.position=c(1,1.1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Median scaling")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.6) )

### Median fold change

library(ggfortify)
data_auto_PCA_FC<-cbind(group=label,data_xinxueguan_FC)
data_auto_PCA_FC<-as.data.frame(data_auto_PCA_FC)
data_auto_PCA_FC$group<-as.factor(data_auto_PCA_FC$group)
PC_FC<-prcomp(scale(data_xinxueguan_FC),retx=T)

s6<-autoplot(PC_FC,data=data_auto_PCA_FC,colour="group")+scale_colour_discrete(name="",label=c("AS","UA","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Median fold change")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )

###Combat
library(ggfortify)
data_auto_PCA_combat<-cbind(group=group,t(data_xinxueguan_eb))
data_auto_PCA_combat<-as.data.frame(data_auto_PCA_combat)
data_auto_PCA_combat$group<-as.factor(data_auto_PCA_combat$group)
PC_combat<-prcomp(scale(t(data_xinxueguan_eb)),retx=T)


s7<-autoplot(PC_combat,data=data_auto_PCA_combat,colour="group")+scale_colour_discrete(name="",label=c("AS","UA","QC"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(0.3,"cm"),legend.key.height = unit(0.3,"cm"))+
  labs(title="Combat")+
  theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )

### Put them together
bitmap(file="PCA_zong.jpeg",type="jpeg",res=800,height=7,width=12)
library(grid)
pushViewport(viewport(layout=grid.layout(nrow=2,ncol=4)))
print(s1,vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
print(s2,vp=viewport(layout.pos.row = 1,layout.pos.col = 2))
print(s5,vp=viewport(layout.pos.row = 1,layout.pos.col = 3))

print(s6,vp=viewport(layout.pos.row = 1,layout.pos.col = 4))
print(s3,vp=viewport(layout.pos.row = 2,layout.pos.col = 1))
print(s4,vp=viewport(layout.pos.row = 2,layout.pos.col = 2))
print(s7,vp=viewport(layout.pos.row = 2,layout.pos.col = 3))


dev.off()

####################### Pearson correlation coefficients
#################  Heatmaps
### Original
rownames(data_xinxueguan_QC)<-paste("QC",1:40,sep="")
cor_original<-cor(t(data_xinxueguan_QC))
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
bitmap(file="heatmap_original.jpeg",type="jpeg",res=1000)
autoplot(cor_original,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:40,sep=""))+
  scale_y_discrete(limits=paste("QC",1:40,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)
dev.off()



### SVA 
rownames(data_xinxueguan_sva_QC)<-paste("QC",1:40,sep="")
cor_sva<-cor(t(data_xinxueguan_sva_QC))
bitmap(file="heatmap_msva.jpeg",type="jpeg",res=1000)
autoplot(cor_sva,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:40,sep=""))+
  scale_y_discrete(limits=paste("QC",1:40,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)
dev.off()


###### Quantile
rownames(data_xinxueguan_quantile_QC)<-paste("QC",1:40,sep="")
cor_quantile<-cor(t(data_xinxueguan_quantile_QC))
bitmap(file="heatmap_quantile.jpeg",type="jpeg",res=1000)
autoplot(cor_quantile,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:40,sep=""))+
  scale_y_discrete(limits=paste("QC",1:40,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)
dev.off()

#### Regression
rownames(data_xinxueguan_regression_QC)<-paste("QC",1:40,sep="")
cor_regression<-cor(t(data_xinxueguan_regression_QC))
bitmap(file="heatmap_regression.jpeg",type="jpeg",res=1000)
autoplot(cor_regression,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:40,sep=""))+
  scale_y_discrete(limits=paste("QC",1:40,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)
dev.off()


#### Median scaling
rownames(data_xinxueguan_median_QC)<-paste("QC",1:40,sep="")
cor_median<-cor(t(data_xinxueguan_median_QC))
bitmap(file="heatmap_median.jpeg",type="jpeg",res=1000)
autoplot(cor_median,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:40,sep=""))+
  scale_y_discrete(limits=paste("QC",1:40,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)
dev.off()


#### Median fold change
rownames(data_xinxueguan_FC_QC)<-paste("QC",1:40,sep="")
cor_FC<-cor(t(data_xinxueguan_FC_QC))
bitmap(file="heatmap_FC.jpeg",type="jpeg",res=1000)
autoplot(cor_FC,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:40,sep=""))+
  scale_y_discrete(limits=paste("QC",1:40,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)
dev.off()



#### Combat
rownames(data_xinxueguan_eb_QC)<-paste("QC",1:40,sep="")
cor_combat<-cor(t(data_xinxueguan_eb_QC))
bitmap(file="heatmap_combat1.jpeg",type="jpeg",res=1000)
autoplot(cor_combat,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
  theme(axis.text.x=element_text(angle=90))+
  scale_x_discrete(limits=paste("QC",1:40,sep=""))+
  scale_y_discrete(limits=paste("QC",1:40,sep=""))+
  theme(axis.text=element_text(face="bold"))+
  guides(fill=F)

dev.off()




######## Pearson correlation coefficient of selected QCs
##### Original
data_xinxueguan_QC_tran<-t(data_xinxueguan_QC[c(6,13,25,30,33,36),])
original_cor<-cor(data_xinxueguan_QC_tran)

bitmap(file="scatterPlot_Original.jpeg",type="jpeg",res=1000)
pairs(~QC6+QC13+QC25+QC30+QC33+QC36,data=log(data_xinxueguan_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Original.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))

dev.off()

#### MSVA
data_xinxueguan_QC_tran<-t(data_xinxueguan_sva_QC[c(6,13,25,30,33,36),])
original_cor<-cor(data_xinxueguan_QC_tran)

bitmap(file="scatterPlot_MSVA.jpeg",type="jpeg",res=1000)
pairs(~QC6+QC13+QC25+QC30+QC33+QC36,data=log(data_xinxueguan_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_MSVA.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

#### Quantile
data_xinxueguan_QC_tran<-t(data_xinxueguan_quantile_QC[c(6,13,25,30,33,36),])
original_cor<-cor(data_xinxueguan_QC_tran)

bitmap(file="scatterPlot_Quantile.jpeg",type="jpeg",res=1000)
pairs(~QC6+QC13+QC25+QC30+QC33+QC36,data=log(data_xinxueguan_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Quantile.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()


### Regression
data_xinxueguan_QC_tran<-t(data_xinxueguan_regression_QC[c(6,13,25,30,33,36),])
original_cor<-cor(data_xinxueguan_QC_tran)

bitmap(file="scatterPlot_Regression.jpeg",type="jpeg",res=1000)
pairs(~QC6+QC13+QC25+QC30+QC33+QC36,data=log(data_xinxueguan_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Regression.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

### Median scaling
data_xinxueguan_QC_tran<-t(data_xinxueguan_median_QC[c(6,13,25,30,33,36),])
original_cor<-cor(data_xinxueguan_QC_tran)

bitmap(file="scatterPlot_Median.jpeg",type="jpeg",res=1000)
pairs(~QC6+QC13+QC25+QC30+QC33+QC36,data=log(data_xinxueguan_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Median.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

### Median fold change
data_xinxueguan_QC_tran<-t(data_xinxueguan_FC_QC[c(6,13,25,30,33,36),])
original_cor<-cor(data_xinxueguan_QC_tran)

bitmap(file="scatterPlot_FC.jpeg",type="jpeg",res=1000)
pairs(~QC6+QC13+QC25+QC30+QC33+QC36,data=log(data_xinxueguan_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_FC.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()

###### Combat
data_xinxueguan_QC_tran<-t(data_xinxueguan_eb_QC[c(6,13,25,30,33,36),])
original_cor<-cor(data_xinxueguan_QC_tran)

bitmap(file="scatterPlot_Combat.jpeg",type="jpeg",res=1000)
pairs(~QC6+QC13+QC25+QC30+QC33+QC36,data=log(data_xinxueguan_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
dev.off()

bitmap(file="COR_Combat.jpeg",type="jpeg",res=1000)
library(corrgram)
corrgram(original_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2,col.regions=colorRampPalette(c("grey0", "grey9", "grey25", "black")))
dev.off()



