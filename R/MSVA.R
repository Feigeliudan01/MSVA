#' @title MSVA
#'
#' @description A Novel Approach to Removing Unwanted Variation in Metabolomics Data
#'             Based on Surrogate Variable Analysis
#'
#' @param data The rows of the data represent the samples and the columns of the data
#' represent the variables. The variables should contain "group", "batch" and "injection.order".
#' And other variables are metabolites. The values of variable "group" should be "0", "1" and "QC" with
#' "0" represents control groups and "1" represents disease groups.
#'
#'@param seed
#'@param plot Whether the plots are generated or not.
#'@param QC.picked Choose some QC samples to plot scatterplotMatrix and correlation diagrams so that the performances of the
#' methods can be visualized.
#' @return label_msva data_calibration
#'
#'
#' @examples
#'
#' @export MSVA

########################### MSVA R function ##################
###### MSVA
MSVA<-function(data,seed=1234,plot=T,QC.picked=c(1,2,3,4)){
  library(sva)
  library(grid)
  library(graphics)
  library(corrgram)
  set.seed(seed)
  if (sum(!is.na(match(colnames(data),c("group","batch","injection.order"))))!=3){
    stop("group,batch and injection.order should be one of the colnames of the data")
  }
  if (!is.character(data$group)){
    stop("The variable of group should be character")
  }
  if (sum(!is.na(match(levels(factor(data$group)),c("0","1","QC"))))!=3){
    stop("The levels of groups should be 0,1 and QC")
  }
  data$group[data$group=="QC"]<-2
  group<-data$group
  batch<-data$batch
  injection_order<-data$injection.order
  data_msva<-data[,-which(colnames(data) %in% c("group","batch","injection.order"))]

  if (length(levels(factor(batch)))>1) {
    mod = model.matrix(~as.factor(group)+as.factor(batch),data=data_msva)
    mod0 = model.matrix(~as.factor(batch),data=data_msva)
    svobj<-sva(t(data_msva),mod=mod,mod0=mod0)
    data_model<-cbind(data_msva,svobj$sv,label=as.matrix(group),batch=as.matrix(batch))
    colnames(data_model)[-c(1:ncol(data_msva))]<-c(paste("PC",1:svobj$n.sv,sep=""),"label","batch")
    #### For QCs
    data_model_QC<-data_model[group==2,]
    data_model_sample<-data_model[group!=2,]
    ###MSVA of QC
    data_msva_QC<-matrix(NA,nrow(data_model_QC),ncol(data_msva))
    for (i in (1:ncol(data_msva))){
      cat(paste("#######","QC calibration",i,"###########\n"))
      y<-colnames(data_model_QC)[i]
      formula<-paste(y,"~",sep="")
      formula1<-paste(formula,paste("PC",1:svobj$n.sv,sep=""),"+batch",sep="")
      model_lm<-lm(formula1,data=as.data.frame(data_model_QC))
      data_msva_QC[,i]<-data_model_QC[,i]-predict(model_lm)+mean(data_model_QC[,i],na.rm=T)
    }
    ###MSVA of sample
    data_msva_sample<-matrix(NA,nrow(data_model_sample),ncol(data_msva))
    for (i in (1:ncol(data_msva))){
      cat(paste("#######","sample calibration",i,"###########\n"))
      y<-colnames(data_model_sample)[i]
      formula<-paste(y,"~",sep="")
      formula1<-paste(formula,paste("PC",1:svobj$n.sv,sep=""),"+batch",sep="")
      model_lm<-lm(formula1,data=as.data.frame(data_model_sample))
      data_msva_sample[,i]<-data_model_sample[,i]-predict(model_lm,newdata=data_model_sample[,c(i,(ncol(data_msva)+1):ncol(data_model_sample))])+mean(data_model_sample[,i],na.rm=T)
    }
  }
  if (length(levels(factor(batch)))==1) {
    mod = model.matrix(~as.factor(group),data=data_msva)
    mod0 = model.matrix(~1,data=data_msva)
    svobj<-sva(t(data_msva),mod=mod,mod0=mod0)
    data_model<-cbind(data_msva,svobj$sv,label=as.matrix(group),batch=as.matrix(batch))
    colnames(data_model)[-c(1:ncol(data_msva))]<-c(paste("PC",1:svobj$n.sv,sep=""),"label","batch")
    #### For QCs
    data_model_QC<-data_model[group==2,]
    data_model_sample<-data_model[group!=2,]
    ###MSVA of QC
    data_msva_QC<-matrix(NA,nrow(data_model_QC),ncol(data_msva))
    for (i in (1:ncol(data_msva))){
      cat(paste("#######","QC",i,"###########\n"))
      y<-colnames(data_model_QC)[i]
      formula<-paste(y,"~",sep="")
      formula1<-paste(formula,paste("PC",1:svobj$n.sv,sep=""),sep="")
      model_lm<-lm(formula1,data=as.data.frame(data_model_QC))
      data_msva_QC[,i]<-data_model_QC[,i]-predict(model_lm)+mean(data_model_QC[,i],na.rm=T)
    }
    ###MSVA of sample
    data_msva_sample<-matrix(NA,nrow(data_model_sample),ncol(data_msva))
    for (i in (1:ncol(data_msva))){
      cat(paste("#######","sample",i,"###########\n"))
      y<-colnames(data_model_sample)[i]
      formula<-paste(y,"~",sep="")
      formula1<-paste(formula,paste("PC",1:svobj$n.sv,sep=""),sep="")
      model_lm<-lm(formula1,data=as.data.frame(data_model_sample))
      data_msva_sample[,i]<-data_model_sample[,i]-predict(model_lm,newdata=data_model_sample[,c(i,(ncol(data_msva)+1):ncol(data_model_sample))])+mean(data_model_sample[,i],na.rm=T)
    }
  }
  injection_order_msva<-c(injection_order[group!="2"],injection_order[group=="2"])
  group_msva<-c(group[group!="2"],group[group=="2"])
  data_calibration_msva<-rbind(data_msva_sample,data_msva_QC)

  if (plot){
    cat(paste("### Ploting PCA score plots #####\n"))
    pca_MSVA<-PCA_ScorePlot(data=data_calibration_msva,id=group_msva,main="MSVA")
    pca_original<-PCA_ScorePlot(data=data_msva,id=group,main="RAW")

    dir.create("plot")
    bitmap(file="plot/PCA_scoreplot_Raw.jpeg",res=1000,type="jpeg")
    print(pca_original)
    dev.off()

    bitmap(file="plot/PCA_scoreplot_msva.jpeg",res=1000,type="jpeg")
    print(pca_MSVA)
    dev.off()

    cat(paste("### Ploting heatmaps of Pearson correlatio coefficients #####\n"))
    ##raw data
    data_raw_QC<-data_msva[group=="2",]
    injection_order_QC<-injection_order[group=="2"]

    data_raw_QC_order<-data_raw_QC[order(injection_order_QC),]

    rownames(data_raw_QC_order)<-paste("QC",1:nrow(data_raw_QC_order),sep="")

    bitmap(file="plot/RAW_Heatmap.jpeg",res=1000,type="jpeg")
    print(heatmap_QC(data_raw_QC_order))
    dev.off()

    ##msva data
    data_msva_QC<-data_calibration_msva[group_msva==2,]
    injection_order_qc_msva<-injection_order_msva[group_msva==2]
    data_msva_QC_order<-data_msva_QC[order(injection_order_qc_msva),]
    rownames(data_msva_QC_order)<-paste("QC",1:nrow(data_msva_QC_order),sep="")

    bitmap(file="plot/msva_Heatmap.jpeg",res=1000,type="jpeg")
    print(heatmap_QC(data_msva_QC_order))
    dev.off()

    cat(paste("### Ploting batplot of RSDS of QCs #####\n"))
    RAW_barplot<-barplot_RSD(data_QC=data_raw_QC_order,main="RAW")
    MSVA_barplot<-barplot_RSD(data_QC=data_msva_QC_order,main="MSVA")

    bitmap(file="plot/batplot_of_RSD.jpeg",res=1000,type="jpeg")
    pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1)))
    print(RAW_barplot,vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
    print(MSVA_barplot,vp=viewport(layout.pos.row = 2,layout.pos.col = 1))
    dev.off()

    cat(paste("## Ploting ScatterplotMatrix and correlation diagram for picked QC samples ##\n"))

    ####Raw data
    data_raw_QC_tran<-t(data_raw_QC_order[QC.picked,])
    raw_cor<-cor(data_raw_QC_tran)

    formula<-paste("QC",QC.picked,sep="",collapse = "+")
    formula1<-paste("~",formula,sep="")
    bitmap(file="plot/scatterPlot_raw.jpeg",type="jpeg",res=1000)
    pairs(x=log(data_raw_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
    dev.off()

    bitmap(file="plot/COR_raw.jpeg",type="jpeg",res=1000)
    corrgram(raw_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2)
    dev.off()

    ### MSVA
    data_msva_QC_tran<-t(data_msva_QC_order[QC.picked,])
    msva_cor<-cor(data_msva_QC_tran)

    bitmap(file="plot/scatterPlot_msva.jpeg",type="jpeg",res=1000)
    pairs(x=log(data_msva_QC_tran),pch=16,gap=0.5,font.labels = 2,cex=0.8)
    dev.off()

    bitmap(file="plot/COR_msva.jpeg",type="jpeg",res=1000)
    corrgram(msva_cor,type="cor",upper.panel =panel.pie,lower.panel = panel.cor,font.labels = 2)
    dev.off()
  }
  group_msva[group_msva=="2"]<-"QC"
  return(list(label_msva=group_msva,data_calibration=data_calibration_msva))
}

#### PCA score plot
PCA_ScorePlot<-function(data,id,label=c("control","disease","QC"),main=""){
  library(ggfortify)
  data_PCA<-cbind(group=id,data)
  data_PCA<-as.data.frame(data_PCA,stringsAsFactors = F)
  data_PCA$group<-as.factor(data_PCA$group)
  PC<-prcomp(scale(data),retx=T)
  X<-autoplot(PC,data=data_PCA,colour="group")+scale_colour_discrete(name="",label=label)+
    theme(legend.position=c(1,1),legend.justification = c(1,1),legend.key.width=unit(2.5,"cm"),legend.key.height = unit(1,"cm"))+
    labs(title=main)+
    theme(plot.title =element_text(face="bold",size=rel(1.5),hjust=0.5) )
  return(X)
}

### heatmaps of Pearson correlation coefficients
heatmap_QC<-function(data_QC){
  if (is.null(rownames(data_QC))){
    stop("data should have rownames")
  }
  cor_qc<-cor(t(data_QC))
  library(RColorBrewer)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  X<-autoplot(cor_qc,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.85,1),name="r value")+labs(x="",y="")+
    theme(axis.text.x=element_text(angle=90))+
    scale_x_discrete(limits=paste("QC",1:nrow(data_QC),sep=""))+
    scale_y_discrete(limits=paste("QC",1:nrow(data_QC),sep=""))+
    theme(axis.text=element_text(face="bold"))+
    guides(fill=F)
  return(X)
}

### Barplot of RSDs of QC samples
barplot_RSD<-function(data_QC,main=""){
  RSD<-apply(data_QC,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
  data<-RSD
  RSD_data<-RSD_group(data)
  RSD_data<-as.data.frame(RSD_data)
  RSD_data<-cbind(RSD_data,color_group(data,0.30))
  colnames(RSD_data)<-c("x","col")
  X<-ggplot(data=RSD_data,aes(x=x,fill=col))+geom_bar()+scale_fill_manual(values=c("grey","red"))+guides(fill=F)+
    scale_x_discrete(limits=c(1,2,3,4,5,6,7),labels=c("0-15%","15-30%","30-45%","45-60%","60-75%","75-90%",">90%"))+labs(x="%RSD",y="Peak Number",title=main)+
    theme(axis.text.y=element_text(face="bold",size=rel(1.2)),axis.text.x=element_text(angle=45,face="bold",size=rel(1.2)),axis.title=element_text(size=rel(1.2),face="bold"),plot.title=element_text(size=rel(1.5),face="bold",hjust=0.5))
  return(X)
}


###RSD分组分7组
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

##RSD 分三组
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

###颜色分组
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
