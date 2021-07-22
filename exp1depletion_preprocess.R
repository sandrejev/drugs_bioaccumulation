library(plyr)
library(readr)
library(ggplot2)
library(xlsx)
library(parallel)
source("functions.R")


get.dataclean <- function(data.combined){
  data.clean <- data.combined[,c("Plate","GroupName","drug","DrugRatio","Status","Well","Batch","Replicate","IS","Species.x","Plate_order","growth","maxOD")] #"max.blank","max.time","mu.blank","mu.finish","lag","auc"
  data.clean <- data.clean[!data.clean$Plate=="BuniformisA",]
  data.clean$DrugRatio[which(data.clean$DrugRatio>9)] <- Inf #######Attention!!! limit introduced!!!
  
  ####implicitely kicking out everything which does grow in control if not using all.x=T
  data.clean <- plyr::ddply(data.clean, .(GroupName,Plate), mutate, CtrlMean= rep(mean(DrugRatio[which(Status=="ctrl")],na.rm=T),length(DrugRatio)))
  
  data.clean$NormedCtrl <- data.clean$DrugRatio/data.clean$CtrlMean
  
  batchmean.ctrls <- plyr::ddply(data.clean,c("Batch","GroupName"), function(z) mean(z$DrugRatio[which(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0 & z$growth=="no growth")],na.rm=T))
  colnames(batchmean.ctrls)[3] <- "BatchMean_ctrls"
  data.clean <- merge(data.clean,batchmean.ctrls, all.x=T)
  data.clean$DiffToBatchMean <- (data.clean$DrugRatio-data.clean$BatchMean_ctrls)/data.clean$BatchMean_ctrls
  
  allmean.ctrls <- plyr::ddply(data.clean[which(data.clean$Status=="ctrl" & is.finite(data.clean$DrugRatio) & data.clean$DrugRatio>0 & data.clean$growth=="no growth" & !data.clean$Batch==3),],c("GroupName"), function(z) mean(z$DrugRatio,na.rm=T))
  colnames(allmean.ctrls)[2] <- "AllMean_ctrls"
  data.clean <- merge(data.clean,allmean.ctrls, all.x=T)
  data.clean$DiffToAllMean <- (data.clean$DrugRatio-data.clean$AllMean_ctrls)/data.clean$AllMean_ctrls
  
  batchmedian.ctrls <- plyr::ddply(data.clean,c("Batch","GroupName"), function(z) median(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0 & z$growth=="no growth")],na.rm=T))
  colnames(batchmedian.ctrls)[3] <- "BatchMedian_ctrls"
  data.clean <- merge(data.clean,batchmedian.ctrls, all.x=T)
  data.clean$DiffToBatchMedian <- (data.clean$DrugRatio-data.clean$BatchMedian_ctrls)/data.clean$BatchMedian_ctrls
  
  sd.ctrls <- plyr::ddply(data.clean, c("Batch","GroupName"), function(z) sd(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0 & z$growth=="no growth")],na.rm=T))
  colnames(sd.ctrls)[3] <- "SD_ctrls"
  data.clean <- merge(data.clean,sd.ctrls, all.x=T)
  
  n.ctrls <- plyr::ddply(data.clean, c("Batch","GroupName"), function(z) length(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0 & z$growth=="no growth")]))
  colnames(n.ctrls)[3] <- "N_ctrls"
  data.clean <- merge(data.clean,n.ctrls, all.x=T)
  
  sdall.ctrls <- plyr::ddply(data.clean, c("GroupName"), function(z) sd(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0 & z$growth=="no growth" & !z$Batch==3)],na.rm=T))
  colnames(sdall.ctrls)[2] <- "SDall_ctrls"
  data.clean <- merge(data.clean,sdall.ctrls, all.x=T)
  
  nall.ctrls <- plyr::ddply(data.clean[(data.clean$Status=="ctrl" & is.finite(data.clean$DrugRatio) & data.clean$DrugRatio>0 & data.clean$growth=="no growth" & !data.clean$Batch==3),], c("GroupName"), function(z) length(z$DrugRatio))
  colnames(nall.ctrls)[2] <- "Nall_ctrls"
  data.clean <- merge(data.clean,nall.ctrls, all.x=T)
  
  data.clean$Norm.Max <- data.clean$DiffToBatchMean/(data.clean$maxOD)
  data.clean$Norm.Max.Median <- data.clean$DiffToBatchMedian/(data.clean$maxOD)
  #data.clean$Norm.Time <- data.clean$DiffToBatchMean/data.clean$max.time
  
  ###kich out things
  data.clean <- data.clean[which(is.finite(data.clean$DrugRatio)),]
  
  data.ctrls <- data.clean[which(data.clean$Status=="ctrl" & data.clean$DrugRatio>0 & is.finite(data.clean$DrugRatio) & data.clean$growth=="no growth"),]
  
  batch.means <- plyr::ddply(data.ctrls, .(Batch,GroupName), summarize, BatchMean= mean(DrugRatio,na.rm=T))
  data.clean <- merge(data.clean,batch.means, all.x=T)
  data.clean$DiffToBatchMean <- (data.clean$DrugRatio-data.clean$BatchMean)/data.clean$BatchMean  
  
  sd.diff.ctrls <- plyr::ddply(data.ctrls, c("Batch","GroupName"), summarize, SD.diff_ctrls=sd(DiffToBatchMean,na.rm=T))
  data.clean <- merge(data.clean,sd.diff.ctrls, all.x=T)
  n.diff.ctrls <- plyr::ddply(data.ctrls[which(!is.na(data.ctrls$DiffToBatchMean)),], c("Batch","GroupName"), summarize, N.diff_ctrls=length(DiffToBatchMean))
  data.clean <- merge(data.clean,n.diff.ctrls, all.x=T)
  
  sd.ctrls <- plyr::ddply(data.ctrls, c("Batch","GroupName"), summarize, SD.ctrls=sd(DrugRatio,na.rm=T))
  data.clean <- merge(data.clean,sd.ctrls, all.x=T)
  n.ctrls <- plyr::ddply(data.ctrls[which(!is.na(data.ctrls$DrugRatio)),], c("Batch","GroupName"), summarize, N.ctrls=length(DrugRatio))
  data.clean <- merge(data.clean,n.ctrls, all.x=T)
  
  mean.ctrls <- plyr::ddply(data.ctrls,.(Plate,GroupName), summarize, Mean.Ctrl=mean(DrugRatio,na.rm=T))
  data.clean <- (merge(data.clean,mean.ctrls, all.x=T))
  
  
  #####for calculating diftoown, replace DrugRatio in ctrls with growth by mean batch value
  data.dirty <- data.clean
  data.dirty$DrugRatio[which(data.dirty$Status=="ctrl" & data.dirty$growth=="growth")] <- data.dirty$BatchMean[which(data.dirty$Status=="ctrl" & data.dirty$growth=="growth")]
  
  difftoown <- plyr::ddply(data.dirty, .(Plate,GroupName), summarize, 
                           DiffToOwn= (mean(DrugRatio[which(Status=="sample")],na.rm=T)-mean(DrugRatio[which(Status=="ctrl")],na.rm=T))/mean(DrugRatio[which(Status=="ctrl")],na.rm=T),
                           DirtyCtrlMean = mean(DrugRatio[which(Status=="ctrl")],na.rm=T))
  # DiffToOwnRatio= DrugRatio/mean(DrugRatio[which(Status=="sample")]))
  data.clean <- merge(data.clean,difftoown,all.x=T)
  
  data.clean$dummy <- rep(1,dim(data.clean)[1])
  data.clean <- plyr::ddply(data.clean, .(Plate,GroupName), transform, dummy = ifelse(!(sum(DrugRatio[which(Status=="ctrl")])>0),0,dummy) )
  
  data.clean$DiffCtrlSample <- (data.clean$DrugRatio-data.clean$DirtyCtrlMean)/data.clean$DirtyCtrlMean
  
  return(list(data.clean,data.ctrls))
}

get.datalong <- function(data.norm, code){
  data.area = data.norm %>%
    dplyr::mutate(Status=dplyr::case_when(
      grepl("ctrl", SampleName) ~ "ctrl",
      grepl("GMM", SampleName) ~ "GMM",
      T ~ "sample")) %>%
    dplyr::select(SampleName, GroupName, Extraction, Status, WellReplicate, dplyr::matches("Area"))
  colnames(data.area) <- sapply(strsplit(colnames(data.area),".",fixed=T), function(z) as.character(z[[1]]))
  data.area = data.area %>%
    reshape2::melt(id=c("SampleName", "Extraction", "GroupName","Status", "WellReplicate"), variable.name="Plate", value.name="DrugRatio") %>%
    dplyr::mutate(SampleName.long=paste0("Depletion-", Plate, "-", Extraction, "-", GroupName, "-", Status, "-", WellReplicate))
  
  # No code for: Buniformis65A1, Buniformis66A1, LlactisC1, RtorquesC1
  data.long = data.area %>%
    dplyr::inner_join(code, by="Plate") %>%
    dplyr::group_by(Batch, GroupName) %>%
    dplyr::mutate(BatchMean_ctrls=mean(DrugRatio[(Status=="ctrl" & is.finite(DrugRatio) & DrugRatio>0)], na.rm=T)) %>%
    dplyr::mutate(DiffToBatchMean=(DrugRatio-BatchMean_ctrls)/BatchMean_ctrls) %>%
    dplyr::mutate(BatchMedian_ctrls=median(DrugRatio[(Status=="ctrl" & is.finite(DrugRatio) & DrugRatio>0)],na.rm=T)) %>%
    dplyr::mutate(DiffToBatchMedian=(DrugRatio-BatchMedian_ctrls)/BatchMedian_ctrls) %>%
    dplyr::group_by(Batch, GroupName) %>%
    dplyr::mutate(Mean_ctrls=mean(DrugRatio[(Status=="ctrl" & is.finite(DrugRatio) & DrugRatio>0)],na.rm=T)) %>%
    dplyr::mutate(DiffToMean=(DrugRatio-Mean_ctrls)/Mean_ctrls) %>%
    dplyr::mutate(Median_ctrls=median(DrugRatio[(Status=="ctrl" & is.finite(DrugRatio) & DrugRatio>0)],na.rm=T)) %>%
    dplyr::mutate(DiffToMedian=(DrugRatio-Median_ctrls)/Median_ctrls) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(as.character(Plate),GroupName,Status) %>%
    dplyr::group_by(Batch, GroupName, Plate, Status) %>%
    dplyr::mutate(Well=1:length(Status)) %>%
    dplyr::group_by(Batch, GroupName) %>%
    dplyr::mutate(SD_ctrls=sd(DrugRatio[(Status=="ctrl" & is.finite(DrugRatio) & DrugRatio>0)],na.rm=T)) %>%
    dplyr::mutate(N_ctrls=length(DrugRatio[(Status=="ctrl" & is.finite(DrugRatio) & DrugRatio>0)])) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(as.character(Plate), GroupName, Status, Well)
  
  
  
  return(list(data.long,data.area))
  
  
  
  
  # data.area = data.norm[,c(1,2,grep("Area",colnames(data.norm)))]
  # colnames(data.area) <- sapply(strsplit(colnames(data.area),".",fixed=T), function(z) as.character(z[[1]]))
  # colnames(data.area)[1:2] <- c("SampleName","GroupName")
  # status <- rep(c("sample","sample","sample","ctrl","ctrl","GMM"),19)
  # data.area$Status <- status
  
  
  
  # data.area <- data.area[,-1]
  # data.area <- data.area[,c(1,dim(data.area)[2],2:(dim(data.area)[2]-1))]
  # data.area <- reshape2::melt(data.area,id=c("GroupName","Status"),variable.name="Plate",value.name="DrugRatio")
  
  # data.long <- merge(data.area,code,by="Plate")
  batchmean.ctrls <- plyr::ddply(data.long,c("Batch","GroupName"), function(z) mean(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0)],na.rm=T))
  colnames(batchmean.ctrls)[3] <- "BatchMean_ctrls"
  data.long <- merge(data.long,batchmean.ctrls)
  data.long$DiffToBatchMean <- (data.long$DrugRatio-data.long$BatchMean_ctrls)/data.long$BatchMean_ctrls
  
  batchmedian.ctrls <- plyr::ddply(data.long,c("Batch","GroupName"), function(z) median(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0)],na.rm=T))
  colnames(batchmedian.ctrls)[3] <- "BatchMedian_ctrls"
  data.long <- merge(data.long,batchmedian.ctrls)
  data.long$DiffToBatchMedian <- (data.long$DrugRatio-data.long$BatchMedian_ctrls)/data.long$BatchMedian_ctrls
  
  mean.ctrls <- plyr::ddply(data.long,"GroupName", function(z) mean(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0)],na.rm=T))
  colnames(mean.ctrls)[2] <- "Mean_ctrls"
  data.long <- merge(data.long,mean.ctrls)
  data.long$DiffToMean <- (data.long$DrugRatio-data.long$Mean_ctrls)/data.long$Mean_ctrls
  
  median.ctrls <- plyr::ddply(data.long,"GroupName", function(z) median(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0)],na.rm=T))
  colnames(median.ctrls)[2] <- "Median_ctrls"
  data.long <- merge(data.long,median.ctrls)
  data.long$DiffToMedian <- (data.long$DrugRatio-data.long$Median_ctrls)/data.long$Median_ctrls
  
  data.long <- data.long[order(as.character(data.long$Plate),data.long$GroupName,data.long$Status),]
  data.long <- plyr::ddply(data.long, .(Batch, GroupName,Plate,Status), mutate, Well= 1:length(Status))
  data.long <- data.long[order(as.character(data.long$Plate),data.long$GroupName,data.long$Status,data.long$Well),]
  
  sd.ctrls <- plyr::ddply(data.long, c("Batch","GroupName"), function(z) sd(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0)],na.rm=T))
  colnames(sd.ctrls)[3] <- "SD_ctrls"
  data.long <- merge(data.long,sd.ctrls)
  
  n.ctrls <- plyr::ddply(data.long, c("Batch","GroupName"), function(z) length(z$DrugRatio[(z$Status=="ctrl" & is.finite(z$DrugRatio) & z$DrugRatio>0)]))
  colnames(n.ctrls)[3] <- "N_ctrls"
  data.long <- merge(data.long,n.ctrls)
  
  return(list(data.long,data.area))
}

peak.common.nomalize = function(df) 
{
  df[is.na(df[,c("Area")]),c("Area","Height")] <- 0 # TODO: I think this should be NA
  df[which(df$Area=="NA"),c("Area","Height")] <- NA # TODO: I think this should be NA
  
  df_sub1 = df %>%
    dplyr::filter(!is.na(SampleName)) %>%
    dplyr::mutate(Area=as.numeric(Area), Height=as.numeric(Height)) %>%
    dplyr::mutate(Drug=substr(Name,1,5)) %>%
    dplyr::mutate(GroupName=substr(group,3,6)) %>%
    dplyr::mutate(Plate_Martina_Position=substr(SampleName,1,1)) %>%
    dplyr::group_by(Drug, SampleName, Vial) %>%
    dplyr::mutate(
      Area_IS=ifelse(any(grepl("caffe",Name) & !is.na(Area)), unique(Area[grepl("caffe",Name) & !is.na(Area)])*IS_Factor, NA_real_),
      Height_IS=ifelse(any(grepl("caffe",Name) & !is.na(Height)), unique(Height[grepl("caffe",Name) & !is.na(Height)])*IS_Factor, NA_real_)
    ) %>%
    dplyr::group_by(Drug, Plate_Martina_Position) %>%
    dplyr::mutate(
      Area=ifelse(grepl("GMM", SampleName), Area, Area-Area[grepl("GMM", SampleName) & grepl("peak", Name)]),
      Height=ifelse(grepl("GMM", SampleName), Height, Height-Height[grepl("GMM", SampleName) & grepl("peak", Name)])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Area_norm=Area/Area_IS, Height_norm=Height/Height_IS) %>%
    dplyr::filter(grepl("peak", Name)) %>%
    dplyr::select(Plate_Martina, SampleName, Area_norm, Height_norm, GroupName, Area, Height, Area_IS, Height_IS)
  
  # Buniformis65A1, Buniformis65B
  # data.long %>% dplyr::filter(grepl("metf", GroupName) & grepl("ctrl", SampleName) & grepl("Buniformis65", Plate)) %>% View()
  # df %>% dplyr::filter(grepl("metf", group))
  # df_sub1 %>% dplyr::filter(grepl("metf", GroupName))
  # data.clean.old %>%
  #   dplyr::filter(Status=="ctrl" & drug=="metformin" & Plate=="Buniformis65B" & GroupName=="metf")
  # data.clean.new %>%
  #   dplyr::filter(Status=="ctrl" & drug=="metformin" & Plate=="Buniformis65B" & GroupName=="metf")
  
  # table(df_sub1.bck$Area_norm/df_sub1$Area_norm)
  return(df_sub1)
}

peak.IS.normalize = function(df) # a is number of measurements per drug across # of plates e.g. a = 6 measurements*2 plates = 12
{
  dff <<- df
  # df = raw.IS.s[[9]]
  
  
  df <- df[which(!is.na(df$Name)),]
  if(dim(df)[1] > 300){
    a <- 12
  }
  else { a <- 6}
  print(a)
  df[is.na(df[,c("Area")]),c("Area","Height")] <- 0 # TODO: I think this should be NA
  df[which(df$Area=="NA"),c("Area","Height")] <- NA # TODO: I think this should be NA
  
  df_sub1 = df %>%
    dplyr::mutate(Drug=substr(Name,1,5)) %>%
    dplyr::mutate(GroupName=substr(group,3,6)) %>%
    dplyr::mutate(Plate_Martina_Position=substr(SampleName,1,1)) %>%
    dplyr::group_by(Drug, SampleName, Vial) %>%
    dplyr::mutate(
      Area_IS=round(Area[grep("caffe",Name)],0),
      Height_IS=round(Height[grep("caffe",Name)],0)
    ) %>%
    dplyr::group_by(Drug, Plate_Martina_Position) %>%
    dplyr::mutate(
      Area=ifelse(grepl("GMM", SampleName), Area, Area-Area[grepl("GMM", SampleName) & grepl("peak", Name)]),
      Height=ifelse(grepl("GMM", SampleName), Height, Height-Height[grepl("GMM", SampleName) & grepl("peak", Name)])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Area_norm=Area/Area_IS, Height_norm=Height/Height_IS) %>%
    dplyr::filter(grepl("peak", Name)) %>%
    dplyr::select(Plate_Martina, SampleName, Area_norm, Height_norm, GroupName, Area, Height, Area_IS, Height_IS)
  return(df_sub1)
  
  df$Area = as.double(df$Area)
  df$Height = as.double(df$Height)
  df <- df[order(df$Name,df$SampleName,df$Vial),]
  #df_sub <- df[!df$group %in% c("1_tolm","1_rani","1_mont","1_metr","2_tolm","2_rani","2_mont","2_metr"),]
  df_sub <- df
  Plate_Martina <- as.character(df_sub$Plate_Martina[grep("caffe",df_sub$Name)])
  SampleName <- as.character(df_sub$SampleName[grep("caffe",df_sub$Name)])
  GroupName <- as.character(substr(df_sub$group[grep("caffe",df_sub$Name)],3,6))
  IS_area <- round(df_sub$Area[grep("caffe",df_sub$Name)]/2,0)
  #Area_norm <- df_sub$Area[grep("caffe",df_sub$Name)+a]/IS_area
  IS_height <- round(df_sub$Height[grep("caffe",df_sub$Name)]/2,0)
  #Height_norm <- df_sub$Height[grep("caffe",df_sub$Name)+a]/IS_height
  df_sub$Drug <- substr(df_sub$Name,1,5)
  if(a==12){
    Area_corrected <- plyr::ddply(df_sub, .(Drug), function(z) c(z$Area[(grep("peak", z$Name)[1:3])]-z$Area[grep("GMM", z$SampleName)[3]],
                                                                 z$Area[(grep("peak", z$Name)[4:6])],
                                                                 z$Area[(grep("peak", z$Name)[7:9])]-z$Area[grep("GMM", z$SampleName)[4]],
                                                                 z$Area[(grep("peak", z$Name)[10:12])]))
    Height_corrected <- plyr::ddply(df_sub, .(Drug), function(z) c(z$Height[(grep("peak", z$Name)[1:3])]-z$Height[grep("GMM", z$SampleName)[3]],
                                                                   z$Height[(grep("peak", z$Name)[4:6])],
                                                                   z$Height[(grep("peak", z$Name)[7:9])]-z$Height[grep("GMM", z$SampleName)[4]],
                                                                   z$Height[(grep("peak", z$Name)[10:12])]))
  } else{
    Area_corrected <- plyr::ddply(df_sub, .(Drug), function(z) { zz<<-z; c(z$Area[(grep("peak", z$Name)[1:3])]-z$Area[grep("GMM", z$SampleName)[2]],
                                                                           z$Area[(grep("peak", z$Name)[4:6])]) })
    Height_corrected <- plyr::ddply(df_sub, .(Drug), function(z) c(z$Height[(grep("peak", z$Name)[1:3])]-z$Height[grep("GMM", z$SampleName)[2]],
                                                                   z$Height[(grep("peak", z$Name)[4:6])]))
  }
  
  
  Area_corrected = as.vector(t(Area_corrected[,-1])) 
  Height_corrected = as.vector(t(Height_corrected[,-1]))
  
  Area_norm <- Area_corrected/IS_area
  Height_norm <- Area_corrected/IS_height
  
  data <- data.frame(
    Plate_Martina=Plate_Martina, 
    SampleName=SampleName, 
    Area_norm=as.numeric(Area_norm), 
    Height_norm=as.numeric(Height_norm), 
    GroupName=GroupName,
    Area=Area_corrected,
    Height = Height_corrected,
    Area_IS=IS_area,
    Height_IS=IS_height
  )
  data <- data[order(data$SampleName),]
  
  if(sum(!is.na(data$Area_norm[grep("GMM",data$SampleName)])!=0)){
    df.negCtrlpos <- data.frame(SampleName=data$SampleName[grep("GMM",data$SampleName)[!is.na(data$Area_norm[grep("GMM",data$SampleName)])]],Area_norm=data$Area_norm[grep("GMM",data$SampleName)[!is.na(data$Area_norm[grep("GMM",data$SampleName)])]])
    negCtrlpos <- df.negCtrlpos[df.negCtrlpos$Area_norm!=0,]
    print("negative ctrl is positive! check negCtrlpos")  
  }
  return(data)
}

peak.read = function(path, remove.empty=T, verbose=T)
{
  if(verbose) {
    writeLines(paste0("[Reading] '", path, "'"))
  }
  
  df = read.xlsx(path,1,startRow=1,stringsAsFactors=F)[1:456,] %>%
    dplyr::mutate(group = substr(SampleName,1,6)) %>%
    dplyr::arrange(Name,Channel.Description,Sample.Set.Name,SampleName,Vial)
  
  if(verbose) {
    writeLines(paste0("[Done] '", path, "'"))
  }
  
  if (length(unique(df$Sample.Set.Name))>1){
    print("Data processed, but more than 1 sample set detected")
    
  }
  
  return(df)
}

peak.normalize = function(df) 
{
  dff <<- df
  df <- df[which(!is.na(df$Name)),]
  if(dim(df)[1] > 300){
    a <- 12
  }else { a <- 6}  # a is number of measurements per drug across # of plates e.g. a = 6 measurements*2 plates = 12
  print(a)
  print(df$Sample.Set.Name[1])
  df[is.na(df[,c("Area")]),c("Area","Height")] <- 0
  df[which(df$Area=="NA"),c("Area","Height")] <- NA
  df$Area = as.double(df$Area)
  df$Height = as.double(df$Height)
  df <- df[order(df$Name,df$SampleName,df$Vial),]
  #df_sub <- df[!df$group %in% c("1_tolm","1_rani","1_mont","1_metr","2_tolm","2_rani","2_mont","2_metr"),]
  df_sub <- df
  Plate_Martina <- as.character(df_sub$Plate_Martina[grep("caffe",df_sub$Name)])
  SampleName <- as.character(df_sub$SampleName[grep("caffe",df_sub$Name)])
  GroupName <- as.character(substr(df_sub$group[grep("caffe",df_sub$Name)],3,6))
  df_sub$Drug <- substr(df_sub$Name,1,5)
  if(a==12){
    Area_corrected <- plyr::ddply(df_sub, .(Drug), function(z) c(z$Area[(grep("peak", z$Name)[1:3])]-z$Area[grep("GMM", z$SampleName)[3]],
                                                                 z$Area[(grep("peak", z$Name)[4:6])],
                                                                 z$Area[(grep("peak", z$Name)[7:9])]-z$Area[grep("GMM", z$SampleName)[4]],
                                                                 z$Area[(grep("peak", z$Name)[10:12])]))
    Height_corrected <- plyr::ddply(df_sub, .(Drug), function(z) c(z$Height[(grep("peak", z$Name)[1:3])]-z$Height[grep("GMM", z$SampleName)[3]],
                                                                   z$Height[(grep("peak", z$Name)[4:6])],
                                                                   z$Height[(grep("peak", z$Name)[7:9])]-z$Height[grep("GMM", z$SampleName)[4]],
                                                                   z$Height[(grep("peak", z$Name)[10:12])]))
  } else{
    Area_corrected <- plyr::ddply(df_sub, .(Drug), function(z) c(z$Area[(grep("peak", z$Name)[1:3])]-z$Area[grep("GMM", z$SampleName)[2]],
                                                                 z$Area[(grep("peak", z$Name)[4:6])]))
    Height_corrected <- plyr::ddply(df_sub, .(Drug), function(z) c(z$Height[(grep("peak", z$Name)[1:3])]-z$Height[grep("GMM", z$SampleName)[2]],
                                                                   z$Height[(grep("peak", z$Name)[4:6])]))
  }
  
  Area_corrected <- as.vector(t(Area_corrected[,-1])) 
  Area_IS = df_sub$Area[grep("caffe",df_sub$Name)]
  Height_corrected = as.vector(t(Height_corrected[,-1]))
  Height_IS = df_sub$Height[grep("caffe",df_sub$Name)]
  
  Area_norm <- Area_corrected/Area_IS
  Height_norm <- Height_corrected/Height_IS
  
  data <- data.frame(
    Plate_Martina=Plate_Martina, 
    SampleName=SampleName, 
    Area_norm=as.numeric(Area_norm), 
    Height_norm=as.numeric(Height_norm), 
    GroupName=GroupName,
    Area=Area_corrected,
    Height = Height_corrected,
    Area_IS=Area_IS,
    Height_IS = Height_IS
  )
  data <- data[order(data$SampleName),]
  
  if(sum(!is.na(data$Area_norm[grep("GMM",data$SampleName)])!=0)){
    df.negCtrlpos <- data.frame(SampleName=data$SampleName[grep("GMM",data$SampleName)[!is.na(data$Area_norm[grep("GMM",data$SampleName)])]],Area_norm=data$Area_norm[grep("GMM",data$SampleName)[!is.na(data$Area_norm[grep("GMM",data$SampleName)])]])
    negCtrlpos <- df.negCtrlpos[df.negCtrlpos$Area_norm!=0,]
    print("negative ctrl is positive! check negCtrlpos")  
  }
  return(data)
}


exp12read_raw = function() {
  files.peaks = list(
    LlactisA_StandardA= "F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140517_140509_p2Llactis_p4standard.xlsx",
    BvulgatusA_one = "F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140518_140509_p3Bvulgatus.xlsx",
    StandardB_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140520_140509_p4standard.xlsx",
    EcoliiAi1A_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140529_140527_p1_E.coliiAi1.xlsx",
    LplantarumA_CramosumA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140530_140527_p2_p3_L.plantarum_C.ramosum.xlsx",
    CsaccharoA_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140531_140528_p4_C.saccharolyticum.xlsx",
    CboltaeA_LgasseriA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140601_140528_p5_p6_C.boltae_L.gasseri.xlsx",
    BuniformisA_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140602_140527_p7_B.uniformis.xlsx",
    LparacaseiA_RgnavusA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140603_140530_p1_p2_L.paracasei_R.gnavus.xlsx",
    BlongumA_BfragilisA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/140619_140605_140530_p3_p4_Bifido.longum_B.fragilis.xlsx",
    CsaccharoA2_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140916_140528_p4Csaccharo_splitedUPLC.xlsx",
    BuniformisA2_freezetest="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140917_140528_p7Buniformis_140917_pfreezetest.xlsx",
    EcoliiAi1A2_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140919_140527_p1EcoliiAi1.xlsx",
    
    Fnucleatumaerob_Fnucleatumanaerob="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/150116_141210_p4_Fnucleatumnucleatum_p4anaerob.xlsx",
    Buniformisaerob_Buniformisanaerob="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/150217_141210_p1_Buniformis_p1anaerob.xlsx",
    Cbolteaeaerob_Cbolteaeanaerob="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/150217_141210_p3_Cbolteae_p3anaerob.xlsx",
    Banimalislactisaerob_Banimalislactisanaerob="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/150218_141210_p2_Banimalislactis_p2anaerob.xlsx",
    
    BfragilisB_BthetaB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150724_150722_p1_p2_Bfragilis_Bthetataomicron.xlsx",
    BuniformisB_BvulgatusB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150728_150722_p3_p4_Buniformis_Bvulgatus.xlsx",
    BanimalislactisB_BlonguminfantisB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150805_150722_p5_p6_Banimalislactis07_Blonguminfantis.xlsx",
    BlongumB_CcomesB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150805_150722_p7_p8_Blongumlongum_Ccomes.xlsx",
    EcoliiAi1B_FnucleatumB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150807_150729_p3_p4_EcoliiAi1_Fnucleatum.xlsx",
    CbolteaeC_EcoliED1aA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150810_150729_p1_p2_Cbolteae_EcoliED1a.xlsx",
    LparacaseiB_LplantarumB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150811_150729_p5_p6_Lparacasei_Lplantarum.xlsx",
    RgnavusD_RtorquesB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150813_150729_p7_p8_Rgnavus_Rtorques.xlsx",
    ErectaleA_LlactisB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150816_150805_p1_p2_Erectale_Llactis.xlsx",
    EcoliED1aB_FnucleatumC="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/151012_150817_150805_p3_p4_EcoliED1a_Fnucleatum.xlsx",
    
    coldA_warmA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160216_160212_150911_p3cold_p2warm.xlsx",
    Buniformis66A1_Buniformis65A1="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160222_160215_160212_p1_Buniforis66_p2_Buniformis65.xlsx",
    RtorquesC1_LlactisC1="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160222_160217_160212_p3_Rtorques_p4_Llactis.xlsx",
    ElentaA_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160330_160321_160219_p1_Elenta.xlsx",
    ElentaB_SsalivariusB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160405_160402_160119_p2_p3_Elenta_Ssalivarius.xlsx",
    
    Buniformis65B_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160415_160406_160219_p4_Buniformis65.xlsx",
    #coldB_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160418_160219_p10cold.xlsx",
    #warmB_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/160419_160219_p11warm.xlsx",
    
    MixNoA_MixDegradA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/161112_161024_160219_p7_p8_MixNo_MixDegrad.xlsx",
    RtorquesC2_LlactisC2="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/161115_161025_160212_p3_p4_Rtorques_Llactis.xlsx",
    Buniformis66A2_Buniformis65A2="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/161121_161027_160212_p1_p2_Buniformis66_Buniformis65.xlsx",
    MixDegradB_MixNoB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/161122_161019_160212_p7_p8_MixDegrad_MixNo.xlsx",
    Buniformis66B_ErectaleB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/161122_161020_160219_p5_p6_Buniformis66_Erectale.xlsx"
  )
  
  files.IS = list(
    BlonguminfantisA_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140919_140604_p3Blonguminfantis.xlsx",
    CramosumB_CboltaeB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140920_140822_p1Cramosum_p2Cboltae.xlsx",
    CcomesA_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140921_140822_p8Ccomes.xlsx",
    LgasseriB_CsaccharoB="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140922_140822_p3Lgasseri_p4Csaccharo.xlsx",
    RgnavusB_SsalivariusA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_140923_140822_p5Rgnavus_p6Ssalivarius.xlsx",
    destBvulgatus_destBlongumlongum="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_141006_140924_p3Bvulgatus_p4Blongumlongum_destructio.xlsx",
    destBfragilis_destEcoliiAi1="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_141007_140924_p5Bfragilis_p6EcoliiAi1_destructio.xlsx",
    destLgasseri_destCsaccharo="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141010_141008_140924_p7Lgasseri_p8Csaccharo_destructio.xlsx",
    CramosumA2_CsaccharoB2="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141110_141108_140528_p3Cramosum_p4Csaccharo.xlsx",
    CboltaeA2_LgasseriA2="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141114_141110_140528_p5_C.boltae_140527_p6_L.gasseri.xlsx",
    RtorquesA_FnucleatumA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141111_141104_140604_p4_R.torques_p5_F.nucleatumnucleatum.xlsx",
    BanimalislactisA_BthetaA="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141111_141105_140604_p6_B.animalislactis_p8_B.thetaiotaomicron.xlsx",
    RgnavusC_one="F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141111_141107_140604_p2_R.gnavus.xlsx"
  )
  
  raw.peaks = lapply(names(files.peaks), function(z) peak.read(files.peaks[[z]]) %>% dplyr::mutate(Plate_Martina=z))
  raw.IS = lapply(names(files.IS), function(z) peak.read(files.IS[[z]]) %>% dplyr::mutate(Plate_Martina=z))

  
  code = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/141021_platetospecies_encoding.tsv", na="")
  depletion_plate_map = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/uplc_plata_mapping.tsv")


  # names(raw.peaks) = names(files.peaks)
  # names(raw.IS) = names(files.IS)
  options(warn=1)
  normalized.peaks = lapply(raw.peaks, function(z) {
    zz<<-z
    # z = z %>% dplyr::left_join(IS_factor_map, by="Plate_Martina") %>% dplyr::mutate(IS_Factor=tidyr::replace_na(IS_Factor, 1L))
    peak.common.nomalize(z %>% dplyr::mutate(IS_Factor=1L))
  })
  normalized.IS = lapply(raw.IS, function(z) {
    # z = z %>% dplyr::left_join(IS_factor_map, by="Plate_Martina") %>% dplyr::mutate(IS_Factor=tidyr::replace_na(IS_Factor, 1L))
    peak.common.nomalize(z %>% dplyr::mutate(IS_Factor=ifelse(grepl("2", z$Plate_Martina[1]), 1, 0.5)))
    })
  
  
  # TODO: cleanup
  # dim(do.call(rbind, normalized))
  # 
  # which(sapply(normalized.IS, function(z) any(grepl("CboltaeA2_LgasseriA2|CramosumA2_CsaccharoB2", z$Plate_Martina))))
  # depletion_plate_map %>%
  #   dplyr::filter(grepl("CramosumA2|CsaccharoB2|LgasseriA2|CboltaeA2", Plate)) %>%
  #   dplyr::arrange(Plate_Martina)
  # data.clean.old %>% dplyr::filter(grepl("CramosumA2", Plate) & GroupName=="metf" & Status=="ctrl") %>% dplyr::select(DrugRatio)
  # normalized.IS[[9]] %>% dplyr::filter(grepl("CramosumA2_CsaccharoB2", Plate_Martina) & GroupName=="metf" & grepl("ctrl", SampleName)) %>% dplyr::select(DrugRatio)
  # do.call(rbind, normalized) %>% dplyr::filter(grepl("CramosumA2_CsaccharoB2", Plate_Martina) & GroupName=="metf" & grepl("^1.*ctrl", SampleName))
  
  
  # load("F:/home/andrejev/Workspace/170610_save/computational stuff/Drug Degradation analysis/161122_processingUPLCpeaks_data.norm.long"); data.long.old = data.long
  load(file="F:/home/andrejev/Workspace/170610_save/computational stuff/Drug Degradation analysis/data.norm.161122.RData"); data.norm.old = data.norm
  
  normalized_df = do.call(rbind, append(normalized.peaks, normalized.IS)) %>%
    dplyr::mutate(Plate_Martina_Number=as.numeric(gsub("_.*", "", SampleName))) %>% 
    dplyr::inner_join(depletion_plate_map, by=c("Plate_Martina", "Plate_Martina_Number")) %>%
    dplyr::group_by(Plate, SampleName, GroupName) %>%
    dplyr::mutate(WellReplicate=1:n()) %>%
    dplyr::ungroup()
  
  # do.call(rbind, append(normalized.peaks, normalized.IS)) %>%
  #   dplyr::mutate(Plate_Martina_Number=as.numeric(gsub("_.*", "", SampleName))) %>% 
  #   dplyr::mutate(SampleName=gsub("^[^_]*_", "1_", SampleName)) %>% 
  #   dplyr::filter(Plate_Martina=="BanimalislactisA_BthetaA" & grepl("leva_ctrl", SampleName))
  # 
  # normalized_df %>% 
  #   dplyr::filter(Plate_Martina=="BanimalislactisA_BthetaA" & grepl("leva_ctrl", SampleName))
  
  
  data.export = normalized_df %>%
    dplyr::mutate(Extraction="Supernatant") %>%
    dplyr::mutate(Status=dplyr::case_when(
      grepl("ctrl", SampleName) ~ "ctrl",
      grepl("GMM", SampleName) ~ "GMM",
      T ~ "sample"))
  readr::write_tsv(data.export, "C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp1depletion/peaks_exported.tsv")

  
  # CsaccharoB2, CramosumA1, CboltaeA2 - values are different
  normalized_df.wide = normalized_df %>%
    dplyr::filter(!(Plate %in% c("BuniformisA","destBfragilis","destBvulgatus","destLgasseri","freezetest","destEcoliiAi1","destBlongumlongum","destCsaccharo"))) %>%
    dplyr::mutate(Extraction="Supernatant") %>%
    dplyr::mutate(Plate_Area=paste0(Plate, ".Area")) %>%
    reshape2::dcast(SampleName + Extraction + GroupName + WellReplicate ~ Plate_Area, value.var="Area_norm") #%>%
  
  
  data.get = get.datalong(normalized_df.wide, code)
  data.long = data.get[[1]]
  
  
  # names(normalized) = sapply(normalized, function(z) z$Plate_Martina[1])
  # normalized_df = do.call(cbind.data.frame, normalized)
  # normalized_df = normalized_df[,grep("Area_norm",colnames(normalized_df))]
  # all.names <- sapply(strsplit(colnames(normalized_df),".",fixed=T), function(z) as.character(z[[1]]))
  # all.names <- unlist(strsplit(all.names,"_",fixed=T))
  # 
  # normalized_df.wide <- cbind(normalized_df[1:114,],normalized_df[115:228,])
  # all.names <- paste(c(all.names[seq(1,length(all.names),2)],all.names[seq(2,length(all.names),2)]),".Area", sep="")
  # colnames(normalized_df.wide) <- all.names
  # normalized_df.wide = cbind(normalized[[1]][1:114,c(1,5)], normalized_df.wide) #now with sample and drug info but without height info
  # normalized_df.wide = normalized_df.wide[,!grepl("one",colnames(normalized_df.wide))]
  # normalized_df.wide = normalized_df.wide[, !colnames(normalized_df.wide) %in% c("BuniformisA.Area","destBfragilis.Area","destBvulgatus.Area","destLgasseri.Area","freezetest.Area","destEcoliiAi1.Area","destBlongumlongum.Area","destCsaccharo.Area")]
  # 
  # # bacteria_area = colnames(normalized_df.wide[,3:length(normalized_df.wide)])
  # # drugs.real = c("Acetaminophen","Aripiprazole","Digoxin","Donepezil","Duloxetine","Ezetimibe","Levamisole","Loperamide","Metformin","Metronidazole","Montelukast","Ranitidine","Roflumilast","Rosiglitazone","Rosuvastatin","Simvastatin","Sulfasalazine","Tenofovir","Tolmetin")
  # code = read.csv("F:/home/andrejev/Workspace/170610_save/experiments/UPLC data/141021_platetospecies_encoding.csv",header=T,stringsAsFactors=F)
  # data.get = get.datalong(normalized_df.wide, code)
  # data.long = data.get[[1]]
  
  
  # The growth that is used here is growth that Martina calculated per each plate. 
  # In all the papers the average growth for species/drug replicates is used
  data.growth.tsv = read.table("data/exp0growth/2016-11-28_curves_annotation.tab", stringsAsFactors=F, na.strings="", sep="\t", quote="", header=T)
  load(file="F:/home/andrejev/Workspace/170610_save/experiments/Growth Curves/161123own_data.growth.RData") # data.growth.old = data.growth
  data.growth = data.growth %>%
    dplyr::filter(!grepl("destruction", File)) %>%
    dplyr::mutate(Species=gsub("subsp. ","",Species)) %>%
    dplyr::arrange(Species, File) %>%
    dplyr::group_by(Species) %>%
    dplyr::mutate(Replicate=match(File, unique(File))) %>%
    dplyr::ungroup()
  
  GMMgrowth = data.growth %>%
    dplyr::filter(Status=="GMM" & Well!=3) %>%
    dplyr::distinct(growth, Plate) %>% 
    dplyr::group_by(Plate) %>%
    dplyr::filter(!(n()>1 & growth)) %>%
    dplyr::ungroup() %>%
    tibble::add_row(growth=F, Plate=c("warmA", "coldA")) %>%
    dplyr::rename(growth.GMM="growth")
  
  data.secondctrlwell = data.growth %>%
    dplyr::filter(Status=="ctrl") %>%
    dplyr::mutate(Well=2) %>% 
    dplyr::anti_join(data.growth %>% dplyr::distinct(Plate,GroupName,Status,Well))
  
  data.combined = data.long %>%
    dplyr::left_join(data.growth %>% dplyr::full_join(data.secondctrlwell), by=c("Plate","GroupName","Status","Well")) %>%
    dplyr::inner_join(GMMgrowth, by="Plate") %>%
    dplyr::mutate(growth=ifelse(Status=="GMM", growth.GMM, growth)) %>%
    dplyr::mutate(growth=ifelse(growth, "growth", "no growth")) %>%
    dplyr::select(-growth.GMM)
  
  
  #data.clean contains data corrected for no growth in ctrls, extreme increase of drug concentration(limit to 9fold max)
  # possibly kick out plates where GMM_ctrl didnt grow - eliminate noise from those plates, but reduces ctrl samples
  data.clean.new = get.dataclean(data.combined)[[1]]
  # save(data.clean, file="170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction")
  
  load("data/exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData"); data.clean.old = data.clean
  x1 = data.clean.old %>%
    dplyr::inner_join(data.clean.new, by=c("Plate", "GroupName", "Species.x", "IS", "Status", "Batch"))
  
  p = ggplot(x1) +
    geom_abline(intercept=0, slope=1) +
    geom_point(aes(DrugRatio.x, DrugRatio.y, color=Plate), size=0.1) +
    facet_wrap(~Batch.x) +
    coord_cartesian(xlim=c(0, 4), ylim=c(0, 4))
  plotly::ggplotly(p)
  
  xx2 = data.clean.old %>%
    dplyr::anti_join(data.clean.new, by=c("Plate", "GroupName", "Species.x", "Status"))
  dim(xx2)
  View(xx2)
  
  table(xx2$Status, xx2$drug)
  
  xx1 = data.clean.new %>%
    dplyr::anti_join(data.clean.old, by=c("Plate", "GroupName", "Species.x", "IS", "Status"))
  
  #
  # TESTING
  #
  
  x2 = x1 %>%
    # dplyr::filter(Batch==2 & !grepl("CramosumA2|CsaccharoB2|LgasseriA2|CboltaeA2", Plate)) %>%
    dplyr::filter(DrugRatio.y/DrugRatio.x > 1.5) 
  
  ggplot(x1) +
    geom_density(aes(x=DrugRatio.y))
  
  ggplot(x1 %>% dplyr::filter(DrugRatio.y>1)) +
    geom_density(aes(x=DrugRatio.y/DrugRatio.x, color=!grepl("CramosumA2|CsaccharoB2|LgasseriA2|CboltaeA2", Plate))) +
    geom_vline(xintercept=1:2) +
    coord_cartesian(xlim=c(0,5))
  
  plot(density(x2$DrugRatio.x/x2$DrugRatio.y, na.rm=T))
  
  x2 = x1 %>%
    dplyr::filter(grepl("CramosumA2|CsaccharoB2|LgasseriA2|CboltaeA2", Plate)) %>%
    dplyr::select(GroupName:Plate_order, dplyr::matches("^DrugRatio"), dplyr::matches("^growth"), dplyr::matches("^N_ctrls"), dplyr::matches("^CtrlMean"), dplyr::matches("^NormedCtrl"), dplyr::matches("^BatchMean_ctrl"), dplyr::matches("^DiffToBatchMean"), dplyr::matches("^N_ctrls"))
  
  # Outliers: CramosumA2, LgasseriA2, CboltaeA2
  # CboltaeA2_LgasseriA2 -> LgasseriA2, CboltaeA2
  # CramosumA2_CsaccharoB2 -> CramosumA2, CsaccharoB2
  
  plot(sort(x2$NormedCtrl.x), y=sort(x2$NormedCtrl.y))
  
  ggplot(x2) +
    geom_abline(intercept=0, slope=1) +
    geom_point(aes(NormedCtrl.x, NormedCtrl.y, color=grepl("CramosumA2|CsaccharoB2|LgasseriA2|CboltaeA2", Plate)), size=0.1) +
    facet_wrap(~Batch)
  
  colnames(data.clean.new)
  col = "Well"
  plot(data.clean.new[[col]], data.clean.old[[col]])
  table(data.clean.new[[col]]==data.clean.old[[col]])
  
  data.clean.new %>%
    dplyr::filter(Plate=="CboltaeA2") %>%
    dplyr::filter(Species.x=="Clostridium bolteae") %>%
    dplyr::filter(grepl("dulo|acet|rosi|teno|ezet|metf|rosu", GroupName)) %>%
    dplyr::select(GroupName:Plate_order, dplyr::matches("^DrugRatio"), dplyr::matches("^growth"), dplyr::matches("^N_ctrls"), dplyr::matches("^CtrlMean"), dplyr::matches("^NormedCtrl"), dplyr::matches("^BatchMean_ctrl"), dplyr::matches("^DiffToBatchMean"), dplyr::matches("^N_ctrls"))
  
  data.clean.old %>%
    dplyr::filter(Plate=="CboltaeA2") %>%
    dplyr::filter(Species.x=="Clostridium bolteae") %>%
    dplyr::filter(grepl("dulo|acet|rosi|teno|ezet|metf|rosu", GroupName)) %>%
    dplyr::select(GroupName:Plate_order, dplyr::matches("^DrugRatio"), dplyr::matches("^growth"), dplyr::matches("^N_ctrls"), dplyr::matches("^CtrlMean"), dplyr::matches("^NormedCtrl"), dplyr::matches("^BatchMean_ctrl"), dplyr::matches("^DiffToBatchMean"), dplyr::matches("^N_ctrls"))
  
  ggplot(x2) +
    geom_abline(intercept=0, slope=1) +
    geom_point(aes(DrugRatio.x, DrugRatio.y, color=Plate)) +
    facet_wrap(~Batch)
  
}