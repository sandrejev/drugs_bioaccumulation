library(Mfuzz)
library(limma)
library(gplots)
library(RColorBrewer)
library(MSnbase)
library(vsn)
library(impute)
library(imputeLCMD)
library(ggplot2)
library(stringi)
library(stringr)
source("functions.R")


####
#### PART 1
####
# Add MNAT (Missing not at random) for pairwise comparison
#read in the whole dataset
dataset<-read.delim("data/exp3proteomics/Data_all.txt", as.is = TRUE)
#read in NA table (make sum of NA per sample group in a seperate file)
na.count<-read.delim("data/exp3proteomics/NA_counts.txt", as.is = TRUE)
# define number of samples total inclusing control
ns=2
#define number of replicates per sample
nrep=4
p<-seq(1,(nrep*ns),by=nrep)


##### Control vs Du  ####
control<-((grep("Control", colnames(dataset))[1]-1)/nrep)+1
controls<-dataset[,c(p[control]:(p[control]+nrep-1))]
sample.num<-seq(1,ns,by=1)
sample.num.noctr<-sample.num[-control]
# Compared to Du ----------------------------------------------------

i=1

na.count2<-na.count
samples<-dataset[,(p[sample.num.noctr[i]]:(p[sample.num.noctr[i]]+nrep-1))]
data<-cbind(controls,samples)
name<-colnames(samples)[1] 
pdf(paste0("exp3proteomics_", name, "_vs_Control.pdf", sep=""))
#log2 transform all the intensities
data1<-log2(data)
boxplot(data1, main= "Original data distribution")

# define MNAR (Missing not at random). Must have a 3+ difference in number of NAs in sample in control.
na.count2$randna<-ifelse((na.count[,control] == 0 & na.count[, sample.num.noctr[i]] >= 3) | (na.count[,control] >= 3&na.count[, sample.num.noctr[i]] == 0)| 
                           (na.count[,control] == 1 & na.count[, sample.num.noctr[i]] == 4)|(na.count[,control] == 4&na.count[, sample.num.noctr[i]] == 1),"FALSE", "TRUE")

#removes all proteins which are not sufficiently identified in both samplegroup
data2<-data1[na.count[,control]<=1|na.count[, sample.num.noctr[i]]<=1,]

#get the same dataframe dimesion as in data
na.count.1<-na.count2[na.count[,control]<=1|na.count[, sample.num.noctr[i]]<=1,]

# number of proteins that have Missing not at random (MNAR) values (e.g 4/0, 3/0, 0/4, 0/3)
barplot(c(
  nrow(na.count.1),
  nrow(na.count.1[na.count.1$randna=="TRUE"&(na.count.1[,control] >= 1 & na.count.1[, sample.num.noctr[i]] >= 1),]), 
  nrow(na.count.1[na.count.1$randna=="FALSE",])), 
  names.arg=c("TOTAL","MAR","MNAR"), main="Number of proteins with missing values")
dev.off()
# comparison with "usual" NA removal scheme that keep only proteins with one missing value per sample group, 
# you can see how many proteins you will gain by adding missing values (MNAR)

# write out expression data file
write.table(data2, "data/exp3proteomics/exprsFile.txt", sep="\t", quote=FALSE)

# write out feature data file (it must have the same number of rows and rownames as data)
fdata<-na.count.1

write.table(fdata, "data/exp3proteomics/fdataFile.txt", sep="\t", quote=FALSE)

# write out phenotype data file (it must have the same number of rows as number of columns of data and rownames as data)

pdata<-data.frame(row.names=colnames(data), sample.name=c("Con_1", "Con_2", "Con_3", "Con_4", "Du_1", "Du_2", "Du_3", "Du_4"))
write.table(pdata, "data/exp3proteomics/pdataFile.txt", sep="\t", quote=FALSE)


# generates "res" in values (this contains all the information from the 3 necessary files (expression data, feature file, phenotype file))
res<-readMSnSet(exprsFile="data/exp3proteomics/exprsFile.txt", featureDataFile="data/exp3proteomics/fdataFile.txt", phenoDataFile="data/exp3proteomics/pdataFile.txt", sep="\t", header=TRUE)

#correlation original data
correlation = cor(exprs(res), method="pearson",use="pairwise.complete.obs")
pdf("reports/exp3proteomics_correlation_original_data.pdf", width=11, height=11)
heatmap.2(correlation,trace = "none",density.info="none", cexRow=0.9, cexCol=0.9, col=brewer.pal(11, "PRGn"), 
          colsep=c(seq(1,ncol(correlation),1)),
          rowsep=c(seq(1,ncol(correlation),1)),
          sepcolor='white',sepwidth=c(0.0125,0.02), main="correlation_original_data")
dev.off()

#data imputation 
res3 <- impute(res, method = "mixed",  randna = fData(res)$randna, mar = "knn", mnar = "MinDet")
res.norm<-normalise(res3, "quantiles")
boxplot(exprs(res.norm), main= "Data distribution after normalization")

#correlation
correlation = cor(exprs(res.norm), method="pearson",use="pairwise.complete.obs")
pdf("reports/exp3proteomics_correlation_after_imputation_normalization.pdf", width=11, height=11)
heatmap.2(correlation,trace = "none",density.info="none", cexRow=0.9, cexCol=0.9, col=brewer.pal(11, "PRGn"), 
          colsep=c(seq(1,ncol(correlation),1)),
          rowsep=c(seq(1,ncol(correlation),1)),
          sepcolor='white',sepwidth=c(0.0125,0.02), main = "correlation_after_imputation_normalization")
dev.off()

# extract matrix normalized and with data imputed
write.table(exprs(res.norm), paste("data/exp3proteomics/", name, "_vs_Du_imputation_quan_nor.txt", sep=""), sep="\t", quote=FALSE)

uniprot = readr::read_delim("data/exp3proteomics/uniprot_C_saccharolyticum_enzymes.tab", "\t")
new.data = as.data.frame(exprs(res.norm)) %>%
  dplyr::mutate(proteinId = rownames(.)) %>%
  dplyr::left_join(uniprot, by=c("proteinId"="Entry")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    pvalue=t.test(c(Control_1, Control_2, Control_3, Control_4), c(Du_1, Du_2, Du_3, Du_4))$p.value,
    FC=log2(2^mean(c(Du_1, Du_2, Du_3, Du_4))/2^mean(c(Control_1, Control_2, Control_3, Control_4)))) %>%
  data.frame() %>%
  dplyr::mutate(
    p.adjst=p.adjust(pvalue, method="BH"), 
    p.log10=-log10(p.adjst), 
    group=dplyr::case_when(p.log10>1 & FC>2 ~ "enriched", T ~ "not sign."),
    text=dplyr::case_when(group=="enriched" & grepl("^(2.4.2.7|1.17.1.4|6.3.3.1|1.3.1.14|4.1.1.23|2.7.1.11|1.3.)", EC.number)~EC.number, T~""))


pdf("reports/exp3proteomics_volcano_plot.pdf", width=11, height=11)
ggplot(new.data, aes(x=FC,y=p.log10, color=group)) +
  geom_point() +
  scale_color_manual(values=c("not sign."="darkgrey", "enriched"="#893033")) +
  geom_hline(yintercept=-log10(0.1), linetype="dotted") +
  geom_vline(xintercept=2, linetype="dotted") +  
  ggrepel::geom_text_repel(aes(label=text), size=8, show.legend=F) +
  scale_x_continuous(limits=c(-9,9), breaks=seq(-8,8,2)) +
  labs(x="log2 fold change of imputated intensity", y="-log10 of FDR adjusted p-value", color="") +
  myTheme
dev.off()
#  scale_color_manual(values=c("darkgrey","black","dodgerblue4","firebrick","chartreuse4","darkorchid4"))
