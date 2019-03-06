dir.create("reports", showWarnings=F)
library(gplots)
library(MSnbase)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
source("functions.R")

exp3proteomics.analyze = function()
{
  ####
  #### PART 1
  ####
  # Add MNAT (Missing not at random) for pairwise comparison
  #read in the whole dataset
  dataset = read.delim("data/exp3proteomics/Data_all.txt", as.is = TRUE)
  #read in NA table (make sum of NA per sample group in a seperate file)
  na.count = read.delim("data/exp3proteomics/NA_counts.txt", as.is = TRUE)
  # define number of samples total inclusing control
  ns=2
  #define number of replicates per sample
  nrep=4
  p = seq(1,(nrep*ns),by=nrep)
  
  
  ##### Control vs Du  ####
  control = ((grep("Control", colnames(dataset))[1]-1)/nrep)+1
  controls = dataset[,c(p[control]:(p[control]+nrep-1))]
  sample.num = seq(1,ns,by=1)
  sample.num.noctr = sample.num[-control]
  # Compared to Du ----------------------------------------------------
  
  i=1
  
  na.count2 = na.count
  samples = dataset[,(p[sample.num.noctr[i]]:(p[sample.num.noctr[i]]+nrep-1))]
  data = cbind(controls,samples)
  name = colnames(samples)[1] 
  #log2 transform all the intensities
  data1 = log2(data)
  
  # define MNAR (Missing not at random). Must have a 3+ difference in number of NAs in sample in control.
  na.count2$randna = ifelse((na.count[,control] == 0 & na.count[, sample.num.noctr[i]] >= 3) | (na.count[,control] >= 3&na.count[, sample.num.noctr[i]] == 0)| 
                             (na.count[,control] == 1 & na.count[, sample.num.noctr[i]] == 4)|(na.count[,control] == 4&na.count[, sample.num.noctr[i]] == 1),"FALSE", "TRUE")
  
  #removes all proteins which are not sufficiently identified in both samplegroup
  data2 = data1[na.count[,control]<=1|na.count[, sample.num.noctr[i]]<=1,]
  
  #get the same dataframe dimesion as in data
  na.count.1 = na.count2[na.count[,control]<=1|na.count[, sample.num.noctr[i]]<=1,]
  
  # number of proteins that have Missing not at random (MNAR) values (e.g 4/0, 3/0, 0/4, 0/3)
  pdf(paste0("reports/exp3proteomics_", name, "_vs_Control.pdf", sep=""))
  boxplot(data1, main= "Original data distribution")
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
  fdata = na.count.1
  
  write.table(fdata, "data/exp3proteomics/fdataFile.txt", sep="\t", quote=FALSE)
  
  # write out phenotype data file (it must have the same number of rows as number of columns of data and rownames as data)
  
  pdata = data.frame(row.names=colnames(data), sample.name=c("Con_1", "Con_2", "Con_3", "Con_4", "Du_1", "Du_2", "Du_3", "Du_4"))
  write.table(pdata, "data/exp3proteomics/pdataFile.txt", sep="\t", quote=FALSE)
  
  
  # generates "res" in values (this contains all the information from the 3 necessary files (expression data, feature file, phenotype file))
  res = readMSnSet(exprsFile="data/exp3proteomics/exprsFile.txt", featureDataFile="data/exp3proteomics/fdataFile.txt", phenoDataFile="data/exp3proteomics/pdataFile.txt", sep="\t", header=TRUE)
  
  #correlation original data
  correlation = cor(exprs(res), method="pearson",use="pairwise.complete.obs")
  pdf("reports/exp3proteomics_correlation_original_data.pdf", width=11, height=11)
  heatmap.2(correlation,trace = "none",density.info="none", cexRow=0.9, cexCol=0.9, col=RColorBrewer::brewer.pal(11, "PRGn"), 
            colsep=c(seq(1,ncol(correlation),1)),
            rowsep=c(seq(1,ncol(correlation),1)),
            sepcolor='white',sepwidth=c(0.0125,0.02), main="correlation_original_data")
  dev.off()
  
  #data imputation 
  res3  =  MSnbase::impute(res, method = "mixed",  randna = fData(res)$randna, mar = "knn", mnar = "MinDet")
  res.norm = normalise(res3, "quantiles")
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
  final_data = as.data.frame(exprs(res.norm)) %>%
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
      is_significant=ifelse(p.log10>1 & FC>2, 1, 0),
      group=dplyr::case_when(
        is_significant>0 & grepl("^(2.4.2.7|1.17.1.4|6.3.3.1|1.3.1.14|4.1.1.23|2.7.1.11|1.3.)", EC.number) ~ "annotated",
        is_significant>0 ~ "enriched", 
        T ~ "not sign."),
      text=ifelse(group=="annotated", EC.number, ""))
  
  final_data.export = final_data %>% 
    dplyr::mutate(is_significant=ifelse(is_significant>0, "Yes", "No")) %>%
    dplyr::select(
      proteinId, EC.number, is_significant, FC, pvalue, padjst=p.adjst,
      Control_1, Control_2, Control_3, Control_4, Du_1, Du_2, Du_3, Du_4,
      Pathway, Protein.names, GO.Molecular_Function=Gene.ontology..GO., GO.Biological_Process=Gene.ontology..biological.process.
    )
  readr::write_tsv(final_data.export, "reports/exp3proteomics_volcano_plot_data.tsv", col_names=T, na="")  
  
  pdf("reports/exp3proteomics_volcano_plot.pdf", width=11, height=11)
  ggplot(aes(x=FC,y=p.log10, color=group), data=final_data) +
    geom_point() +
    scale_color_manual(values=c("not sign."="darkgrey", "enriched"="#893033", "annotated"="#241011")) +
    geom_hline(yintercept=-log10(0.1), linetype="dotted") +
    geom_vline(xintercept=2, linetype="dotted") +  
    ggrepel::geom_text_repel(aes(label=text), size=8, show.legend=F) +
    scale_x_continuous(limits=c(-9,9), breaks=seq(-8,8,2)) +
    labs(x="log2 fold change of imputated intensity", y="-log10 of FDR adjusted p-value", color="") +
    myTheme
  dev.off()
  
  #
  # GO terms enrichment
  #
  go.terms = readr::read_delim("data/go_terms.tsv", "\t")
  final_data.go = final_data %>%
    reshape2::melt(measure.vars=c("Gene.ontology..GO."), value.name="go") %>%
    dplyr::filter(!is.na(go)) %>%
    tidyr::separate_rows(go, sep="; ") %>%
    dplyr::mutate(go.name=gsub("(.*) \\[GO:.*", "\\1", go), go.id=gsub(".*(GO:.*)\\].*", "\\1", go)) %>%
    dplyr::inner_join(go.terms, by="go.id") %>%
    dplyr::group_by(go.ontology) %>%
    dplyr::mutate(total_hits=sum(!duplicated(proteinId[is_significant==1])), total_n=length(unique(proteinId))) %>%
    dplyr::group_by(go.ontology, go.name, go.id, total_hits, total_n) %>%
    dplyr::do((function(z) {
      z.ret = data.frame(go_hits=sum(z$is_significant==1), go_n=length(z$proteinId))
      f11 = z.ret$go_hits
      f21 = z.ret$go_n - z.ret$go_hits
      f12 = z$total_hits[1] - z.ret$go_hits[1]
      f22 = z$total_n[1] - z.ret$go_n - z$total_hits[1] + z.ret$go_hits
      t = fisher.test(matrix(c(f11,f21,f12,f22), ncol=2), alternative="two.sided")
      z.ret$pval = t$p.value
      z.ret$odds = t$estimate
      z.ret
    })(.)) %>%
    data.frame() %>%
    dplyr::mutate(padj=p.adjust(pval))
  readr::write_tsv(final_data.export, "reports/exp3proteomics_enrichment.tsv", col_names=T, na="")  
  final_data.go %>% dplyr::filter(pval < 0.05)
}
