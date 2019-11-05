dir.create("reports", showWarnings=F)
library(xlsx)
library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)
library(naniar)
source("functions.R")

peak.read = function(path,len, remove.empty=T)
{
  df = read.xlsx(path,1,startRow=1,stringsAsFactors=F)
  df$group = substr(df$SampleName,1,6)
  df = df[order(df$Name,df$Channel.Description,df$Sample.Set.Name,df$SampleName,df$Vial),]
  if (length(unique(df$Sample.Set.Name))>1){
    print("Data processed, but more than 1 sample set detected")
    return(df[1:len,])
  }
  else {
    return(df[1:len,])
  }
}

exp6transfers.analyze = function()
{
  #
  # Read transfer assay targeted metabolomics data
  #
  data.transfer = peak.read("data/exp6transfers/160916_160907_160823_transferassaydulox.xlsx", 272)
  data.transfer.norm = data.transfer %>%
    dplyr::filter(grepl("peak", Name)) %>%
    dplyr::mutate_at(.vars=dplyr::vars(Area, Height), .funs=dplyr::funs(ifelse(is.na(Area), 0, .))) %>%
    dplyr::mutate(Area=as.double(Area), Height=as.double(Height), Date=as.numeric(factor(Date, sort(unique(Date))))) %>% 
    dplyr::group_by(Date) %>%
    dplyr::mutate(MeanCtrl=mean(Area[Condition=="Control"]), Normed=Area/MeanCtrl) %>% 
    data.frame()
  
  data.T = data.transfer.norm %>%
    dplyr::group_by(Date) %>%
    dplyr::mutate(
      pvalue.Bu=t.test(Normed[Condition=="Control"], Normed[Condition=="+Bu;+D"])$p.value,
      pvalue.No=t.test(Normed[Condition=="Control"], Normed[Condition=="-Bu;+D"])$p.value,
      pvalue.BuNo=t.test(Normed[Condition=="+Bu;+D"], Normed[Condition=="-Bu;+D"])$p.value,
      Sign=rowSums(cbind(pvalue.Bu, pvalue.No, pvalue.BuNo)<0.05))
  
  pdf("reports/exp6transfers_duloxetine_degradation.pdf", width=8, height=5)
  data.plot = data.T %>% 
    dplyr::filter(!(Condition %in% c("+Bu;-D","-Bu;-D","d"))) %>%
    dplyr::mutate(ConditionName=dplyr::case_when(
      Condition=="+Bu;+D" ~ "with B. uniformis",
      Condition=="-Bu;+D" ~ "no B. uniformis",
      Condition=="Control" ~ "Bacteria-free control",
      T ~ "This should not happen"
    ))
  ggplot(data.plot %>% dplyr::filter(Condition=="-Bu;+D")) +
    geom_hline(yintercept = 1, linetype="longdash" ) +
    geom_boxplot(aes(x=as.factor(Date), y=Normed, fill=ConditionName)) +
    ylim(c(0,1.1)) +
    scale_fill_manual(values=c("Bacteria-free control"="red", "no B. uniformis"="chartreuse4", "with B. uniformis"="darkorchid4")) + 
    myTheme +
    theme(axis.ticks.length = unit(0, "lines")) +
    labs(x="Transfer", y= "AUC normalized by mean of controls", color="", fill="Community") +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18))
  dev.off()
  
  #
  # Read transfer assay metagenomics data
  #
  map = readr::read_delim("data/exp6transfers/Sample_map_mod.csv", ",")
  data.seq = readr::read_delim("data/exp6transfers/otu_table_wo_mono_normalized_by_gene_copy_mod.csv", ",") %>%
    dplyr::inner_join(map, by=c("Sample")) %>%
    reshape2::melt(measure.vars=c("B. thetaiotaomicron", "B. uniformis", "E. rectale", "L. gasseri", "R. torques", "S. salivarius"), variable.name="Species", value.name="Percent") %>%
    dplyr::mutate(Bug=gsub("([^;]+);.*", "\\1", Condition), Drug=gsub("[^;]+;([^;]+)", "\\1", Condition))
  
  data.seq_plot = data.seq %>%
    dplyr::filter(Condition %in% c("+Bu;-D","+Bu;+D","-Bu;+D","-Bu;-D")) %>%
    dplyr::filter(!(Bug=="-Bu" & Species=="B. uniformis")) %>%
    dplyr::mutate(Drug=ifelse(Drug=="-D", "- Duloxetine", "+ Duloxetine")) %>%
    dplyr::mutate(Bug=ifelse(Bug=="-Bu", "- B. uniformis", "+ B. uniformis"))  %>%
    dplyr::group_by(Condition, Drug, Bug, Species, Date) %>%
    dplyr::summarise(Percent=mean(Percent)) %>%
    dplyr::group_by(Drug, Bug, Species) %>%
    dplyr::mutate(ODCorrection=0.2/Percent[Date==0]) %>%
    dplyr::group_by(Drug, Bug, Date) %>%
    dplyr::mutate(PercentNorm=100*ODCorrection*Percent/sum(ODCorrection*Percent))
  
  species.colors = c("B. thetaiotaomicron"="#4f91c9", "E. rectale"="#cf1f2b", "L. gasseri"="#6c63b8", "R. torques"="#f57a66", "S. salivarius"="#d4d1e7", "B. uniformis"="#6CBA6F")
  
  pdf("reports/exp6transfers_abundance.pdf", width=12, height=8)
  ggplot(data.seq_plot) +
    geom_bar(aes(x=Date, y=Percent, fill=Species), position="stack", stat="identity") +
    scale_fill_manual(values=species.colors) +
    facet_grid(Bug ~ Drug) +
    labs(y="Species abundance (%)", x="Transfers") +
    myTheme
  
  ggplot(data.seq_plot) +
    geom_bar(aes(x=Date, y=PercentNorm, fill=Species), position="stack", stat="identity") +
    scale_fill_manual(values=species.colors) +
    facet_grid(Bug ~ Drug) +
    labs(y="Species abundance normalized to equal inoculum OD (%)", x="Transfers") +
    myTheme
  dev.off()
  
  pdf("reports/exp6transfers_rectale-Bu.pdf", width=9, height=6)
  classicTheme =  theme_classic(base_size=18) + 
    theme(
      strip.text.x = element_text(size=24),
      strip.background = element_blank(),
      legend.position="bottom", legend.direction="vertical", 
      aspect.ratio=1)
  data.seq_plot.rectale = data.seq_plot %>% dplyr::filter(Bug=="- B. uniformis" & Species=="E. rectale" & Date>0)
  gridExtra::grid.arrange(
    ggplot(data.seq_plot.rectale) +
      geom_point(aes(x=Date, y=Percent, color=Drug)) +
      geom_line(aes(x=Date, y=Percent, color=Drug), linetype="dotted") +
      labs(y="\n\nSpecies abundance (%)", x="Transfers") +
      classicTheme,
    ggplot(data.seq_plot.rectale) +
      geom_point(aes(x=Date, y=PercentNorm, color=Drug)) +
      geom_line(aes(x=Date, y=PercentNorm, color=Drug), linetype="dotted") +
      labs(y="Species abundance\nnormalizedto\nequal inoculum OD (%)", x="Transfers") +
      classicTheme,
    ncol=2
  )
  dev.off()
  
  
  
  
  depleters = c("B. uniformis","B. thetaiotaomicron","S. salivarius")
  nopes = c("E. rectale","R. torques","L. gasseri")
  data.deplete = data.seq %>%
    dplyr::filter(Drug=="+D") %>%
    dplyr::group_by(Date,Condition,Species) %>%
    dplyr::mutate(Rep=1:length(Date)) %>%
    dplyr::group_by(Date,Bug,Rep) %>%
    dplyr::summarise(
      meanDep = sum(Percent[which(Species %in% depleters)]), 
      meanNop=sum(Percent[which(Species %in% nopes)]),
      Condition=dplyr::case_when(Bug[1]=="+Bu"~"+Bu;+D", T~"-Bu;+D"),
      meanDep=meanDep/100) %>%
    data.frame()
  
  pdf("reports/exp6transfers_duloxetine_degradation_bacteria_fraction.pdf", width=8, height=5)
  data.deplete.plot = data.deplete %>%
    dplyr::mutate(ConditionName=dplyr::case_when(
      Condition=="+Bu;+D" ~ "with B. uniformis",
      Condition=="-Bu;+D" ~ "no B. uniformis",
      Condition=="Control" ~ "Bacteria-free control",
      T ~ "This should not happen"
    ))
  ggplot(data.deplete.plot %>% dplyr::filter(Condition=="-Bu;+D"), aes(x=as.factor(Date), y=meanDep, fill=ConditionName)) +
    geom_boxplot() +
    ylim(c(0,1.1)) +
    scale_fill_manual(values=c("Bacteria-free control"="red", "no B. uniformis"="chartreuse4", "with B. uniformis"="darkorchid4"))  +
    myTheme +
    labs(x="Transfer", y= "Fraction of duloxetine depleting bacteria", fill="Community")
  dev.off()
}
