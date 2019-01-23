library(multtest)
library(ggplot2)
library(plyr)
library(LSD)
library(boot)
library(ade4)
library(DESeq2)
library(vsn)
library(gplots)
library(RColorBrewer)
library(matrixStats)
library(MASS)
library(vegan)
library(locfit)
library(stringr)
library(sva)
library(limma)
library(corpcor)
library(tidyr)
library(reshape2)
library(dplyr)
library(stringr)
library(magrittr)
library(purrr)
library(readr)
library(mutoss)
library(qvalue)
library(xcms)
library(rcdk)
source("functions.R")

#ntop = 500
# "278.2/1089" amitriptyline extra
# "195.1/224" caffeine extra
# "195.1/222" caffeine lys
# "278.2/1093" amitriptyline lys
# 'dulox peak
# 298.1/1050 lys
# 298.1/1040 extra'

log.na = function(x) log10(ifelse(x>0, x, NA))
log.zero = function(x) log10(ifelse(x>0, x, 1))
glog2 = function(x) ((asinh(x)-log(2))/log(2))

fun_uncertaintest = function(data_interest, pair, dulox_peak) {
  calcunc = function(m1,m2,sd1,sd2){ 2*sqrt(((m2/m1)-(m2/(m1+sd1)))^2+((m2/m1)-((m2+sd2)/m1))^2) }
  
  data_pointestimate = data_interest %>%
    dplyr::group_by(Peak) %>%
    dplyr::summarize(Var=var(Intensity_amitrip), SD=sd(Intensity_amitrip), CV=SD/mean(Intensity_amitrip))
  
  # 
  # Between conditions coefficient of variation
  #
  peaks_biovar = data_pointestimate %>% dplyr::filter(CV>0.2) %>% .$Peak
  if(!(dulox_peak %in% peaks_biovar)) peaks_biovar = c(peaks_biovar, dulox_peak)
  
  data_sum = data_interest %>%
    dplyr::filter(Peak %in% peaks_biovar) %>%
    dplyr::group_by(Treatment, Peak) %>%
    dplyr::summarize(Mean=mean(Intensity_amitrip), SD=sd(Intensity_amitrip),CV=sd(Intensity_amitrip)/mean(Intensity_amitrip)) %>%
    data.frame()
  
  data_sum %>%
    dplyr::filter(Treatment==pair[1]) %>%
    dplyr::select(Peak, Mean.ctrl=Mean, SD.ctrl=SD) %>%
    dplyr::inner_join(data_sum %>% dplyr::filter(Treatment==pair[2]) %>% dplyr::select(Peak, Mean.treated=Mean,SD.treated=SD), by="Peak") %>%
      dplyr::filter(SD.treated<0.2 | SD.ctrl<0.2) %>%
      dplyr::group_by(Peak) %>%
      dplyr::mutate(FC=ifelse(Mean.treated>Mean.ctrl,Mean.treated/Mean.ctrl,-Mean.ctrl/Mean.treated), 
      Ufc=ifelse(Mean.treated>Mean.ctrl,calcunc(Mean.treated,Mean.ctrl,SD.treated,SD.ctrl), calcunc(Mean.ctrl,Mean.treated,SD.ctrl,SD.treated)),
      pvalue=t.test2(Mean.ctrl, Mean.treated, SD.ctrl, SD.treated, 6, 6,equal.variance=TRUE)$p_value) %>%
    data.frame() %>%
    dplyr::mutate(
      FCmin=(1/(1-Ufc)), padjust=p.adjust(pvalue, method="BH"),
      logFC=ifelse(FC<0,-log10(abs(FC)),log10(FC))) %>%
    data.frame()
}

load("data/exp3metabolomics/clickxset.fillden.noInj1.RData")
load("data/exp3metabolomics/extraxset.fillden.noInj1.RData")
load("data/exp3metabolomics/lysxset.fillden.noInj1.RData")

KEGG_metabolites = read.csv("data/exp3metabolomics/KEGG_metbolites_info.csv", stringsAsFactors=F)
KEGG_metabolites = read.csv("data/exp3metabolomics/KEGG_metabolites_bugs.csv", stringsAsFactors=F) %>%
  dplyr::select(KEGG.id=METABOLITES, Cs=csh, Bu=bth) %>%
  reshape2::melt(id.vars="KEGG.id", measure.vars=c("Cs", "Bu"), variable.name="species.code") %>%
  dplyr::filter(value>0) %>%
  dplyr::select(-value) %>%
  dplyr::inner_join(KEGG_metabolites, by="KEGG.id") %>%
  dplyr::filter(Monoisotopic.mass != 0) %>%
  dplyr::mutate(H_Monoisotopic.mass=Monoisotopic.mass+1.0078250321, ACN_Monoisotopic.mass=Monoisotopic.mass+42.033823) 

#
# Load experiments data and peaks
#
data_long = data.frame()
peak_long = data.frame()
datasets.all = list("Extracellular"=extraxset.fillden, "Lysate"=lysxset.fillden, "Click"=clickxset.fillden)
for(d in names(datasets.all)) {
  # Load raw data
  data_raw.d = as.data.frame(xcms::groupval(datasets.all[[d]], method="maxint", intensity="into", value="into"))
  
  # Peaks data
  peak_long.d = as.data.frame(datasets.all[[d]]@groups) %>%
    dplyr::mutate(Peak=rownames(data_raw.d)) %>%
    dplyr::mutate(mzlow=mzmed-(mzmed*0.000005), mzhigh=mzmed+(mzmed*0.000005)) %>%
    tidyr::crossing(KEGG_metabolites) %>% 
    dplyr::mutate(donor_ACN=ACN_Monoisotopic.mass>=mzlow & ACN_Monoisotopic.mass<=mzhigh, donor_H=H_Monoisotopic.mass>=mzlow & H_Monoisotopic.mass<=mzmax) %>%
    dplyr::filter(donor_ACN | donor_H) %>%
    dplyr::rename(BLK="Blank") %>%
    reshape2::melt(measure.vars=intersect(c("BLK", "Bu", "BuDulox", "Cs", "CsDulox", "Dulox", "Ctrl", "Heat", "Lys"), names(.)), variable.name="Treatment") %>%
    dplyr::filter(value != 0) %>%
    dplyr::select(-value) %>%
    dplyr::mutate(Dataset=d) %>%
    data.frame()
  peak_long = rbind(peak_long, peak_long.d)
  
  data.d = data_raw.d %>%
    dplyr::mutate(Peak=make.unique(rownames(.))) %>%
    dplyr::mutate_at(dplyr::vars(-dplyr::matches("Peak")), dplyr::funs(round)) %>%
    setNames(gsub("^\\d+_", "", colnames(.))) %>%
    tibble::column_to_rownames("Peak") %>%
    data.matrix()
  
  data_long.d = reshape2::melt(data.d) %>%
    dplyr::select(Peak=Var1, Sample=Var2, Intensity=value) %>%
    tidyr::separate("Sample", into=c("Treatment","State","BioRep","TechRep"), extra="drop", remove=F) %>%
    dplyr::mutate(species.code=ifelse(grepl("(Bu|Cs)", Treatment), gsub(".*(Bu|Cs).*", "\\1", Treatment), NA)) %>%
    dplyr::mutate(Dataset=d, Sample=as.character(Sample), Peak=as.character(Peak), species.code=as.character(species.code)) %>%
    dplyr::select(Dataset, Peak, Sample, Intensity, Treatment, State, species.code, BioRep, TechRep)
  
  data_long = rbind(data_long, data_long.d)
  
  # Mean/SD plot 
  vsn::meanSdPlot(glog2(data.d %>% data.matrix()))
}
peak_long.f = peak_long %>% dplyr::select(Peak, species.code, Treatment, Dataset, Name, mzlow, mzhigh, donor_ACN, donor_H, npeaks)

#
# Normalized
#
data_long.norm = data_long %>%
  dplyr::filter(Treatment!="BLK" & State!="heat") %>%
  dplyr::group_by(Sample, species.code) %>%
  dplyr::do((function(z) {
    denominator_peak = dplyr::case_when(z$Dataset[1]=="Lysate"~"278.2/1093", z$Dataset[1]=="Extracellular"~"278.2/1089")
    denominator = z$Intensity[grepl(denominator_peak, z$Peak)]
    z$IntensityLog10=log10(z$Intensity)
    z$Intensity_amitrip = z$Intensity/denominator
    z$IntensityLog10_amitrip=log10(z$Intensity_amitrip)
    z
  })(.))

#
# Replicates matching table
#
data_long.norm_techrep = data_long.norm %>% 
  dplyr::filter(TechRep=="inj2" & Dataset!="Click") %>%
  dplyr::inner_join(data_long.norm %>% dplyr::filter(TechRep=="inj3"), by=c("Dataset", "Peak","Treatment","State")) %>%
  dplyr::mutate(group="Technical replicates") %>%
  data.frame()
data_long.norm_biorep = data_long.norm %>%
  dplyr::filter(BioRep=="rep2" & Dataset!="Click") %>%
  dplyr::inner_join(data_long.norm %>% dplyr::filter(BioRep=="rep3"), by=c("Dataset", "Peak","Treatment","State")) %>%
  dplyr::mutate(group="Biological replicates") %>%
  data.frame()
data_long.norm_rep = rbind(data_long.norm_techrep, data_long.norm_biorep)

png("reports/exp3proteomics_replicated_correlation.png", width = 2048, height=2048)
data_long.norm_cor = data_long.norm_rep %>% 
  dplyr::group_by(group, Dataset) %>%
  dplyr::summarise(cor=cor(IntensityLog10_amitrip.x, IntensityLog10_amitrip.y, method="spearman", use="pairwise.complete.obs"))
ggplot(data_long.norm_rep, aes(x=IntensityLog10_amitrip.x, y=IntensityLog10_amitrip.y)) +
  geom_point(alpha=0.1,color="black", size=0.1) +
  coord_equal() +
  scale_x_continuous(limits =c(-8,2), breaks=seq(-8,2,2)) +
  scale_y_continuous(limits =c(-8,2), breaks=seq(-8,2,2)) +
  labs(x="Rep. 1, Log10 Peak Intensity normalized by Amitriptyline", y="Rep. 2, Log10 Peak Intensity normalized by Amitriptyline") +
  geom_text(aes(x=-4, y=1, label=round(cor, 3)), data=data_long.norm_cor, colour="black", size=10,inherit.aes=FALSE, parse=FALSE) +
  facet_grid(group~Dataset) +
  myTheme +
  theme(strip.text = element_text(size=48)) +
  theme(axis.text = element_text(size=32)) 
dev.off()


#
# Calculate fold changes of different peaks AUC and corresponding p-values
#
data_amitrip.effect = data_long.norm %>% 
  dplyr::filter(!is.na(species.code)) %>%
  dplyr::group_by(State, Dataset, species.code) %>%
  dplyr::do((function(z) {
    sp = z$species.code[1]
    dulox_peak = ifelse(z$State[1]=="lys", "298.1/1050", "298.1/1040")
    
    z.drug = data_long.norm %>% dplyr::filter(Treatment=="Dulox" & Dataset==z$Dataset[1] & Peak %in% z$Peak) %>% data.frame()
    z.drug = rbind(z.drug, z %>% dplyr::filter(Treatment==paste0(sp, "Dulox")) %>% data.frame())
    ret.drug = fun_uncertaintest(data_interest=z.drug, pair=c("Dulox", paste0(sp, "Dulox")), dulox_peak=dulox_peak) %>% dplyr::mutate(Control="Dulox") %>% data.frame()
    ret.bug = fun_uncertaintest(data_interest=z, pair=c(sp, paste0(sp, "Dulox")), dulox_peak=dulox_peak) %>% dplyr::mutate(Control=sp) %>% data.frame()
    rbind(ret.bug, ret.drug)
  })(.))
data_amitrip.sign_effect = data_amitrip.effect %>% filter(padjust<0.05, abs(FC)>FCmin) %>% data.frame()

#
# 
#
data_amitrip.sign_effect.sum = data_amitrip.sign_effect %>%
  dplyr::filter(Control!="Dulox") %>%
  setNames(gsub("^((?!Dataset|Peak|State|species).*)$", "\\1.bug", names(.), perl=T)) %>%
  dplyr::inner_join(data_amitrip.sign_effect %>%
    dplyr::filter(Control=="Dulox") %>%
    setNames(gsub("^((?!Dataset|Peak|State|species).*)$", "\\1.drug", names(.), perl=T)),
    by=c("Dataset", "State", "Peak", "species.code")) %>%
  dplyr::mutate(FC.both_same=dplyr::between(logFC.bug/logFC.drug, 0.9, 1.1)) %>%
  dplyr::left_join(peak_long.f %>% dplyr::group_by(Dataset, species.code, Peak) %>% dplyr::summarise(MetaboliteNames=paste(unique(Name), collapse=",")), by=c("Dataset", "species.code", "Peak")) %>%  #Treatment.drug=Treatment, 
  dplyr::mutate(Group.detailed=dplyr::case_when(
    species.code=="Bu" & State=="lys" & Peak=="674.5/1499" ~ "extracellular and lysate",
    species.code=="Bu" & State=="lys" & Peak=="298.1/1050" ~ "duloxetine",
    species.code=="Bu" & State=="extra" & Peak=="674.5/1496" ~ "extracellular and lysate",
    species.code=="Bu" & State=="extra" & Peak=="298.1/1040" ~ "duloxetine",
    species.code=="Bu" & State=="extra" & !is.na(MetaboliteNames) ~ "annotated",
    
    species.code=="Cs" & State=="lys" & FC.bug>10 & FC.drug>10 ~ "extracellular and lysate", 
    species.code=="Cs" & State=="lys" & Peak=="298.1/1050" ~ "duloxetine",
    species.code=="Cs" & State=="extra" & FC.bug>10 & FC.drug>10 ~ "extracellular and lysate", 
    species.code=="Cs" & State=="extra" & Peak=="298.1/1040" ~ "duloxetine",
    species.code=="Cs" & State=="extra" & !is.na(MetaboliteNames) ~ "annotated",
    T ~ "unique"
  ), Group=ifelse(Group.detailed=="duloxetine", Group.detailed, "other"))

  # TODO: Why some specific?
  # result_both.Cs.extra <- result_both.Cs.extra[c(4,5,13,15,22,25,30,38,39,47,48,52,54,65,71,72,79,82,86,88,89,92,95,97,106,120)]
  # result_both.Cs.lys <- result_both.Cs.lys[c(1,2,4,5,8,10,11,14,15,18:24,26:32,34,37,41)]

  pdf(file="reports/exp3metabolomics_hits_scatterplot.pdf", paper="a4r")
  ggplot(data_amitrip.sign_effect.sum, aes(x=logFC.bug,y=logFC.drug, size=Mean.treated.bug, color=Group)) +
    geom_point() +
    #geom_text(aes(label=MetaboliteNames), data=data_amitrip.sign_effect.sum %>% dplyr::filter(!is.na(MetaboliteNames))) %>%
    scale_color_manual(values=c(other="darkgrey", unique="darkgrey", duloxetine="dodgerblue4", "extracellular and lysate"="firebrick", annotated="chartreuse4")) +
    scale_x_continuous(limits =c(-2,5), breaks=seq(-2,5,1)) +
    scale_y_continuous(limits =c(-2,5), breaks=seq(-2,5,1)) +
    coord_equal() +
    facet_grid(species.code~State) +
    labs(y="Fold change drug treated bacteria vs. drug control; log10", 
         x="Fold change drug treated bacteria vs. bacteria control; log10", size="Intensity in Treated", color="") +
    myTheme
  dev.off()
  