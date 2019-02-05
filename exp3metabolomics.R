library(multtest)
library(ggplot2)
library(plyr)
library(LSD)
library(boot)
library(ade4)
library(DESeq2)
library(vsn)
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

#
# Load KEGG metabolites from Cs/Bu species
#
KEGG_pathways = read.csv("data/exp3metabolomics/KEGG_pathways.csv", stringsAsFactors=F, na.strings=c("")) %>% 
  reshape2::melt(id.vars=c("KEGGid", "PathwayName")) %>% 
  dplyr::filter(!is.na(value)) %>%
  dplyr::select(pathway.id=KEGGid, pathway.name=PathwayName, -variable, pathway.compound=value)
KEGG_metabolites_info = read.csv("data/exp3metabolomics/KEGG_metbolites_info.csv", stringsAsFactors=F)
KEGG_metabolites = read.csv("data/exp3metabolomics/KEGG_metabolites_bugs.csv", stringsAsFactors=F) %>%
  dplyr::select(KEGG.id=METABOLITES, Cs=csh, Bu=bth) %>%
  reshape2::melt(id.vars="KEGG.id", measure.vars=c("Cs", "Bu"), variable.name="species.code", value.name="KEGG.present") %>%
  dplyr::mutate(KEGG.present=KEGG.present>0) %>%
  #dplyr::filter(KEGG.present>0) %>%
  dplyr::inner_join(KEGG_metabolites_info, by="KEGG.id") %>%
  dplyr::filter(Monoisotopic.mass != 0) %>%
  dplyr::mutate(H_Monoisotopic.mass=Monoisotopic.mass+1.0078250321, ACN_Monoisotopic.mass=Monoisotopic.mass+42.033823) 
KEGG_pathways.met_count = KEGG_metabolites %>% 
  dplyr::inner_join(KEGG_pathways, by=c("KEGG.id"="pathway.compound")) %>%
  dplyr::group_by(species.code) %>%
  dplyr::mutate(species.size=length(unique(KEGG.id[KEGG.present]))) %>%
  dplyr::group_by(species.code, pathway.id, species.size) %>%
  dplyr::summarise(pathway.size=length(unique(KEGG.id[KEGG.present])))


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
  if(d %in% c("Extracellular", "Lysate")) {
    peak_long.d = as.data.frame(datasets.all[[d]]@groups) %>%
      dplyr::mutate(Peak=rownames(data_raw.d), mzmed.se=0.000005, mzlow=mzmed*(1-mzmed.se), mzhigh=mzmed*(1+mzmed.se)) %>%
      tidyr::crossing(KEGG_metabolites) %>% 
      dplyr::mutate(donor_ACN=ACN_Monoisotopic.mass>=mzlow & ACN_Monoisotopic.mass<=mzhigh, donor_H=H_Monoisotopic.mass>=mzlow & H_Monoisotopic.mass<=mzhigh) %>%
      dplyr::filter(donor_ACN | donor_H) %>%
      dplyr::rename(BLK="Blank") %>%
      reshape2::melt(measure.vars=intersect(c("Bu", "BuDulox", "Cs", "CsDulox", "Blank", "Ctrl", "Heat", "Lys"), names(.)), variable.name="Treatment") %>%
      dplyr::filter(gsub("(Cs|Bu).*", "\\1", Treatment)==species.code) %>% # value != 0 (we don't consider peaks annotated by CellZome, as we calculate our own mapping)
      dplyr::select(-value) %>%
      dplyr::mutate(Dataset=d) %>%
      data.frame()
    peak_long = rbind(peak_long, peak_long.d)
  }
  
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
  # vsn::meanSdPlot(glog2(data.d %>% data.matrix()))
}

#
# Generate single name for each peak
#
peak_long.sum = peak_long %>% 
  dplyr::left_join(KEGG_pathways, by=c("KEGG.id"="pathway.compound")) %>%
  dplyr::mutate(is_phosphorylated=grepl("phosphate", Name)) %>%
  dplyr::arrange(is_phosphorylated, nchar(Name)) %>% 
  dplyr::group_by(Dataset, species.code, Peak) %>% 
  dplyr::summarise(
    pathways=paste(unique(pathway.id), collapse=","),
    metabolite_display=Name[1], 
    metabolite_kegg=paste(unique(KEGG.id), collapse=","), 
    metabolite_name=ifelse(any(!is_phosphorylated), paste(unique(Name[!is_phosphorylated]), collapse=","), paste(unique(Name), collapse=",")))

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
  })(.)) %>%
  dplyr::left_join(peak_long.sum, by=c("Dataset", "species.code", "Peak"))

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
  })(.)) %>%
  data.frame() %>%
  dplyr::left_join(peak_long.sum, by=c("Dataset", "species.code", "Peak")) %>%
  dplyr::mutate(
    State=factor(State, c("lys", "extra")),
    is_significant=ifelse(padjust<0.05 & abs(FC)>FCmin, 1, 0),
    metabolite_display=dplyr::case_when(State=="lys"&Peak=="298.1/1050" |  State=="extra"&Peak=="298.1/1040" ~ "Duloxetine", T ~ metabolite_display)) %>% 
    dplyr::select(-Ufc, -FCmin)


#
# Significally effected peaks
#
data_amitrip.effect.sum = data_amitrip.effect %>%
  dplyr::filter(Control!="Dulox") %>% setNames(gsub("^((?!Dataset|Peak|State|species.code|metabolite_name|pathways|metabolite.*).*)$", "\\1.bug", names(.), perl=T)) %>%
  dplyr::inner_join(data_amitrip.effect %>% dplyr::filter(Control=="Dulox") %>% setNames(gsub("^((?!Dataset|Peak|State|species.code|pathways|metabolite.*).*)$", "\\1.drug", names(.), perl=T)), by=c("Dataset", "State", "Peak", "pathways", colnames(data_amitrip.effect)[grepl("metabolite", colnames(data_amitrip.effect))], "species.code")) %>%
  dplyr::mutate(
    FC.both_same_33=dplyr::between(logFC.bug/logFC.drug, 0.75, 1.33),
    Group=dplyr::case_when(grepl("Duloxetine", metabolite_display) ~ "Duloxetine", T ~ "Other")
  ) %>%
  dplyr::group_by(species.code, State) %>%
  dplyr::mutate(metabolite_display.top10=ifelse(is_significant.bug & is_significant.drug & !is.na(metabolite_display) & FC.both_same_33 & FC.bug*FC.drug >= sort((FC.bug*FC.drug)[is_significant.bug & is_significant.drug & !is.na(metabolite_display) & FC.both_same_33], decreasing=T)[10], metabolite_display, NA)) %>%
  data.frame() %>%
  dplyr::mutate(Group=ifelse(!is.na(metabolite_display.top10), "Annotated", Group)) %>% 
  data.frame()

#
# Perform hypergeometric test to calculate pathway significance
#
data_amitrip.pathway_hits = data_amitrip.effect.sum %>%
  dplyr::filter(Dataset=="Extracellular") %>% # & Treatment %in% c("Cs", "Bu"), data_long.norm
  tidyr::separate_rows(pathways, sep=",") %>%
  tidyr::separate_rows(metabolite_kegg, sep=",") %>%
  
  # Leave only those pathway/metabolite pairs which are present in KEGG database. The reason is that both pathways and metabolites are lists and not all metabolites come from all pathways
  dplyr::inner_join(KEGG_pathways, by=c("metabolite_kegg"="pathway.compound", "pathways"="pathway.id")) %>%
  dplyr::filter(!is.na(pathways) & pathways != "NA") %>%
  dplyr::inner_join(KEGG_pathways.met_count, by=c("pathways"="pathway.id", "species.code")) %>%

  # Calculate number significant/total metabolites in species
  dplyr::group_by(species.code) %>%
  dplyr::mutate(pathway.id=pathways, species.hits=length(unique(metabolite_kegg[is_significant.bug & is_significant.drug])), species.detected_size=length(unique(metabolite_kegg))) %>%
  
  # Calculate number significant metabolites in pathways
  dplyr::group_by(species.code, pathway.id) %>%
  dplyr::mutate(pathway.detected_size=length(unique(metabolite_kegg)), pathway.hits=length(unique(metabolite_kegg[is_significant.bug & is_significant.drug]))) %>%
  
  # Perform hypergeometric/fisher test  
  dplyr::group_by(species.code, pathway.id, pathway.name, pathway.hits, pathway.size, pathway.detected_size, species.detected_size, species.size, species.hits) %>% 
  dplyr::summarise(pvalue.hyper=phyper(q=pathway.hits[1], m=pathway.size[1], n=species.size[1]-pathway.size[1], k=species.hits[1], lower.tail=F)) %>%
  data.frame() %>%
  
  # Multiple testing adjustment
  dplyr::filter(pathway.detected_size > 1) %>%
  #dplyr::group_by(species.code) %>%
  dplyr::mutate(padjst.hyper=p.adjust(pvalue.hyper), is_significant=ifelse(padjst.hyper<0.05, "Y", "N")) %>% 
  dplyr::arrange(species.code, dplyr::desc(pathway.size), dplyr::desc(pathway.hits/pathway.size))

# Show hits
# data_amitrip.pathway_hits %>% 
#   dplyr::filter(species.code=="Cs" & pvalue.hyper<0.05) %>%
#   dplyr::mutate(sign=ifelse(padjst.hyper<=0.05, "Y", "")) %>%
#   dplyr::select(pathway.name, pvalue.hyper, sign, pathway.hits, pathway.detected_size, pathway.size, species.hits, species.size, species.detected_size, pvalue.fisher, f11, f21, f12, f22) %>% 
#   data.frame() %>%
#   View()


data_amitrip.effect.sum_export = data_amitrip.effect.sum %>% 
  dplyr::mutate(Species=dplyr::case_when(species.code=="Sc"~"Clostridium saccharolyticum", species.code=="Bu"~"Bacteroides uniformis", T~"No species")) %>%
  dplyr::select(
    Dataset, Species, Peak, metabolite_name.display=metabolite_display, metabolite_name.all=metabolite_name, metabolite_kegg,
    
    is_significant.treated_vs_species=is_significant.bug,
    is_significant.treated_vs_dulox=is_significant.drug,
    FC.treated_vs_species=FC.bug,
    FC.treated_vs_dulox=FC.drug,
    
    pvalue.treated_vs_specie=pvalue.bug,
    padjust.treated_vs_specie=padjust.bug,
    pvalue.treated_vs_dulox=pvalue.drug,
    padjust.treated_vs_dulox=padjust.drug,
    
    Mean.treated_species=Mean.treated.bug,
    SD.treated_species=SD.treated.bug,
    Mean.dulox_control=Mean.ctrl.drug,
    SD.dulox_control=SD.ctrl.drug,
    Mean.species_control=Mean.ctrl.bug,
    SD.species_control=SD.ctrl.bug
  )

  readr::write_tsv(data_amitrip.effect.sum_export, "reports/exp3metabolomics_hits_scatterplot_data.tsv", col_names=T, na="")  

  pdf(file="reports/exp3metabolomics_hits_scatterplot.pdf", width=11, height=11)
  ggplot(aes(x=logFC.bug, y=logFC.drug, size=Mean.treated.bug), data=data_amitrip.effect.sum %>% dplyr::filter(Group=="Other" & is_significant.bug & is_significant.drug)) +
    geom_point(aes(color=Group)) +
    geom_point(aes(color=Group), data=data_amitrip.effect.sum %>% dplyr::filter(Group!="Other" & is_significant.bug & is_significant.drug)) +
    ggrepel::geom_text_repel(aes(label=metabolite_display.top10), data=data_amitrip.effect.sum %>% dplyr::filter(!is.na(metabolite_display.top10) & is_significant.bug & is_significant.drug), size=4) +
    scale_color_manual(values=c(Other="darkgrey", Annotated="#000000", Duloxetine="#8A3134")) +
    scale_x_continuous(limits =c(-2,5), breaks=seq(-2,5,1)) +
    scale_y_continuous(limits =c(-2,5), breaks=seq(-2,5,1)) +
    scale_size_continuous(breaks=c(0.4, 1.6)) +
    coord_equal() +
    facet_grid(species.code~State) +
    labs(y="Fold change drug treated bacteria vs. drug control; log10", 
         x="Fold change drug treated bacteria vs. bacteria control; log10", size="Intensity in Treated", color="") +
    myTheme
  dev.off()
  
  
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
  