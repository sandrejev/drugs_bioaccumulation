dir.create("reports", showWarnings=F)
library(ggplot2)
library(vsn)
library(vegan)
library(stringr)
library(tidyr)
library(reshape2)
library(dplyr)
library(stringr)
library(readr)
library(xcms)
source("functions.R")

#ntop = 500
# "278.2/1089" amitriptyline extra
# "195.1/224" caffeine extra
# "195.1/222" caffeine lys
# "278.2/1093" amitriptyline lys
# 'dulox peak
# 298.1/1050 lys
# 298.1/1040 extra'

fun_uncertaintest = function(data_interest, pair, dulox_peak) {
  #
  # K Ortmayr et al. Uncertainty budgeting in fold change determination and implications for non-targeted metabolomics studies in model systems
  # https://pubs.rsc.org/en/content/articlehtml/2016/an/c6an01342b
  #
  # Vinaixa, M. et al., 2012. A Guideline to Univariate Statistical Analysis for LC/MS-Based Untargeted Metabolomics-Derived Data. Metabolites, 2(4), pp.775
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


exp3metabolomics.analyze = function()
{
  load("data/exp3metabolomics/clickxset.fillden.noInj1.RData")
  load("data/exp3metabolomics/extraxset.fillden.noInj1.RData")
  load("data/exp3metabolomics/lysxset.fillden.noInj1.RData")
  
  bugs = readr::read_delim("data/bug_map.tsv", "\t")
  
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
    data_raw.xcms = xcms::groupval(datasets.all[[d]], method="maxint", intensity="into", value="into")
    data_raw.d = data.frame(data_raw.xcms, check.names=F) %>%
      dplyr::mutate(Peak=rownames(data_raw.xcms)) %>%
      data.frame(check.names=F)
    
    
    # Peaks data
    if(d %in% c("Extracellular", "Lysate")) {
      # This is used to prefilter before interval joining with dplyr because crossing with interval joining is very slow
      mzmed.se = 0.000005
      ACN.lb = KEGG_metabolites$ACN_Monoisotopic.mass*(1-mzmed.se)
      ACN.ub = KEGG_metabolites$ACN_Monoisotopic.mass*(1+mzmed.se)
      H.lb = KEGG_metabolites$H_Monoisotopic.mass*(1-mzmed.se)
      H.ub = KEGG_metabolites$H_Monoisotopic.mass*(1+mzmed.se)
      
      peak_long.d = as.data.frame(datasets.all[[d]]@groups) %>%
        dplyr::mutate(Peak=data_raw.d$Peak) %>%
        dplyr::filter(data.table::inrange(mzmed, ACN.lb, ACN.ub) | data.table::inrange(mzmed, H.lb, H.ub)) %>%
        dplyr::mutate(mzlow=mzmed*(1-mzmed.se), mzhigh=mzmed*(1+mzmed.se)) %>%
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
      dplyr::mutate(Peak=make.unique(Peak)) %>%
      dplyr::mutate_at(dplyr::vars(-dplyr::matches("Peak")), list(round)) %>%
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
      zz<<-z
      denominator_peak = dplyr::case_when(z$Dataset[1]=="Lysate"~"278.2/1093", z$Dataset[1]=="Extracellular"~"278.2/1089")
      denominator = z$Intensity[grepl(denominator_peak, z$Peak)]
      z$IntensityLog10=log10(z$Intensity)
      z$Intensity_amitrip = z$Intensity/denominator
      z$IntensityLog10_amitrip=log10(z$Intensity_amitrip)
      z
    })(.)) %>%
    data.frame() %>%
    dplyr::left_join(peak_long.sum, by=c("Dataset", "species.code", "Peak")) %>%
    dplyr::mutate(metabolite_name=dplyr::case_when(State=="lys" & Peak=="298.1/1050" |  State=="extra" & Peak=="298.1/1040" ~ "Duloxetine", T ~ metabolite_name)) %>%
    dplyr::mutate(metabolite_display=dplyr::case_when(State=="lys" & Peak=="298.1/1050" |  State=="extra" & Peak=="298.1/1040" ~ "Duloxetine", T ~ metabolite_display))
  
  
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
      is_significant=ifelse(pvalue<0.05 & abs(FC)>1.1, 1, 0)) %>%
    # is_significant=ifelse(padjust<0.05 & abs(FC)>FCmin, 1, 0)) %>% 
    dplyr::select(-Ufc, -FCmin)
  
  
  #
  # Significally effected peaks
  #
  data_amitrip.effect.sum = data_amitrip.effect %>%
    dplyr::filter(Control!="Dulox") %>% 
    setNames(gsub("^((?!Dataset|Peak|State|species.code|metabolite_name|pathways|metabolite.*).*)$", "\\1.bug", names(.), perl=T)) %>%
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
    dplyr::inner_join(bugs %>% dplyr::select(species.code, species.short), by="species.code") %>%
    # dplyr::mutate(is_significant.bug=pvalue.bug<0.05 & abs(FC.bug)>1.1, is_significant.drug=pvalue.drug<0.05 & abs(FC.bug)>1.1) %>%
    dplyr::filter(Dataset=="Extracellular") %>% # & Treatment %in% c("Cs", "Bu"), data_long.norm
    tidyr::separate_rows(pathways, sep=",") %>%
    tidyr::separate_rows(metabolite_kegg, sep=",") %>%
    
    # Leave only those pathway/metabolite pairs which are present in KEGG database. The reason is that both pathways and metabolites are lists and not all metabolites come from all pathways
    dplyr::inner_join(KEGG_pathways, by=c("metabolite_kegg"="pathway.compound", "pathways"="pathway.id")) %>%
    dplyr::filter(!is.na(pathways) & pathways != "NA") %>%
    
    # Calculate number significant/total metabolites in species
    dplyr::group_by(Dataset, species.short) %>%
    dplyr::mutate(pathway.id=pathways, species.hits=length(unique(Peak[is_significant.bug & is_significant.drug])), species.detected_size=length(unique(Peak))) %>%
    
    # Calculate number significant metabolites in pathways
    dplyr::group_by(Dataset, species.short, pathway.id) %>%
    dplyr::mutate(pathway.detected_size=length(unique(Peak)), pathway.hits=length(unique(Peak[is_significant.bug & is_significant.drug]))) %>%
    
    # Perform hypergeometric/fisher test  
    dplyr::group_by(Dataset, species.short, pathway.id, pathway.name, pathway.hits, pathway.detected_size, species.detected_size, species.hits) %>% 
    dplyr::do((function(z){
      m = matrix(c(
        z$pathway.hits[1],
        z$pathway.detected_size[1]-z$pathway.hits[1],
        z$species.hits[1]-z$pathway.hits[1],
        z$species.detected_size[1]-z$species.hits[1]-z$pathway.detected_size[1] + z$pathway.hits[1]), ncol=2)
      z.test = fisher.test(m)
      data.frame(pvalue=z.test$p.value, odds=unname(z.test$estimate))
    })(.)) %>%
    data.frame() %>%
    
    # Multiple testing adjustment
    dplyr::mutate(
      is_small_pathway=ifelse(pathway.detected_size < 2, "Yes", "No"),
      is_significant=dplyr::case_when(is_small_pathway=="Yes" ~ "No (to few metabolites detected in pathway)", pvalue<0.1~"Yes", T~"No")
    ) %>%
    dplyr::select(species.short, pathway.id, pathway.name, pvalue, odds, is_significant, pathway.hits, pathway.detected_size, species.hits, species.detected_size) %>%
    dplyr::arrange(dplyr::desc(is_significant), pvalue>0.1, odds<1, pvalue)
  
  
  
  #
  # Export pathway enrichment data
  #
  readr::write_tsv(data_amitrip.pathway_hits, "reports/exp3metabolomics_pathway_enrichment.tsv", col_names=T, na="")  
  
  
  
  #
  # Export figure 2c
  #
  data_amitrip.effect.sum_export = data_amitrip.effect.sum %>% 
    dplyr::mutate(Species=dplyr::case_when(
      species.code=="Cs"~"Clostridium saccharolyticum", 
      species.code=="Bu"~"Bacteroides uniformis", 
      T~"No species")) %>%
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
  
  
  png("reports/exp3metabolomics_replicated_correlation.png", width = 2048, height=2048, type="cairo-png")
  data_long.norm_cor = data_long.norm_rep %>% 
    dplyr::group_by(group, Dataset) %>%
    dplyr::summarise(cor=cor(IntensityLog10_amitrip.x, IntensityLog10_amitrip.y, method="spearman", use="pairwise.complete.obs"))
  ggplot(data_long.norm_rep, aes(x=IntensityLog10_amitrip.x, y=IntensityLog10_amitrip.y)) +
    geom_point(color="#00000080", size=0.5) +
    coord_equal() +
    scale_x_continuous(limits =c(-8,2), breaks=seq(-8,2,2)) +
    scale_y_continuous(limits =c(-8,2), breaks=seq(-8,2,2)) +
    labs(x="Rep. 1, Log10 Peak Intensity normalized by Amitriptyline", y="Rep. 2, Log10 Peak Intensity normalized by Amitriptyline") +
    geom_text(aes(x=-4, y=1, label=round(cor, 3)), data=data_long.norm_cor, colour="black", size=10,inherit.aes=FALSE, parse=FALSE) +
    facet_grid(group~Dataset) +
    myTheme +
    theme(strip.text=element_text(size=48), axis.text=element_text(size=32), axis.title=element_text(size=48))
  dev.off()
  
  pdf(file="reports/exp3metabolomics_duloxetine_degradation.pdf", width=11, height=5)
  data_long.norm_dulox = data_long.norm %>% 
    dplyr::filter(metabolite_name=="Duloxetine" & (is.na(species.code) | species.code=="Cs")) %>%
    dplyr::mutate(Subset=dplyr::case_when(
      Treatment=="Cs" ~ "Not exposed to duloxetine",
      Treatment=="CsDulox" ~ "Exposed to duloxetine",
      Treatment=="Dulox" ~ "Duloxetine control",
      T ~ "This should not hapen"
    )) %>%
    dplyr::filter(Dataset=="Extracellular" & Treatment %in% c("CsDulox", "Dulox"))
  ggplot(data_long.norm_dulox) +
    geom_boxplot(aes(x=Subset, y=Intensity_amitrip), fill="#999999") +
    labs(y="Peak intensity normalized by Amitriptyline", x="") +
    myTheme +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18))
  dev.off()
}

