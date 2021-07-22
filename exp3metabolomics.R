library(dplyr)
library(ggplot2)
library(ggridges)
library(xcms)
source("functions.R")

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
  KEGG_pubchem = na.omit(unique(KEGG_metabolites$PubChem.id[KEGG_metabolites$PubChem.id!="n.a"]))
  KEGG_pubchem_smiles = readr::read_csv("data/exp3metabolomics/KEGG_metbolites_smiles.csv")
  sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/CanonicalSMILES,IsomericSMILES,InChIKey,InChI/CSV", paste(KEGG_pubchem[1:100], collapse=","))
  
  
  
  #
  # Load experiments data and peaks
  #
  data_long = data.frame()
  all_peak_long = data.frame()
  peak_long = data.frame()
  datasets.all = list("Extracellular"=extraxset.fillden, "Lysate"=lysxset.fillden, "Click"=clickxset.fillden)
  for(d in names(datasets.all)) {
    # Load raw data
    data_raw.xcms = xcms::groupval(datasets.all[[d]], method="maxint", intensity="into", value="into")
    data_raw.d = data.frame(data_raw.xcms, check.names=F) %>%
      dplyr::mutate(Peak=rownames(data_raw.xcms)) %>%
      data.frame(check.names=F)
    
    all_peak_long.d = as.data.frame(datasets.all[[d]]@groups) %>%
      dplyr::mutate(Peak=data_raw.d$Peak)
    
    # Peaks data
    if(d %in% c("Extracellular", "Lysate")) {
      # This is used to prefilter before interval joining with dplyr because crossing with interval joining is very slow
      mzmed.se = 0.000005
      ACN.lb = KEGG_metabolites$ACN_Monoisotopic.mass*(1-mzmed.se)
      ACN.ub = KEGG_metabolites$ACN_Monoisotopic.mass*(1+mzmed.se)
      H.lb = KEGG_metabolites$H_Monoisotopic.mass*(1-mzmed.se)
      H.ub = KEGG_metabolites$H_Monoisotopic.mass*(1+mzmed.se)
      
      peak_long.d = all_peak_long.d %>%
        dplyr::filter(data.table::inrange(mzmed, ACN.lb, ACN.ub) | data.table::inrange(mzmed, H.lb, H.ub)) %>%
        dplyr::mutate(mzlow=mzmed*(1-mzmed.se), mzhigh=mzmed*(1+mzmed.se)) %>%
        tidyr::crossing(KEGG_metabolites) %>% 
        dplyr::mutate(donor_ACN=ACN_Monoisotopic.mass>=mzlow & ACN_Monoisotopic.mass<=mzhigh, donor_H=H_Monoisotopic.mass>=mzlow & H_Monoisotopic.mass<=mzhigh) %>%
        dplyr::mutate(adduct=dplyr::case_when(donor_H~"[M+H]+", donor_ACN~"[M+ACN+H]+", T~NA_character_)) %>%
        dplyr::filter(donor_ACN | donor_H) %>%
        dplyr::rename(BLK="Blank") %>%
        reshape2::melt(measure.vars=intersect(c("Bu", "BuDulox", "Cs", "CsDulox", "Blank", "Ctrl", "Heat", "Lys"), names(.)), variable.name="Treatment") %>%
        dplyr::filter(gsub("(Cs|Bu).*", "\\1", Treatment)==species.code) %>% # value != 0 (we don't consider peaks annotated by CellZome, as we calculate our own mapping)
        dplyr::select(-value) %>%
        dplyr::mutate(Dataset=d) %>%
        data.frame()
      peak_long = rbind(peak_long, peak_long.d) %>% 
        dplyr::mutate(ChEBI.id=ifelse(ChEBI.id=="n.a", NA_character_, ChEBI.id)) %>%
        dplyr::mutate(PubChem.id=ifelse(PubChem.id=="n.a", NA_character_, ChEBI.id))
    }
    
    
    data.d = data_raw.d %>%
      dplyr::mutate(Peak=make.unique(Peak)) %>%
      dplyr::mutate_at(dplyr::vars(-dplyr::matches("Peak")), list(round)) %>%
      setNames(gsub("^\\d+_", "", colnames(.))) %>%
      tibble::column_to_rownames("Peak") %>%
      data.matrix()
    
    all_peak_long.d = all_peak_long.d[,c("Peak", "mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax")] %>%
      dplyr::mutate(Peak=make.unique(Peak), Dataset=d)
    all_peak_long = rbind(all_peak_long, all_peak_long.d)
    
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
      adduct.display=adduct[1],
      adduct=ifelse(any(!is_phosphorylated), paste(unique(adduct[!is_phosphorylated]), collapse=","), paste(unique(adduct), collapse=",")),
      
      metabolite_formula.display=Sum.formula[1],
      metabolite_formula=ifelse(any(!is_phosphorylated), paste(unique(Sum.formula[!is_phosphorylated]), collapse=","), paste(unique(Sum.formula), collapse=",")),
      metabolite_pubchem.display=PubChem.id[1], 
      metabolite_pubchem=ifelse(any(!is_phosphorylated), paste(unique(PubChem.id[!is_phosphorylated]), collapse=","), paste(unique(PubChem.id), collapse=",")),
      metabolite_kegg.display=KEGG.id[1], 
      metabolite_kegg=ifelse(any(!is_phosphorylated), paste(unique(KEGG.id[!is_phosphorylated]), collapse=","), paste(unique(KEGG.id), collapse=",")),
      metabolite_chebi.display=ChEBI.id[1], 
      metabolite_chebi=ifelse(any(!is_phosphorylated), paste(unique(ChEBI.id[!is_phosphorylated]), collapse=","), paste(unique(ChEBI.id), collapse=",")), 
      metabolite_name.display=Name[1], 
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
    dplyr::mutate(metabolite_name.display=dplyr::case_when(State=="lys" & Peak=="298.1/1050" |  State=="extra" & Peak=="298.1/1040" ~ "Duloxetine", T ~ metabolite_name.display)) %>%
    dplyr::mutate(metabolite_chebi=dplyr::case_when(State=="lys" & Peak=="298.1/1050" |  State=="extra" & Peak=="298.1/1040" ~ "36795", T ~ metabolite_chebi)) %>%
    dplyr::mutate(metabolite_chebi.display=dplyr::case_when(State=="lys" & Peak=="298.1/1050" |  State=="extra" & Peak=="298.1/1040" ~ "36795", T ~ metabolite_chebi.display))
  
  
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
      # zz = data_long.norm %>% dplyr::filter(species.code=="Bu" & Dataset=="Extracellular" & State=="extra")

      zz<<-z
      sp = z$species.code[1]
      dulox_peak = ifelse(z$State[1]=="lys", "298.1/1050", "298.1/1040")
      
      z.drug = data_long.norm %>% dplyr::filter(Treatment=="Dulox" & Dataset==z$Dataset[1] & Peak %in% z$Peak) %>% data.frame()
      z.drug = rbind(z.drug, z %>% dplyr::filter(Treatment==paste0(sp, "Dulox")) %>% data.frame())
      ret.drug = fun_uncertaintest(data_interest=z.drug, pair=c("Dulox", paste0(sp, "Dulox")), dulox_peak=dulox_peak) %>% dplyr::mutate(Control="Dulox") %>% data.frame()
      ret.bug = fun_uncertaintest(data_interest=z, pair=c(sp, paste0(sp, "Dulox")), dulox_peak=dulox_peak) %>% dplyr::mutate(Control=sp) %>% data.frame()
      ret = rbind(ret.bug, ret.drug)
      
      ret
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
    setNames(gsub("^((?!Dataset|Peak|State|species.code|metabolite_name|pathways|metabolite|adduct.*).*)$", "\\1.bug", names(.), perl=T)) %>%
    dplyr::inner_join(data_amitrip.effect %>% dplyr::filter(Control=="Dulox") %>% setNames(gsub("^((?!Dataset|Peak|State|species.code|pathways|metabolite|adduct.*).*)$", "\\1.drug", names(.), perl=T)), by=c("Dataset", "State", "Peak", "pathways", colnames(data_amitrip.effect)[grepl("metabolite|adduct", colnames(data_amitrip.effect))], "species.code")) %>%
    dplyr::mutate(
      FC.both_same_33=dplyr::between(logFC.bug/logFC.drug, 0.75, 1.33),
      Group=dplyr::case_when(grepl("Duloxetine", metabolite_name.display) ~ "Duloxetine", T ~ "Other")
    ) %>%
    dplyr::group_by(species.code, State) %>%
    dplyr::mutate(metabolite_name.display.top10=ifelse(is_significant.bug & is_significant.drug & !is.na(metabolite_name.display) & FC.both_same_33 & FC.bug*FC.drug >= sort((FC.bug*FC.drug)[is_significant.bug & is_significant.drug & !is.na(metabolite_name.display) & FC.both_same_33], decreasing=T)[10], metabolite_name.display, NA)) %>%
    data.frame() %>%
    dplyr::mutate(Group=ifelse(!is.na(metabolite_name.display.top10), "Annotated", Group)) %>% 
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
    dplyr::inner_join(all_peak_long, by=c("Dataset", "Peak"))
  
  data_amitrip.effect.sum_supplementary = data_amitrip.effect.sum_export %>%
    dplyr::select(
      Dataset, Species, Peak,
      mzmed, mzmin, mzmax, rtmed, rtmin, rtmax,
      metabolite_name.display, metabolite_name.all=metabolite_name, metabolite_kegg,
      
      
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
  readr::write_tsv(data_amitrip.effect.sum_supplementary, "reports/exp3metabolomics_hits_scatterplot_data.tsv", col_names=T, na="")
  
  
  #
  # Export MAF table for Metabolights
  #
  metabolights_maf.metabolites = data_amitrip.effect.sum_export %>% 
    dplyr::mutate(metabolite_pubchem.display=as.numeric(metabolite_pubchem.display)) %>%
    dplyr::filter(!is.na(metabolite_name.display)) %>%
    dplyr::arrange(is.na(metabolite_name.display)) %>%
    dplyr::distinct(metabolite_name.display, .keep_all=T)
  
  metabolights_maf.metraw = data_long %>% 
    dplyr::filter(Dataset != "Click") %>%
    dplyr::inner_join(metabolights_maf.metabolites %>% dplyr::distinct(Peak, metabolite_name.display), by="Peak") %>%
    dplyr::group_by(metabolite_name.display, Sample) %>%
    dplyr::summarise(Intensity=max(Intensity, na.rm=T)) %>%
    reshape2::dcast(metabolite_name.display ~ Sample, value.var="Intensity") %>%
    setNames(paste0("M2_", colnames(.)))
  
  metabolights_maf.export = metabolights_maf.metabolites %>%
    dplyr::left_join(KEGG_pubchem_smiles, by=c("metabolite_pubchem.display"="CID")) %>% 
    dplyr::left_join(metabolights_maf.metraw, by=c("metabolite_name.display"="M2_metabolite_name.display")) %>%
    dplyr::mutate(
      M1_database_identifier=paste0("CHEBI:", metabolite_pubchem.display),
      M1_chemical_formula=metabolite_formula.display,
      M1_smiles=CanonicalSMILES,
      M1_inchi=InChI,
      M1_metabolite_identification=metabolite_name.display,
      M1_mass_to_charge=1,
      M1_fragmentation=mzmed,
      M1_modifications=adduct.display,
      M1_charge="",
      M1_retention_time=rtmed,
      M1_taxid="",
      M1_species="",
      M1_database="CHEBI",
      M1_database_version="",
      M1_reliability="",
      M1_uri="",
      M1_search_engine="",
      M1_search_engine_score="",
      M1_smallmolecule_abundance_sub="",
      M1_smallmolecule_abundance_stdev_sub="",
      M1_smallmolecule_abundance_std_error_sub="") %>%
    # dplyr::select(dplyr::starts_with("M2_"))
    dplyr::select(dplyr::starts_with("M1_"), dplyr::starts_with("M2_")) %>%
    setNames(gsub("^(M1|M2)_", "", colnames(.)))
  readr::write_tsv(metabolights_maf.export, "data/exp3metabolomics/metaboligths/m_MTBLS1757_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv", col_names=T, na="")
  
  
  pdf(file="reports/exp3metabolomics_hits_scatterplot.pdf", width=11, height=11)
  ggplot(aes(x=logFC.bug, y=logFC.drug, size=Mean.treated.bug), data=data_amitrip.effect.sum %>% dplyr::filter(Group=="Other" & is_significant.bug & is_significant.drug)) +
    geom_point(aes(color=Group)) +
    geom_point(aes(color=Group), data=data_amitrip.effect.sum %>% dplyr::filter(Group!="Other" & is_significant.bug & is_significant.drug)) +
    ggrepel::geom_text_repel(aes(label=metabolite_name.display.top10), data=data_amitrip.effect.sum %>% dplyr::filter(!is.na(metabolite_name.display.top10) & is_significant.bug & is_significant.drug), size=4) +
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

mapping.samples2raw = function() {
  load("data/exp3metabolomics/clickxset.fillden.noInj1.RData")
  load("data/exp3metabolomics/extraxset.fillden.noInj1.RData")
  load("data/exp3metabolomics/lysxset.fillden.noInj1.RData")
  
  sample.annotations = dplyr::bind_rows(
    clickxset.fillden@phenoData %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(experiment="click", sample_i=1:n()),
    lysxset.fillden@phenoData %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(experiment="lysates", sample_i=1:n()),
    extraxset.fillden@phenoData %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(experiment="extra", sample_i=1:n())
  ) %>% 
    dplyr::mutate(path=paste0("data/exp3metabolomics/raw/", experiment, "/", class, "/", sample, ".mzXML")) %>%
    dplyr::mutate(raw_exists=file.exists(path))
  
  metabolights.samples = sample.annotations %>%
    dplyr::filter(experiment!="click") %>%
    dplyr::mutate(`Sample Name`=sample, `Source Name`="Cultivated bacteria", `Protocol REF`="Sample collection") %>%
    dplyr::mutate(
      `Characteristics[Organism]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"experimental blank",
        grepl("^Cs", class) | class %in% c("Heat", "Lys") ~ "Clostridium saccharolyticum",
        grepl("^Bu", class) | class %in% c("Heat", "Lys") ~ "Bacteroides uniformis"),
      `Term Source REF[Organism]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"MTBLS",
        grepl("^(Cs|Bu)", class) | class %in% c("Heat", "Lys") ~ "NCBITaxon"),
      `Term Accession Number[Organism]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218",
        grepl("^Cs", class) | class %in% c("Heat", "Lys") ~ "http://purl.obolibrary.org/obo/NCBITaxon_610130",
        grepl("^Bu", class) ~ "http://purl.obolibrary.org/obo/NCBITaxon_411479")) %>%
    dplyr::mutate(
      `Characteristics[Organism part]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"experimental control",
        grepl("^Cs|Bu", class) | class %in% c("Heat", "Lys") ~ "whole organism"),
      `Term Source REF[Organism part]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"MTBLS",
        grepl("^(Cs|Bu)", class) | class %in% c("Heat", "Lys") ~ "CARO"),
      `Term Accession Number[Organism part]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218",
        grepl("^Cs|Bu", class) | class %in% c("Heat", "Lys") ~ "http://purl.obolibrary.org/obo/CARO_0000064")) %>%
    dplyr::mutate(
      `Characteristics[Variant]`="normal variability",
      `Term Source REF[Variant]`="PATO",
      `Term Accession Number[Variant]`="http://purl.obolibrary.org/obo/PATO_0045074") %>%
    dplyr::mutate(
      `Characteristics[Sample type]`=dplyr::case_when(
        grepl("^(Ctrl|Dulox)$", class)~"control sample role",
        grepl("^(Bu|Cs)$", class)~"reference reference sample",
        grepl("^(Blank)$", class)~"blank role",
        T~"experiment sample role"),
      `Term Source REF[Sample type]`="AFO",
      `Term Accession Number[Sample type]`=dplyr::case_when(
        grepl("^(Ctrl|Dulox)$", class)~"http://purl.allotrope.org/ontologies/role#AFRL_0000253",
        grepl("^(Bu|Cs)$", class)~"http://purl.allotrope.org/ontologies/role#AFRL_0000258",
        grepl("^(Blank)$", class)~"http://purl.allotrope.org/ontologies/role#AFRL_0000259",
        T~"http://purl.allotrope.org/ontologies/role#AFRL_0000171")) %>%
    dplyr::mutate(
      `Factor[Extraction]`=dplyr::case_when(
        experiment=="lysates" | class %in% c("Lys") ~ "cell lysate",
        experiment=="extra" | class %in% c("Heat") ~ "Extracellular",
        class %in% c("Blank", "Ctrl")~"Extracellular"),
      `Term Source REF[Extraction]`="OBI",
      `Term Accession Number[Extraction]`=dplyr::case_when(
        experiment=="lysates" | class %in% c("Lys") ~ "http://purl.obolibrary.org/obo/OBI_1000036",
        experiment=="extra" | class %in% c("Heat") ~ "http://purl.obolibrary.org/obo/OBI_1000023",
        class %in% c("Blank", "Ctrl")~"http://purl.obolibrary.org/obo/OBI_1000023")
    ) %>%
    dplyr::mutate(
      `Factor[Drug]`=dplyr::case_when(
        class %in% c("Blank")~"no drug",
        class %in% c("Heat","Ctrl","Lys", "Bu", "Cs") | grepl("Dulox", class) ~ "Duloxetine"),
      `Term Source REF[Drug]`=dplyr::case_when(
        class %in% c("Blank")~"MTBLS",
        class %in% c("Heat","Ctrl","Lys", "Bu", "Cs") | grepl("Dulox", class) ~ "CHEBI"),
      `Term Accession Number[Drug]`=dplyr::case_when(
        class %in% c("Blank")~"http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218", 
        class %in% c("Heat","Ctrl","Lys", "Bu", "Cs") | grepl("Dulox", class) ~ "http://purl.obolibrary.org/obo/CHEBI_36796")) %>%
    dplyr::mutate(
      `Factor[Organism]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"no bacteria",
        grepl("^Cs", class)| class %in% c("Heat", "Lys") ~ "Clostridium saccharolyticum",
        grepl("^Bu", class) ~ "Bacteroides uniformis"),
      `Term Source REF[Organism]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"MTBLS",
        grepl("^(Cs|Bu)", class) | class %in% c("Heat", "Lys") ~ "NCBITaxon"),
      `Term Accession Number[Organism]`=dplyr::case_when(
        class %in% c("Blank", "Ctrl", "Dulox")~"http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218",
        grepl("^Cs", class) | class %in% c("Heat", "Lys") ~ "http://purl.obolibrary.org/obo/NCBITaxon_610130",
        grepl("^Bu", class) ~ "http://purl.obolibrary.org/obo/NCBITaxon_411479")) %>%
    dplyr::select(`Sample Name`, `Source Name`, `Protocol REF`, `Source Name`, dplyr::matches("Characteristics|Factor|Term Source REF|Term Accession Number")) %>%
    data.frame(check.names=F) %>%
    setNames(gsub("(Term Source REF|Term Accession Number)\\[.*", "\\1", names(.)))
  readr::write_tsv(metabolights.samples, path="data/exp3metabolomics/metaboligths/s_MTBLS1757.txt",  na="")
  # coupled with Thermo Fisher Q-Exactive plus HRMS
  
  
  metabolights.assays = sample.annotations %>%
    dplyr::filter(experiment!="click") %>%
    dplyr::distinct(sample, experiment) %>%
    dplyr::mutate(ms_name=gsub("^[0-9]+_", "", sample)) %>%
    dplyr::mutate(
      `Sample Name`=sample,
      `Protocol REF[Extraction]`="Extraction",
      `Parameter Value[Post Extraction]`="",
      `Parameter Value[Derivatization]`="",
      `Extract Name`=dplyr::case_when(
        grepl("lys", experiment)~"Lysate",
        grepl("extr", experiment)~"Extracellular"),
      
      `Protocol REF[Chromatography]`="Chromatography",
      `Parameter Value[Chromatography Instrument]`="Vanquish UHPLC system",
      `Term Source REF[Chromatography Instrument]`="OBI",
      `Term Accession Number[Chromatography Instrument]`="http://purl.obolibrary.org/obo/OBI_0000485",
      `Parameter Value[Autosampler model]`="",
      `Parameter Value[Column model]`="ACQUITY UPLC HSS T3 column (2.1 x 100 mm, 1.8 ÃÂµm; Waters)",
      `Parameter Value[Column type]`="reverse phase",
      `Parameter Value[Guard column]`="",
      
      `Labeled Extract Name`="",
      `Label`="",
      `Term Source REF`="",
      `Term Accession Number`="",
      
      `Protocol REF[Mass Spectrometry]`="Mass spectrometry",
      `Parameter Value[Scan polarity]`="positive",
      `Parameter Value[Scan m/z range]`="60 to 800 m/z",
      `Parameter Value[Instrument]`="Q Exactive Plus Hybrid Quadrupol-Orbitrap",
      `Term Source REF[Instrument]`="MS",
      `Term Accession Number[Instrument]`="http://purl.obolibrary.org/obo/MS_1000451",
      `Parameter Value[Ion source]`="electrospray ionization",
      `Term Source REF[Ion source]`="MS",
      `Term Accession Number[Ion source]`="http://purl.obolibrary.org/obo/MS_1000073",
      
      `Parameter Value[Mass analyzer]`="Hybrid quadrupol-orbitrap",
      `Term Source REF[Mass analyzer]`="MS",
      `Term Accession Number[Mass analyzer]`="http://purl.obolibrary.org/obo/MS_1000451",
      `MS Assay Name`=ms_name,
      `Raw Spectral Data File`=paste0(`Sample Name`, ".mzXML"),      
      
      `Protocol REF[Data transformation]`="Data transformation",
      `Normalization Name`="Amitriptyline reference",
      `Derived Spectral Data File`="",
      
      `Protocol REF[Metabolite identification]`="Metabolite identification",
      `Data Transformation Name`="",
      `Metabolite Assignment File`="m_MTBLS1757_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
    ) %>%
    dplyr::select(-sample, -experiment) %>%
    data.frame(check.names=F) %>%
    setNames(gsub("(Term Source REF|Term Accession Number|Protocol REF)\\[.*", "\\1", names(.)))
  
  # readr::write_tsv(metabolights.assays, path="C:\\Users\\sandrejev\\Desktop\\a_MTBLS1757_LC-MS_positive_reverse-phase_metabolite_profiling.txt",  na="")
  readr::write_tsv(metabolights.assays, path="data/exp3metabolomics/metaboligths/a_MTBLS1757_LC-MS_positive_reverse-phase_metabolite_profiling.txt",  na="")
}

export.raw_msms = function() 
{
  load("data/exp3metabolomics/lysxset.fillden.noInj1.RData")
  load("data/exp3metabolomics/extraxset.fillden.noInj1.RData")
  load("data/exp3metabolomics/xfrag.lys.rda")
  load("data/exp3metabolomics/xfrag.extra.rda")
  # extraxset.fillden@filepaths = paste0("data/exp3metabolomics/raw_flat/", basename(extraxset.fillden@filepaths))
  # xfrag.extra <- xcmsFragments(extraxset.fillden)
  # lysxset.fillden@filepaths = paste0("data/exp3metabolomics/raw_flat/", basename(lysxset.fillden@filepaths))
  # xfrag.lys <- xcmsFragments(lysxset.fillden)
  # save(xfrag.extra, file="data/exp3metabolomics/xfrag.extra.rda")
  # save(xfrag.lys, file="data/exp3metabolomics/xfrag.lys.rda")
  
  
  sample.annotations = dplyr::bind_rows(
    lysxset.fillden@phenoData %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(experiment="lysates", sample_i=1:n(), Dataset="Lysate"),
    extraxset.fillden@phenoData %>% tibble::rownames_to_column("sample") %>% dplyr::mutate(experiment="extra", sample_i=1:n(), Dataset="Extracellular")
  ) %>% 
    dplyr::mutate(path=paste0("data/exp3metabolomics/raw/", experiment, "/", class, "/", sample, ".mzXML")) %>%
    dplyr::mutate(raw_exists=file.exists(path)) %>%
    dplyr::mutate(species_short=dplyr::case_when(
      grepl("^Cs", class)~"C. saccharolyticum",
      grepl("^Bu", class)~"B. uniformis",
      grepl("Blank", class)~"No bug")
    ) %>%
    dplyr::mutate(drug_long=dplyr::case_when(
      grepl("Dulox$", class)~"Duloxetine",
      T~"No drug")
    )
  
  ms1.all = list("Extracellular"=extraxset.fillden, "Lysate"=lysxset.fillden)
  ms2.all = list("Extracellular"=xfrag.extra, "Lysate"=xfrag.lys)
  hits_peaks = readr::read_tsv("reports/exp3metabolomics_hits_scatterplot_data.tsv") %>%
    # dplyr::filter(is_significant.treated_vs_dulox==1 & is_significant.treated_vs_species==1) %>%
    dplyr::filter(!is.na(metabolite_name.display) & metabolite_name.display %in% c("Guanosine", "Inosine", "Adenosine", "Aminoimidazole ribotide", "Xanthine")) %>%
    dplyr::group_by(Dataset) %>%
    dplyr::do((function(z){
      zz<<-z
      xs = ms1.all[[z$Dataset[1]]]
      xs.ms2 = data.frame(ms2.all[[z$Dataset[1]]]@peaks) %>%
        dplyr::mutate(class=xs@phenoData$class[Sample], msLevel, sample_i=Sample, sample_name=rownames(xs@phenoData)[Sample], path=xs@filepaths[Sample])
      groups = xcms::groupval(xs, method="maxint", intensity="into", value="into")
      groupidx = xs@groupidx[rownames(groups) %in% z$Peak]
      groups = groups[rownames(groups) %in% z$Peak,]
      peak2ms1 = do.call(rbind, lapply(1:length(groupidx), function(y) {
        yy<<-y
        data.frame(Dataset=z$Dataset[1], Peak=rownames(groups)[y], Peak.MS1=groupidx[[y]], sample.ms1=rownames(xs@phenoData)[xs@peaks[groupidx[[y]],"sample"]], class.ms1=xs@phenoData[xs@peaks[groupidx[[y]],"sample"],"class"], intensity.ms1=xs@peaks[groupidx[[y]],"into"])
      })) %>% dplyr::filter(Peak %in% z$Peak)
      peak2ms2 = peak2ms1 %>%
        dplyr::inner_join(z %>% dplyr::select(Dataset, Peak, mzmed, mzmin, mzmax, rtmed, rtmin, rtmax, metabolite_name.display, Mean.treated_species.ms1=Mean.treated_species, is_significant.treated_vs_dulox, is_significant.treated_vs_species), by=c("Dataset", "Peak")) %>%
        dplyr::inner_join(xs.ms2 %>% dplyr::filter(msLevel==2) %>% dplyr::select(class, path, sample_i, sample_name, Peak.MS2=peakID, MSnParentPeakID, rt.ms2=rt, mz.ms2=mz, intensity.ms2=intensity), by=c("Peak.MS1"="MSnParentPeakID"))
      peak2ms2
    })(.)) %>%
    dplyr::ungroup()
  
  hits_peaks.ggplot = hits_peaks %>%
    # dplyr::filter(grepl("inj2", sample_name)) %>%
    dplyr::mutate(injection=gsub(".*inj([0-9]).*", "inj\\1", sample_name)) %>%
    dplyr::mutate(replicate=gsub(".*rep([0-9]).*", "rep\\1", sample_name)) %>%
    dplyr::mutate(metabolite_name.display=factor(as.character(metabolite_name.display), unique(as.character(metabolite_name.display)))) %>%
    dplyr::group_by(Dataset, class, sample_name, Peak) %>%
    dplyr::arrange(dplyr::desc(intensity.ms2)) %>%
    dplyr::mutate(fragment=as.character(1:n())) %>%
    dplyr::slice(1:10) %>%
    dplyr::filter(intensity.ms2>10e4)
  breaks.mz = seq(50, round(max(hits_peaks.ggplot$mz.ms2), 1), by=0.1)
  hits_peaks.ggplot = hits_peaks.ggplot %>%
    dplyr::mutate(mz_cut.ms2=cut(mz.ms2, breaks.mz, labels=F), mz_cut.ms2=breaks.mz[mz_cut.ms2]+0.5)
  ggplot(hits_peaks.ggplot) +
    # geom_bar(aes(x=mz_cut.ms2, y=intensity.ms2, fill=Dataset, group=paste(Dataset, fragment)), stat="identity", position="dodge", alpha=0.5) +
    geom_point(aes(x=mz_cut.ms2, y=rt.ms2, size=intensity.ms2, color=Dataset, group=paste(Dataset, fragment)), stat="identity", position="dodge", alpha=0.5) +
    ggrepel::geom_text_repel(aes(x=mz_cut.ms2, y=rt.ms2, color=Dataset, label=paste0(replicate, injection, ":", fragment), group=paste(Dataset, fragment))) +
    facet_grid(class~metabolite_name.display, scale="free_y")
  # with(peak2ms2 %>% dplyr::filter(intensity.mz2>10e6), table(class, metabolite_name.display))
  
  ggplot(hits_peaks.ggplot %>% dplyr::distinct(class, metabolite_name.display, Peak, .keep_all=T)) +
    geom_rect(aes(xmin=mzmin, xmax=mzmax, ymin=rtmin, ymax=rtmax, color=metabolite_name.display), alpha=0.5) + 
    facet_grid(~class, scale="free_y")
  
  sample.peaks = hits_peaks %>% 
    dplyr::arrange(dplyr::desc(is_significant.treated_vs_dulox), dplyr::desc(is_significant.treated_vs_species), dplyr::desc(intensity.ms1)) %>% 
    dplyr::group_by(Dataset, metabolite_name.display, Peak, Peak.MS1) %>%
    dplyr::slice(1)
  
  # d = sample.peaks %>%
  #   dplyr::filter(sample=="150928_CsDulox_lys_rep1_inj3" & metabolite_name.display=="Aminoimidazole ribotide") %>%
  #   dplyr::filter(is_significant.treated_vs_dulox==1 & is_significant.treated_vs_species) %>%
  #   dplyr::filter(intensity.ms1==max(intensity.ms1)) %>%
  #   dplyr::arrange(dplyr::desc(intensity.ms2))
  
  sample.annotations_examples = sample.annotations %>%
    dplyr::filter(experiment!="click") #%>%
  # dplyr::group_by(class) %>%
  # dplyr::filter(class=="Blank" & 1:n()==1 | grepl("inj2", sample) & grepl("rep3", sample)) %>%
  # dplyr::ungroup()
  sample.annotations_examples = sample.annotations %>%
    dplyr::filter(experiment!="click")
  samples.raw = lapply(sample.annotations_examples$sample, function(x){
    path = sample.annotations_examples$path[sample.annotations_examples$sample==x]
    MSnbase::readMSData(path, mode = "onDisk")
  })
  names(samples.raw) = sample.annotations_examples$sample
  
  # loop through each sample
  samples.ms2_report = lapply(1:nrow(sample.annotations_examples), function(x){
    xx<<-x
    sample = sample.annotations_examples$sample[x]
    sample_i = sample.annotations_examples$sample_i[x]
    dataset = sample.annotations_examples$Dataset[x]
    raw = samples.raw[[sample]]
    sample.dataset_peaks = sample.peaks %>% dplyr::filter(Dataset==dataset)
    
    # sample.dataset_peaks = data.frame(ms2.all[[dataset]]@peaks) %>%
    #   dplyr::filter(peakID %in% sample.dataset_peaks.x$Peak.MS1 & Sample==sample_i)
    
    # loop through each peak
    raw_ms2.samples.x = lapply(1:nrow(sample.dataset_peaks), function(y){
      yy<<-y
      raw.f.ms1 = MSnbase::filterRt(raw, c(sample.dataset_peaks$rtmin[y], sample.dataset_peaks$rtmax[y]), msLevel.=1) 
      raw.f.ms1 = MSnbase::filterMz(raw.f.ms1, c(sample.dataset_peaks$mzmin[y], sample.dataset_peaks$mzmax[y]), msLevel.=1)
      raw.f.ms1 = MSnbase::filterMsLevel(raw.f.ms1, 1)
      raw.f.ms1_scan = MSnbase::scanIndex(raw.f.ms1)
      # raw.f.ms1_feature = fData(raw.f.ms1)[sapply(raw.f.ms1_mz, length)>0 & sapply(raw.f.ms1_rt, length)>0 & sapply(raw.f.ms1_scan, length)>0,]
      
      # intensity.ms1 = MSnbase::intensity(raw.f.ms1)
      # mz.ms1 = MSnbase::mz(raw.f.ms1)
      # peaks.ms1 = do.call(rbind, lapply(names(mz.ms1), function(g) {
      #   gg<<-g
      #   if(length(mz.ms1[[g]])==0 | length(intensity.ms1[[g]])==0) return(data.frame())
      #   data.frame(Feature=g, mz=mz.ms1[[g]], intensity=intensity.ms1[[g]], stringsAsFactors=F)
      # })) %>% data.frame()
      
      if(length(raw.f.ms1_scan)==0) { return(data.frame()) }
      ppm = 1e6*max(sample.dataset_peaks$mzmed[y]-sample.dataset_peaks$mzmin[y], sample.dataset_peaks$mzmax[y]-sample.dataset_peaks$mzmed[y])/sample.dataset_peaks$mzmed[y]
      ppm = 50
      raw.f.ms2 = MSnbase::filterPrecursorScan(raw, raw.f.ms1_scan)
      raw.f.ms2 = MSnbase::filterPrecursorMz(raw.f.ms2, sample.dataset_peaks$mzmed[y], ppm)
      raw.f.ms2 = MSnbase::filterMsLevel(raw.f.ms2, 2)
      
      intensity.ms2 = MSnbase::intensity(raw.f.ms2)
      mz.ms2 = MSnbase::mz(raw.f.ms2)
      if(is.null(mz.ms2)) { return(data.frame()) }
      
      # Combine all MS/MS into one dataframe for each MS1+sample
      output.intensity.g = do.call(rbind, lapply(names(mz.ms2), function(g) {
        gg<<-g
        data.frame(Feature=g, mz=mz.ms2[[g]], intensity=intensity.ms2[[g]], stringsAsFactors=F)
      })) %>% data.frame()
      output.intensity = cbind(data.frame(sample.dataset_peaks[y,]), output.intensity.g)
      output.intensity
    })
    
    do.call(rbind, raw_ms2.samples.x)
  })
  samples.ms2_report = do.call(rbind, samples.ms2_report) %>%
    dplyr::mutate(species_short=dplyr::case_when(
      grepl("^Cs", class)~"C. saccharolyticum",
      grepl("^Bu", class)~"B. uniformis",
      grepl("Blank", class)~"No bug")
    ) %>%
    dplyr::mutate(drug_long=dplyr::case_when(
      grepl("Dulox$", class)~"Duloxetine",
      T~"No drug")
    ) %>%
    dplyr::group_by(metabolite_name.display, Dataset, class) %>%
    dplyr::arrange(dplyr::desc(intensity)) %>%
    dplyr::filter(sample_name==sample_name[1] & Peak==Peak[1] & Feature==Feature[1])  %>%
    dplyr::ungroup()
  
  readr::write_tsv(samples.ms2_report, path="reports/exp3metabolomics_ms2_plot.tsv", na="")
  samples.ms2_report = readr::read_tsv("reports/exp3metabolomics_ms2_plot.tsv")
  ms2peaks = samples.ms2_report %>%
    dplyr::group_by(sample_name, Peak, Feature) %>%
    dplyr::mutate(intensity.raw=intensity, intensity=intensity/max(intensity)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(intensity>0.01) %>%
    dplyr::group_by(metabolite_name.display, Dataset, class) %>%
    dplyr::arrange(dplyr::desc(intensity.raw)) %>%
    dplyr::filter(sample_name==sample_name[1] & Peak==Peak[1] & Feature==Feature[1])  %>%
    dplyr::ungroup() %>%
    dplyr::select(Dataset, species_short, drug_long, sample_name, class, metabolite_name.display, Peak, Feature, mzmed.ms1=mzmed, mzmin.ms1=mzmin, mzmax.ms1=mzmax, rtmed.ms1=rtmed, rtmin.ms1=rtmin, rtmax.ms1=rtmax, rt.ms2, mz.ms2=mz, intensity_rel.ms2=intensity)
  
  readr::write_tsv(ms2peaks, path="reports/exp3metabolomics_ms2peaks.tsv", na="")
  ms2peaks = readr::read_tsv("reports/exp3metabolomics_ms2peaks.tsv")
  ms2ref = readr::read_tsv("data/exp3metabolomics/ms2ref.tsv") %>%
    tidyr::crossing(ms2peaks %>% dplyr::distinct(Dataset, species_short)) %>%
    dplyr::mutate(sample_name=paste(source, collisionEnergy), class=sample_name, Peak="reference", Feature="reference")
  ms2ref_our = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp3metabolomics/exp3metabolomics_ms2_v2.tsv") %>%
    reshape2::melt(id.var="mz.ms2", value.name="intensity.ms2", variable.name="sample_name") %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(intensity_rel.ms2=intensity.ms2/max(intensity.ms2, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(intensity_rel.ms2>0.05) %>%
    dplyr::mutate(Feature=paste0("F", 1:n()), Peak=paste0("P", 1:n())) %>%
    dplyr::mutate(collisionEnergy=35, source="internal", class=paste0(source, " (", gsub("-.*", "", sample_name), ")")) %>% 
    dplyr::mutate(metabolite_name.display=stringr::str_to_title(gsub(".*-", "", sample_name))) %>%
    tidyr::crossing(expand.grid(Dataset=c("Extracellular", "Lysate"), species_short=c("C. saccharolyticum", "B. uniformis")))
  
  ms2peaks.ggplot = dplyr::bind_rows(ms2peaks, ms2ref, ms2ref_our) %>%
    dplyr::group_by(species_short, metabolite_name.display, source, collisionEnergy) %>%
    dplyr::do((function(z) {
      dplyr::bind_rows(
        z,
        z %>% dplyr::mutate(mz.ms2=mz.ms2-0.0000001, intensity_rel.ms2=0),
        z %>% dplyr::mutate(mz.ms2=mz.ms2-0.0000002, intensity_rel.ms2=NA_real_),
        z %>% dplyr::mutate(mz.ms2=mz.ms2+0.0000001, intensity_rel.ms2=0),
        z %>% dplyr::mutate(mz.ms2=mz.ms2+0.0000002, intensity_rel.ms2=NA_real_)
      ) %>% dplyr::arrange(mz.ms2)
    })(.)) %>%
    dplyr::group_by(species_short, metabolite_name.display, Dataset) %>%
    dplyr::filter(any(Peak!="reference")) %>%
    dplyr::ungroup() 
  
  pdf(file="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/reports/exp3metabolomics_msms_noMZppm50_3.pdf", width=12, height=5)
  for(s in c("C. saccharolyticum", "B. uniformis")) {
    mz.range = c(70, 160)
    ms2peaks.ggplot.s = ms2peaks.ggplot %>%
      dplyr::filter(species_short==s | class=="Blank") %>%
      dplyr::group_by(metabolite_name.display) %>%
      dplyr::filter(any(grepl("Bu|Cs|Blank", class))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(dplyr::between(mz.ms2, mz.range[1], mz.range[2])) 
    
    ms2peaks.ggplot.s %>% dplyr::group_by(class, Dataset, metabolite_name.display) %>%
      dplyr::summarise(mz.ms2=mz.ms2[which.max(intensity_rel.ms2)], intensity_rel.ms2=max(intensity_rel.ms2, na.rm=T)) %>%
      dplyr::arrange(metabolite_name.display, class, Dataset) %>%
      as.data.frame()

        
    p = ggplot(ms2peaks.ggplot.s) +
      ggridges::geom_ridgeline(aes(x=mz.ms2, y=class, height=intensity_rel.ms2, group=paste(sample_name, Peak, Feature), color=class), scale=0.8) +
      # geom_hline(aes(yintercept=class, color=class), data=ms2peaks.ggplot.s %>% dplyr::distinct(metabolite_name.display, Dataset, class)) +
      facet_grid(metabolite_name.display~Dataset) +
      labs(y="", x="m/z", title=s) +
      scale_x_continuous(minor_breaks=seq(mz.range[1], mz.range[2], 1), breaks=seq(mz.range[1], mz.range[2], 20), limits=mz.range) +
      scale_color_manual(values=c("embl 35"="#333333", "metlin 20"="#333333", "metlin 40"="#333333", "Bu"="#E41A1C", "Cs"="#E41A1C", "BuDulox"="#377EB8", "CsDulox"="#377EB8", "Blank"="#4DAF4A", "internal (std)"="#333333", "internal (dulo)"="#333333")) +
      theme_bw()
    
    print(p) #"#984EA3"
  }
  dev.off()
}

library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)

new_ms2 = function() {
  ms2ref_our = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp3metabolomics/exp3metabolomics_ms2_v2.tsv") %>%
    reshape2::melt(id.var="mz.ms2", value.name="intensity.ms2", variable.name="sample_name") %>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(intensity_rel.ms2=intensity.ms2/max(intensity.ms2, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(intensity_rel.ms2>0.05) %>%
    # dplyr::mutate(Peak=paste0("P", 1:n())) %>%
    dplyr::mutate(class=gsub("-.*", "", sample_name)) %>% 
    dplyr::mutate(metabolite_name.display=stringr::str_to_title(gsub(".*-", "", sample_name))) %>%
    dplyr::mutate(intensity_rel.ms2=ifelse(class=="std", -intensity_rel.ms2, intensity_rel.ms2))
  
  
  pdf(file="msms.pdf", width=12, height=5)
  mz.range = c(70, 160)

  ggplot(ms2ref_our) +
    # geom_line(aes(x=mz.ms2, y=intensity_rel.ms2, group=sample_name, color=class), scale=0.8) +
    geom_segment(aes(x=mz.ms2, xend=mz.ms2, y=0, yend=intensity_rel.ms2, group=sample_name, color=class), scale=0.8) +
    geom_hline(yintercept=0) +
    labs(y="", x="m/z") +
    scale_x_continuous(minor_breaks=seq(mz.range[1], mz.range[2], 1), breaks=seq(mz.range[1], mz.range[2], 20), limits=mz.range) +
    scale_color_manual(values=c("std"="#FF0000", "dulo"="#0000FF")) +
    theme_bw()
   
  dev.off()
}
