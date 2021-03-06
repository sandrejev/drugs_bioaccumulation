dir.create("reports", showWarnings=F)
library(rCharts)
library(networkD3)
library(dplyr)
library(readr)
library(cowplot)
library(xlsx)
library(ggplot2)
library(grid)
source("functions.R")

#
# Differential analysis for degradation dataset
#
ttest.degrad.batch = function(data.degrad){
  meanmax = data.degrad %>%
    dplyr::group_by(Plate, GroupName) %>%
    dplyr::summarise(Mean.Max=mean(maxOD))
  
  ttest.diffSD = data.degrad %>%
    dplyr::group_by(GroupName,Plate,Batch,Species.x) %>%
    dplyr::do(t.test2(mean(.$DrugRatio), .$BatchMean[1], sd(.$DrugRatio),.$SD.ctrls[1], length(.$DrugRatio),.$N.ctrls[1], equal.variance=T)) %>%
    data.frame() %>%
    dplyr::mutate(p_adjust=p.adjust(p_value,method="BH"))
  
  ttest.sameSD = data.degrad %>%
    dplyr::group_by(GroupName, Plate, Batch, Species.x) %>%
    dplyr::do(t.test2(mean(.$DrugRatio),.$BatchMean[1], .$SD.ctrls[1], .$SD.ctrls[1], length(.$DrugRatio), .$N.ctrls[1], equal.variance=T)) %>%
    data.frame() %>%
    dplyr::mutate(p_adjust=p.adjust(p_value,method="BH")) %>%
    dplyr::inner_join(unique(data.degrad %>% dplyr::select(GroupName, Plate, BatchMean, growth)), by=c("GroupName", "Plate")) %>%
    dplyr::mutate(DiffOfMeans = Difference.of.means/BatchMean) %>%
    dplyr::inner_join(meanmax, by=c("GroupName", "Plate")) %>%
    dplyr::mutate(Norm.Max=DiffOfMeans/Mean.Max)
  
  data.get.new = list(ttest.sameSD, ttest.diffSD)
}

preprocess.hits_uplc = function()
{
  padjst.th = 0.05
  MedianDiff.general_th = -0.1
  MedianDiff.bioaccumulation_th = 0.1
  MedianDiffFold.bioaccumulation_th = 0.7
  
  drugs = readr::read_delim("data/drug_map.tsv", "\t") %>% dplyr::filter(!grepl("high concentration", drug.comment))
  bugs = readr::read_delim("data/bug_map.tsv", "\t")
  data = readr::read_delim("data/exp2metabolomics/data.depletionmodeassay_long.csv", ",") 
  data = data %>%
    dplyr::rename("species.long"="Bugslong") %>%
    dplyr::inner_join(drugs, by=c("Drugs"="drug.short2")) %>% 
    dplyr::left_join(bugs %>% dplyr::select(species.short, species.long), by=c("species.long")) %>%
    dplyr::filter(Ctrl != "zero") %>% # I don't know what is 'med',  & is.na(drug.uplc_excluded) #  & Bugs != "med"
    data.frame()
  
  data.sum = data %>%
    dplyr::group_by(drug.long, drug.short, species.long, species.short, Extraction, drug.known_activity, drug.uplc_excluded) %>% 
    dplyr::summarise(MedianDiff=median(DiffCtrlSample[Ctrl=="smpl"]), pvalue=wilcox.test(DiffCtrlSample[Ctrl=="ctrl"], DiffCtrlSample[Ctrl=="smpl"])$p.value, n.ctrl=sum(Ctrl=="ctrl"), n.sample=sum(Ctrl=="smpl")) %>%
    data.frame() %>%
    dplyr::mutate(padjst=p.adjust(pvalue, method="BH"), s=Extraction=="Supernatant") %>%
    dplyr::group_by(drug.long, drug.short, species.short, species.long, drug.known_activity) %>% 
    dplyr::summarise(
      n.super_ctrl=n.ctrl[s], n.super_sample=n.sample[s], n.total_ctrl=n.ctrl[!s], n.total_sample=n.sample[!s],
      padjst.super=padjst[s], padjst.total=padjst[!s], MedianDiff.super=MedianDiff[s],  MedianDiff.total=MedianDiff[!s], 
      MedianDiffFold=MedianDiff.total/MedianDiff.super, 
      MedianDiffDiff=pmin(0, MedianDiff.total)-pmin(0, MedianDiff.super),
      interaction = dplyr::case_when(
        any(!is.na(drug.uplc_excluded)) ~ "Excluded",
        all(MedianDiff>MedianDiff.general_th) | all(padjst>padjst.th) ~ "No activity",
        
        # This
        padjst.super<padjst.th & MedianDiffDiff>=MedianDiff.bioaccumulation_th & MedianDiffFold<MedianDiffFold.bioaccumulation_th~ "Bioaccumulation",
        padjst.super<padjst.th & MedianDiffDiff>=MedianDiff.bioaccumulation_th & MedianDiff.total>-MedianDiff.bioaccumulation_th~ "Bioaccumulation",
        
        padjst.super<padjst.th & any(MedianDiff>MedianDiff.general_th) & MedianDiffDiff>=MedianDiff.bioaccumulation_th & MedianDiffFold<MedianDiffFold.bioaccumulation_th~ "Bioaccumulation",
        padjst.super<padjst.th & any(MedianDiff<MedianDiff.general_th) & MedianDiffDiff>=MedianDiff.bioaccumulation_th & MedianDiff.total>MedianDiff.general_th~ "Bioaccumulation",
        
        
        any(padjst<padjst.th) ~ "Biotransformation",
        T ~ "No activity"
      ))
  
  #
  # Write to file
  #
  readr::write_delim(data.sum, "data/exp2metabolomics/hits.tsv", "\t")
}

number_of_replicates = function()
{
  drug_map = readr::read_delim("data/drug_map.tsv", "\t") %>% dplyr::filter(!grepl("high concentration", drug.comment))
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")
  bugs = readr::read_delim("data/bug_map.tsv", "\t")
  load("data/exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
  
  #
  # Statistics for growth data
  #
  data.growth = read.table("data/exp0growth/2016-11-28_curves_annotation.tab", stringsAsFactors=F, na.strings="", sep="\t", quote="", header=T)
  data.growth_reps = data.growth %>%
    dplyr::mutate(drug.long=gsub("\\|.*", "", ConditionSpecies), solution=gsub(".*\\|", "", ConditionSpecies)) %>%
    dplyr::group_by(ConditionSpecies, Species, File, drug.long, solution, TechnicalReplicates) %>%
    dplyr::do((function(z) {
      z.wells = as.numeric(unlist(strsplit(gsub("[^0-9]+", ",", z$TechnicalReplicates), ",")))
      data.frame(TechnicalReplicatesCount=z.wells[3] - z.wells[2] + 1)
    })(.)) %>% 
    data.frame() %>%
    dplyr::group_by(ConditionSpecies, Species, drug.long, solution) %>%
    dplyr::summarise(n_bio=length(unique(File)), n_tech=sum(TechnicalReplicatesCount)) 
  data.growth_control_reps = data.growth_reps %>% dplyr::filter(solution==drug.long) %>% data.frame()
  data.growth_sample_reps = data.growth_reps %>% dplyr::filter(solution!=drug.long) %>% data.frame()
  data.growth_final_reps = data.growth_sample_reps %>% 
    dplyr::select(-drug.long) %>%
    dplyr::left_join(data.growth_control_reps %>% dplyr::select(-drug.long, -ConditionSpecies), by=c("Species", "solution")) %>%
    setNames(gsub("(\\.x)$", ".sample", names(.), perl=T)) %>%
    setNames(gsub("(\\.y)$", ".control", names(.), perl=T)) %>%
    dplyr::mutate(Species=dplyr::case_when(
      Species=="F. nucleatum subsp. nucleatum"~"F. nucleatum", 
      Species=="Mix Depletion"~"Mix Degrad",
      Species=="Mix NoDepletion"~"Mix No",
      T~Species))
  
  #
  # Statistics for statistical report form
  #
  plates_with_growth = unique(data.clean %>% dplyr::filter(Status=="GMM" & growth=="growth") %>% .$Plate)
  data.clean_techreps = data.clean %>% 
    dplyr::filter(Status!="sample" | Plate %in% plates_with_growth) %>%
    dplyr::inner_join(drug_map %>% dplyr::select(drug.short, drug.uplc_excluded), by=c("GroupName"="drug.short"))  %>%
    dplyr::group_by(Species.x, GroupName, Status, Batch, N.ctrls, drug.uplc_excluded, Replicate) %>%
    dplyr::mutate(TechnicalReplicates=length(unique(Well)), N.ctrls_count=length(unique(N.ctrls)))
  data.clean_sample_techreps = data.clean_techreps %>% 
    dplyr::filter(Status=="sample") 
  
  data.clean_sample_bioreps = data.clean_sample_techreps %>%
    dplyr::group_by(Species.x, GroupName, Status, Batch, drug.uplc_excluded) %>%
    dplyr::summarise(BiologicalReplicates=length(unique(Replicate)), Plates=paste(unique(Plate), collapse=","), N.ctrls_count=length(unique(N.ctrls)))
  
  
  #
  # UPLC experiment
  #
  data.uplc = readr::read_delim("data/exp2metabolomics/data.depletionmodeassay_long.csv", ",") 
  data.uplc_techrep = data.uplc %>% 
    dplyr::filter(Ctrl != "zero") %>% 
    dplyr::mutate(species.long=tidyr::replace_na(Bugslong, ""), Drugs=tidyr::replace_na(as.character(Drugs), "")) %>%
    dplyr::inner_join(drug_map, by=c("Drugs"="drug.short2")) %>% 
    dplyr::inner_join(bug_map %>% dplyr::select(species.short, species.long), by=c("species.long")) %>%
    dplyr::group_by(Extraction, Ctrl, SampleName, Sample.Set.Name,drug.long, species.long, drug.uplc_excluded, drug.known_activity, Replicate) %>%
    dplyr::summarise(TechnicalReplicates=length(Replicate)) %>%
    data.frame()
  
  data.uplc_biorep = data.uplc_techrep %>%
    group_by(Extraction, Ctrl, species.long, drug.long, drug.uplc_excluded, drug.known_activity) %>%
    dplyr::summarise(BiologicalReplicates=length(unique(Replicate)), TechnicalReplicates=sum(TechnicalReplicates))
  
  data.uplc_biorep = data.uplc_techrep %>%
    group_by(Extraction, Ctrl, species.long, drug.long, drug.uplc_excluded, drug.known_activity) %>%
    dplyr::summarise(BiologicalReplicates=length(unique(Replicate)), TechnicalReplicates=sum(TechnicalReplicates))
  
}

exp012depletion.sankey = function()
{
  drug_map = readr::read_delim("data/drug_map.tsv", "\t")
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")
  
  #
  # Load depletion experiment (#1) data for LEFT side of the sankey plot
  #
  load("data/exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
  plates_with_growth = unique(data.clean %>% dplyr::filter(Status=="GMM" & growth=="growth") %>% .$Plate)
  data.degrad = data.clean %>% dplyr::filter(Status=="sample" & growth=="growth" & dummy==1 & Plate %in% plates_with_growth)
  data.get = ttest.degrad.batch(data.degrad) # initial screen
  
  #
  # Calculate edges for LEFT side of the sankey plot and join with UPLC experiment (#2)
  #
  hits_uplc = readr::read_delim("data/exp2metabolomics/hits.tsv", "\t")
  hits_degrad.all = data.get[[1]] %>%
    dplyr::mutate(Species.x=tidyr::replace_na(Species.x, ""), GroupName=tidyr::replace_na(as.character(GroupName), "")) %>%
    dplyr::left_join(data.clean %>% dplyr::select(GroupName, Plate, DiffToOwn) %>% unique(), by=c("GroupName", "Plate")) %>%
    dplyr::group_by(GroupName, Species.x) %>%
    dplyr::summarize(
      p_value.screen=min(p_value),
      p_adjust.screen=min(p_adjust),
      replicates=length(GroupName), 
      nhits=sum(-0.3>DiffToOwn), 
      nhits.plates=paste(Plate[-0.3>DiffToOwn], collapse=","), 
      MeanDiffToOwn=mean(DiffToOwn[-0.3>DiffToOwn]), 
      side="left") %>%
    data.frame()
  
  hits_degrad = hits_degrad.all %>%
    dplyr::mutate(target=as.character(GroupName), source=as.character(Species.x), value=abs(MeanDiffToOwn)) %>%
    dplyr::left_join(hits_uplc %>% dplyr::select(species.long, drug.short, interaction), by=c("source"="species.long", "target"="drug.short")) %>%
    dplyr::inner_join(drug_map %>% dplyr::select(drug.short, drug.known_activity, drug.uplc_excluded), by=c("target"="drug.short")) %>%
    dplyr::mutate(drug.known_activity=dplyr::case_when(
      GroupName=="leva" & grepl("infantis|bolteae|ramosum|comes", Species.x) ~ "Depleted",
      GroupName=="digo" & !grepl("lenta", Species.x)                         ~ NA_character_,
      T ~ drug.known_activity)) %>%
    dplyr::mutate(interaction=dplyr::case_when(
      !is.na(drug.known_activity) ~ "Previously_known",
      !is.na(drug.uplc_excluded)  ~ "Excluded",
      is.na(interaction)          ~ "Untested",
      TRUE                        ~ interaction)) %>%
    dplyr::filter(!is.na(drug.known_activity) & nhits>0 | nhits/replicates>=0.5 | nhits>1)
  
  
  
  #
  # Load growth effect experiment (#0) for RIGHT side of the sankey plot
  #
  hits_growth.all = read.csv("data/exp0growth/curves.rel_annotation_2016-11-28.tab",sep="\t", header=T) %>%
    dplyr::filter(!grepl("pyri", cond.org)) %>% # Pyri is not appearing anywhere in the drug list
    dplyr::mutate(side="right", drug.short=substring(cond.org,1,4), max_log_fold=ifelse(is.na(max_log_fold) | max_log_fold<(-2), -2, max_log_fold), value=pmin(1, abs(max_log_fold))) %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short, drug.known_activity, drug.uplc_excluded), by=c("drug.short"="drug.short")) %>%
    dplyr::mutate(drug.known_activity=dplyr::case_when(
      grepl("digo", drug.short) & !grepl("lenta", org) ~ NA_character_,
      TRUE                                             ~ drug.known_activity)) %>%
    dplyr::mutate(
      interaction=dplyr::case_when(
        max_interaction_type=="none"           ~ "No activity",
        !is.na(drug.uplc_excluded)             ~ "Excluded",
        grepl("Inhibits", drug.known_activity) ~ "Previously known",
        max_interaction_type=="lethal"         ~ "Growth Inhibition",
        max_log_fold<0                         ~ "Growth Inhibition", 
        max_log_fold>0                         ~ "Growth Promotion",
        TRUE                                   ~ "Unexpected Results (this shouldn't happen)")) %>%
    dplyr::select(source=drug.short, target=org, value, interaction, max_log_fold, side, max_interaction_type, max_pvalue, max_tech.pval)
  hits_growth = hits_growth.all %>% dplyr::filter(!interaction %in% c("No activity", "Excluded"))
  
  #
  # Build edges and nodes data.frames
  #
  edges.all = rbind(
    hits_degrad %>% dplyr::select(source, target, value, interaction, side) %>% data.frame(), 
    hits_growth %>% dplyr::select(source, target, value, interaction, side) %>% data.frame()) %>%
    dplyr::mutate(drug=ifelse(side=="left", target, source)) %>%
    dplyr::group_by(drug) %>%
    dplyr::do((function(z){
      if(!any(z$side=="left")) { z = rbind(z, data.frame(source="REMOVE1", target=z$drug[1], value=NA, interaction="REMOVE", side="left", drug=z$drug[1])) }
      if(!any(z$side=="right")) { z = rbind(z, data.frame(source=z$drug[1], target="REMOVE2", value=NA, interaction="REMOVE", side="right", drug=z$drug[1])) }
      z
    })(.))
  edges = edges.all %>% dplyr::filter(!(interaction %in% c("Untested", "No activity", "Excluded")) & !grepl("Mix", source)) %>% data.frame() # Hide everything unrelated to plot
  
  
  nodes = rbind(edges %>% dplyr::filter(side=="left") %>% dplyr::mutate(position="left") %>% dplyr::select(name=source, position), 
                edges %>% dplyr::filter(side=="right") %>% dplyr::mutate(position="right") %>% dplyr::select(name=target, position),
                edges %>% dplyr::mutate(position="middle") %>% dplyr::select(name=drug, position)) %>% unique() %>%
    dplyr::left_join(drug_map, by=c("name"="drug.short")) %>%
    dplyr::left_join(bug_map, by=c("name"="species.short")) %>%
    dplyr::left_join(bug_map, by=c("name"="species.long")) %>%
    dplyr::mutate(
      name.display=dplyr::case_when(!is.na(species.short)~species.short, !is.na(drug.long)~drug.long, T~name),
      group=dplyr::case_when(!is.na(species.group.x)~species.group.x, !is.na(species.group.y)~species.group.y, T~"drug"),
      order=dplyr::case_when(!is.na(drug.order)~drug.order, !is.na(species.group_order.x)~as.double(species.group_order.x+100), !is.na(species.group_order.y)~as.double(species.group_order.y+100))) %>%
    dplyr::arrange(order) %>%
    dplyr::mutate(id=1:n()-1) %>%
    dplyr::select(id, name, name.display, group, order)
  
  edges = edges %>% 
    dplyr::inner_join(nodes %>% dplyr::select(source.id=id, name, order.source=order), by=c("source"="name")) %>% 
    dplyr::inner_join(nodes %>% dplyr::select(target.id=id, name, order.target=order), by=c("target"="name")) %>%
    dplyr::arrange(side, source, target) %>% 
    dplyr::mutate(value=abs(value), interaction=ifelse(is.na(interaction), "Untested", gsub(" ", "_", interaction))) %>%
    dplyr::arrange(order.target, order.source) %>%
    dplyr::select(source, source.id, target, target.id, value, interaction) 
  
  #
  # Specify colors
  #
  species_colourScale = unique(bug_map[,c("species.group", "species.group_color")]) %>% dplyr::filter(!is.na(species.group_color)) %>% data.frame()
  colourScale = c(
    species_colourScale$species.group_color %>% setNames(species_colourScale$species.group),
    "Growth_Inhibition"="#F8523B", "Growth_Promotion"="#6CAA44", 
    "Biotransformation"="#C287C0", "Bioaccumulation"="#4483AA", 
    "Mix"="#FFFFFF", 
    "No_activity"="#FF8A11", "Untested"="#ABABAB", "Excluded"="#ABABAB", "Previously_known"="#AFCFCA",
    "drug"="#27303B", "Remove"="#FFFFFF", "REMOVE"="#FFFFFF"
  )
  colourScaleStr = networkD3::JS(paste0("d3.scaleOrdinal().domain(['", paste(names(colourScale), collapse="', '"), "']).range(['", paste(colourScale, collapse="', '"), "'])"))
  
  #
  # Plot legend for Sankey plot
  #
  pdf(file="reports/exp012_sankey_legend.pdf", paper="a4r")
  legend = cowplot::get_legend(ggplot(data.frame(Group=c("Biotransformation", "Bioaccumulation", "Untested", "Previously_known", "Growth_Inhibition", "Growth_Promotion", "Excluded", "Untested", "No_activity"))) + # names(colourScale)
                                 geom_bar(aes(x=Group, y=1, fill=Group), alpha=0.5, stat="identity") +
                                 scale_fill_manual(values=colourScale))
  grid::grid.newpage()
  grid::grid.draw(legend)
  dev.off()
  
  #
  # Plot sankey
  #
  networkD3::sankeyNetwork(
    Links=edges, 
    Nodes=nodes, 
    Source="source.id", 
    Target="target.id", 
    Value="value", 
    NodeID="name.display",
    NodeGroup="group",
    LinkGroup="interaction",
    fontSize = 20,
    iterations=0,
    colourScale=colourScaleStr)

  
  #
  # Summary (first degree)
  #
  hits_all = hits_degrad.all %>% 
    reshape2::dcast(Species.x ~ GroupName, value.var="replicates") %>%
    reshape2::melt(id.vars="Species.x", value.name="nreplicates", variable.name="drug.short") %>% 
    dplyr::left_join(hits_degrad.all %>% dplyr::select(Species.x, GroupName, nhits, MeanDiffToOwn, p_value.screen, p_adjust.screen), by=c("Species.x", "drug.short"="GroupName")) %>% 
    replace(is.na(.), 0) %>%
    dplyr::select(species.long=Species.x, drug.short, exp1.nhits=nhits, exp1.nreplicates=nreplicates, exp1.diff=MeanDiffToOwn, p_value.screen, p_adjust.screen) %>% 
    data.frame() %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short, drug.long, drug.known_activity, drug.excluded=drug.uplc_excluded), by=c("drug.short")) %>%
    dplyr::left_join(bug_map %>% dplyr::select(species.long, species.short), by=c("species.long")) %>%
    dplyr::left_join(hits_uplc %>% dplyr::select(drug.short, species.long, exp2.padjst_super=padjst.super, exp2.padjst_total=padjst.total, exp2.diff_super=MedianDiff.super, exp2.diff_total=MedianDiff.total, exp2.interaction=interaction), by=c("drug.short", "species.long")) %>%
    dplyr::left_join(hits_growth.all %>% dplyr::select(source, target, growth.effect=interaction, growth.pvalue=max_pvalue, growth.maxod_logfold=value), by=c("drug.short"="source", "species.short"="target")) %>%
    dplyr::filter(!grepl("Mix", species.short)) %>%
    dplyr::mutate(drug.known_activity=dplyr::case_when(
      drug.short=="leva" & grepl("infantis|bolteae|ramosum|comes", species.short) ~ "Depleted",
      drug.short=="digo" & !grepl("lenta", species.short)                         ~ NA_character_,
      T ~ drug.known_activity)) %>%
    dplyr::mutate(exp2.interaction=dplyr::case_when(
      !is.na(drug.known_activity)      ~ "Previously_known",
      !is.na(drug.excluded)            ~ "Excluded",
      is.na(exp2.interaction)          ~ "Untested",
      TRUE                             ~ exp2.interaction)) %>% 
    dplyr::mutate(
      exp1.is_hit=ifelse(!is.na(drug.known_activity) & exp1.nhits>0 | exp1.nhits/exp1.nreplicates>=0.5 | exp1.nhits>1, "Yes", "No"), 
      exp1.replicates=paste0(exp1.nhits, "/", exp1.nreplicates))  %>%
    dplyr::mutate(
      drug.known_activity=ifelse(drug.long=="Digoxin" & !grepl("lenta", species.short), NA_character_, drug.known_activity),
      drug.known_activity=dplyr::case_when(
        drug.long=="Levamisole" & species.short %in% c("B. longum subsp. infantis", "C. bolteae", "C. ramosum", "C. comes") ~ "Depleted",
        T ~ tidyr::replace_na(drug.known_activity, "No")),
      drug.excluded=tidyr::replace_na(drug.excluded, "No"),
      exp2.interaction=tidyr::replace_na(exp2.interaction, "Not tested"),
      exp1.is_hit=ifelse(drug.excluded=="No", exp1.is_hit, "Excluded"), 
      exp1.pvalue=round(p_value.screen, 5),
      exp1.padjst=round(p_adjust.screen, 5),
      exp1.diff=round(exp1.diff, 5), 
      exp2.padjst_super=round(exp2.padjst_super, 5), 
      exp2.padjst_total=round(exp2.padjst_total, 5), 
      exp2.diff_super=round(exp2.diff_super, 5), 
      exp2.diff_total=round(exp2.diff_total, 5),
      growth.effect=tidyr::replace_na(growth.effect, "Not tested"),
      growth.effect=ifelse(drug.excluded=="Yes" & growth.effect != "Not tested", "Excluded", growth.effect),
      growth.pvalue=round(growth.pvalue, 5),
      growth.maxod_logfold=round(growth.maxod_logfold, 5)
    ) %>%
    dplyr::select(species.short, drug.long, drug.known_activity, 
                  exp1.is_hit, exp2.interaction, growth.effect,                             # General HIT/NO-HIT information
                  exp1.diff, exp1.replicates, exp1.pvalue, exp1.padjst,                     # exp1
                  exp2.padjst_super, exp2.padjst_total, exp2.diff_super, exp2.diff_total,   # exp2
                  growth.pvalue, growth.maxod_logfold                                       # growth
    ) 
  
  
  #drug.excluded
  readr::write_tsv(hits_all, "reports/exp012_sankey_data.tsv", col_names=T, na="")  
  
  #
  # Hits overlap between UPLC depletion and growth
  #
  pdf(file="reports/exp012_growth2depletion_venn.pdf", paper="a4r")
  hits_all.venn = hits_all %>%
    dplyr::mutate(label=paste(drug.long,species.short))
  hits_all.venn=list(UPLC=hits_all.venn %>% dplyr::filter(exp1.is_hit=="Yes") %>% .$label , Growth=hits_all.venn %>% dplyr::filter(grepl("Growth", hits_all$growth.effect)) %>% .$label)
  hits_all.no = nrow(hits_all %>% dplyr::filter(exp1.is_hit=="No" & growth.effect=="No activity"))
  grid.newpage()
  grid.draw(VennDiagram::venn.diagram(hits_all.venn, filename=NULL, na="remove", force.unique=T, category.names=c("",""), fill=RColorBrewer::brewer.pal(3,  "Set1")[1:2], sub=paste0("No hit: ", hits_all.no)))
  dev.off()
  
  #
  # Statistics
  #
  hits_all.summary = hits_all %>% 
    reshape2::melt(measure.vars=c("drug.long", "species.short")) %>%
    dplyr::group_by(variable, value) %>% 
    dplyr::summarise(
      exp1=sum(grepl("Yes", exp1.is_hit)),
      exp2.bioaccumulation=sum(grepl("Yes", exp1.is_hit) & grepl("Bioaccumulation", exp2.interaction) & grepl("No", drug.known_activity)),
      exp2.biotransformation=sum(grepl("Yes", exp1.is_hit) & grepl("Biotransformation", exp2.interaction) & grepl("No", drug.known_activity)),
      exp2.known=sum(grepl("Yes", exp1.is_hit) & !grepl("No", drug.known_activity)),
      exp2.excluded=sum(grepl("Yes", exp1.is_hit) & grepl("No", drug.known_activity) & grepl("Excluded", exp2.interaction))) %>%
    data.frame()  %>%
    dplyr::arrange(variable, exp1)  
  readr::write_tsv(hits_all.summary, "reports/exp012_sankey_data_statistics.tsv", col_names=T, na="")  
}


