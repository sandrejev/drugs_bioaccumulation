library(rCharts)
library(networkD3)
library(dplyr)
library(readr)
library(xlsx)
library(ggplot2)

#
# t-test
#
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) {
    # welch-satterthwaite equation
    se = sqrt( (s1^2/n1) + (s2^2/n2) )
    df = ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else {
    # pooled standard deviation, scaled by the sample sizes
    se = sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df = n1+n2-2
  }      
  
  t = (m1-m2-m0)/se 
  dat = data.frame("Difference of means"=m1-m2, "Std Error"=se, "t"=t, "p_value"=2*pt(-abs(t),df))
  return(dat) 
}

ttest.degrad.batch = function(data.degrad){
  meanmax = data.degrad %>%
    dplyr::group_by(Plate, GroupName) %>%
    dplyr::summarise(Mean.Max=mean(maxOD))
  
  ttest.diffSD = data.degrad %>%
    dplyr::group_by(GroupName,Plate,Batch,Species.x) %>%
    dplyr::do(t.test2(mean(.$DrugRatio), .$BatchMean[1], sd(.$DrugRatio),.$SD.ctrls[1], length(.$DrugRatio),.$N.ctrls[1], equal.variance=T)) %>%
    data.frame() %>%
    dplyr::mutate(p_adjust=p.adjust(.$p_value,method="BH"))
  
  ttest.sameSD = data.degrad %>%
    dplyr::group_by(GroupName,Plate,Batch,Species.x) %>%
    dplyr::do(t.test2(mean(.$DrugRatio),.$BatchMean[1], .$SD.ctrls[1], .$SD.ctrls[1], length(.$DrugRatio), .$N.ctrls[1], equal.variance=T)) %>%
    data.frame() %>%
    dplyr::mutate(p_adjust=p.adjust(.$p_value,method="BH")) %>%
    dplyr::inner_join(unique(data.degrad %>% dplyr::select(GroupName, Plate, BatchMean, growth)), by=c("GroupName", "Plate")) %>%
    dplyr::mutate(DiffOfMeans = Difference.of.means/BatchMean) %>%
    dplyr::inner_join(meanmax, by=c("GroupName", "Plate")) %>%
    dplyr::mutate(Norm.Max=DiffOfMeans/Mean.Max)
  
  data.get.new = list(ttest.sameSD, ttest.diffSD)
}

#effect = c(      1, 0.55553, 0.31372,      1, 0.32313, 0.38565, 0.3732, 0.45363, 0.37436, 0.31253, 0.39965, 0.37689, 0.6261, 0.31382, 0.31596, 0.55628, 0.36803, 0.74794, 0.47537, 0.62475, 0.33173, 0.45039, 0.8063, 0.87386, 1.0955, 0.98689, 0.68736, 0.6171, 0.98691, 0.87599, 0.44559, 0.31173, 0.30683, 0.69593,      1, 0.94988, 0.37964, 0.67791,      1,      1, 0.71145, 0.47181, 0.46187,      1,      1,      1, 0.72098,      1, 0.31779, 0.45687, 0.33309, 0.49094, 0.43271, 0.37038, 0.44246, 0.39824,      1, 0.31065, 0.33679, 0.37264, 0.31704, 0.31634, 0.34415, 0.35284, 0.30317, 0.35411, 0.74845, 0.69957, 0.64449, 0.46201, 0.44586, 0.31381, 0.74337, 0.41115,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1, 0.37649, 0.2292, 0.3005, 0.27439, 0.76435, 0.12911, 0.45192,      2, 0.11446,      2, 0.30464,      2,      2, 0.74655,      2,      2,      2,      2,      2,      2,      2,      2,      2,      2,      2, 0.16781, 0.41567, 0.36545, 0.13614, 0.29798, 1.2398, 0.11057, 0.25022)
#source = c("Escherichia coli ED1a", "Fusobacterium nucleatum nucleatum", "Escherichia coli iAi1", "Eggerthella lenta", "Bifidobacterium longum longum", "Bifidobacterium longum infantis", "Bacteroides uniformis", "Bacteroides vulgatus", "Coprococcus comes", "Clostridium ramosum", "Lactobacillus plantarum", "Bifidobacterium animalis lactis", "Clostridium ramosum", "Bifidobacterium animalis lactis", "Bifidobacterium longum longum", "Bifidobacterium longum infantis", "Bacteroides thetaiotaomicron", "Bacteroides uniformis HM716", "Bacteroides uniformis", "Coprococcus comes", "Clostridium ramosum", "Clostridium saccharolyticum", "Eubacterium rectale", "Fusobacterium nucleatum nucleatum", "Bacteroides thetaiotaomicron", "Bacteroides uniformis HM715", "Bacteroides uniformis", "Lactobacillus gasseri", "Bifidobacterium animalis lactis", "Bifidobacterium longum infantis", "Bacteroides thetaiotaomicron", "Bacteroides uniformis HM715", "Bacteroides uniformis HM716", "Bacteroides uniformis", "Bacteroides vulgatus", "Clostridium bolteae", "Coprococcus comes", "Clostridium ramosum", "Escherichia coli ED1a", "Escherichia coli iAi1", "Eggerthella lenta", "Eubacterium rectale", "Fusobacterium nucleatum nucleatum", "Lactococcus lactis", "Lactobacillus paracasei", "Lactobacillus plantarum", "Ruminococcus torques", "Streptococcus salivarius", "Bifidobacterium animalis lactis", "Bifidobacterium longum infantis", "Bacteroides uniformis HM715", "Clostridium bolteae", "Coprococcus comes", "Fusobacterium nucleatum nucleatum", "Streptococcus salivarius", "Escherichia coli iAi1", "Eggerthella lenta", "Lactobacillus gasseri", "Bifidobacterium longum infantis", "Bacteroides vulgatus", "Bifidobacterium animalis lactis", "Bacteroides fragilis", "Bifidobacterium longum longum", "Bacteroides thetaiotaomicron", "Coprococcus comes", "Clostridium ramosum", "Bacteroides thetaiotaomicron", "Bacteroides vulgatus", "Clostridium ramosum", "Clostridium saccharolyticum", "Eggerthella lenta", "Eubacterium rectale", "Lactobacillus plantarum", "Streptococcus salivarius", "Bifidobacterium animalis lactis", "Bacteroides fragilis", "Bifidobacterium longum longum", "Bifidobacterium longum infantis", "Bacteroides thetaiotaomicron", "Bacteroides uniformis HM715", "Bacteroides uniformis HM716", "Bacteroides uniformis", "Bacteroides vulgatus", "Clostridium bolteae", "Coprococcus comes", "Clostridium ramosum", "Clostridium saccharolyticum", "Escherichia coli iAi1", "Eggerthella lenta", "Eubacterium rectale", "Fusobacterium nucleatum nucleatum", "Lactobacillus gasseri", "Lactococcus lactis", "Lactobacillus plantarum", "Ruminococcus torques", "Streptococcus salivarius", "Aripiprazole", "Aripiprazole", "Digoxin", "Digoxin", "Digoxin", "Duloxetine", "Duloxetine", "Duloxetine", "Levamisole", "Loperamide", "Loperamide", "Loperamide", "Loperamide", "Loperamide", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Montelukast", "Ranitidine", "Rosiglitazone", "Rosuvastatin", "Sulfasalazine", "Sulfasalazine", "Tolmetin", "Tolmetin")
#target = c("Acetaminophen", "Acetaminophen", "Aripiprazole", "Digoxin", "Duloxetine", "Duloxetine", "Duloxetine", "Duloxetine", "Duloxetine", "Duloxetine", "Duloxetine", "Ezetimibe", "Ezetimibe", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Levamisole", "Metformin", "Metformin", "Metformin", "Metformin", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Metronidazole", "Montelukast", "Montelukast", "Montelukast", "Montelukast", "Montelukast", "Montelukast", "Montelukast", "Ranitidine", "Ranitidine", "Ranitidine", "Roflumilast", "Roflumilast", "Rosiglitazone", "Rosiglitazone", "Rosiglitazone", "Rosiglitazone", "Rosiglitazone", "Rosiglitazone", "Simvastatin", "Simvastatin", "Simvastatin", "Simvastatin", "Simvastatin", "Simvastatin", "Simvastatin", "Simvastatin", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "Sulfasalazine", "E. coli IAI1", "L. gasseri", "E. coli IAI1", "E. lenta", "R. torques", "C. saccharolyticum", "E. coli IAI1", "E. rectale", "L. lactis", "B. longum subsp. infantis", "E. coli IAI1", "E. lenta", "E. rectale", "L. lactis", "B. fragilis", "B. longum subsp. infantis", "B. thetaiotaomicron", "B. uniformis", "B. uniformis HM715", "B. uniformis HM716", "C. ramosum", "C. saccharolyticum", "E. lenta", "E. rectale", "F. nucleatum subsp. nucleatum", "B. uniformis HM715", "E. coli IAI1", "E. coli IAI1", "B. vulgatus", "C. ramosum", "R. torques", "L. plantarum", "R. torques")
#drugs = c("Loperamide", "Tolmetin", "Donepezil", "Roflumilast", "Ezetimibe", "Rosiglitazone","Ranitidine", "Levamisole", "Montelukast", "Simvastatin", "Digoxin", "Aripiprazole", "Duloxetine", "Sulfasalazine", "Metronidazole")
#species_depletors = c("Fusobacterium nucleatum nucleatum", "Eggerthella lenta", "Bifidobacterium longum infantis", "Bifidobacterium animalis lactis", "Clostridium bolteae", "Clostridium saccharolyticum", "Clostridium ramosum", "Coprococcus comes", "Ruminococcus gnavus", "Escherichia coli iAi1", "Escherichia coli ED1a", "Lactococcus lactis", "Lactobacillus plantarum", "Lactobacillus gasseri", "Lactobacillus paracasei", "Streptococcus salivarius", "Bacteroides uniformis", "Bacteroides uniformis HM715", "Bacteroides uniformis HM716", "Bacteroides thetaiotaomicron", "Bacteroides vulgatus", "Mix Depletion", "Mix No Depletion", "REMOVE1")
#species_influeced = c("L. lactis", "L. plantarum", "L. gasseri", "E. lenta", "B. longum subsp. infantis", "E. coli IAI1", "F. nucleatum subsp. nucleatum", "E. rectale", "R. torques", "C. saccharolyticum", "C. ramosum", "B. uniformis HM716", "B. uniformis HM715", "B. thetaiotaomicron", "B. uniformis", "B. fragilis", "REMOVE2")

drug_map = readr::read_delim("data/drug_map.tsv", "\t")
bug_map = readr::read_delim("data/bug_map.tsv", "\t")

species.groups = c(
  "Fusobacterium nucleatum nucleatum"="Fusobacteria", "F. nucleatum subsp. nucleatum"="Fusobacteria",
  "Eggerthella lenta"="Actinobacteria", "Bifidobacterium animalis lactis"="Actinobacteria", "Bifidobacterium longum infantis"="Actinobacteria", "B. longum subsp. infantis"="Actinobacteria", "B. animalis subsp. lactis BI-07"="Actinobacteria", "E. lenta"="Actinobacteria",
  "Clostridium bolteae"="Firmicutes", "Clostridium saccharolyticum"="Firmicutes", "Clostridium ramosum"="Firmicutes", "Coprococcus comes"="Firmicutes", "Ruminococcus gnavus"="Firmicutes", "R. torques"="Firmicutes", "C. saccharolyticum"="Firmicutes", "E. rectale"="Firmicutes", "C. ramosum"="Firmicutes",
  "Escherichia coli ED1a"="Proteobacteria", "Escherichia coli iAi1"="Proteobacteria", "E. coli IAI1"="Proteobacteria",
  "Lactococcus lactis"="Firmicutes", "Lactobacillus plantarum"="Firmicutes", "L. plantarum"="Firmicutes", "Lactobacillus gasseri"="Firmicutes", "Lactobacillus paracasei"="Firmicutes", "Streptococcus salivarius"="Firmicutes", "L. gasseri"="Firmicutes", "L. lactis"="Firmicutes",
  "Bacteroides uniformis"="Bacteroidetes", "Bacteroides uniformis HM715"="Bacteroidetes", "Bacteroides uniformis HM716"="Bacteroidetes", "Bacteroides thetaiotaomicron"="Bacteroidetes", "Bacteroides vulgatus"="Bacteroidetes", "B. fragilis"="Bacteroidetes", "B. thetaiotaomicron"="Bacteroidetes", "B. uniformis"="Bacteroidetes", "B. uniformis HM715"="Bacteroidetes", "B. uniformis HM716"="Bacteroidetes", "B. vulgatus"="Bacteroidetes",
  "Mix Degrad"="Mix", "Mix No"="Mix", "REMOVE1"="remove", "REMOVE2"="remove")
species.groups = data.frame(group=species.groups) %>% 
  tibble::rownames_to_column("name") %>%
  dplyr::mutate(group.order=as.numeric(factor(group, unique(as.character(group))))) %>%
  data.frame()



load("data/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
plates_with_growth = unique(data.clean %>% dplyr::filter(Status=="GMM" & growth=="growth") %>% .$Plate)
data.degrad = data.clean %>% dplyr::filter(Status=="sample" & growth=="growth" & dummy==1 & Plate %in% plates_with_growth)

data.get <- ttest.degrad.batch(data.degrad) 
hits = data.get[[1]] %>%
  dplyr::inner_join(unique(data.clean[,c("GroupName","Plate","DiffToOwn")])) %>%
  dplyr::filter(-0.3>DiffToOwn & Plate!="CramosumA2" & GroupName!="metf") %>%
  dplyr::group_by(GroupName, Species.x) %>%
  dplyr::summarize(nhits=length(GroupName), MeanDiffToOwn=mean(DiffToOwn), side="left") %>%
  data.frame() %>%
  dplyr::mutate(target=as.character(GroupName), source=as.character(Species.x), value=abs(MeanDiffToOwn)) %>%
  dplyr::filter(nhits>1) 


# Only Rosu -> B. vulgatus doesn't appear in thesis
# From 160109_DataGrowthOwn.R
hits_growth <- read.csv("data/curves.rel_annotation_2016-11-28.tab",sep="\t", header=T) %>%
  dplyr::filter(max_interaction_type!="none" & !grepl("pyri", cond.org)) %>%
  dplyr::mutate(side="right", GroupName=substring(cond.org,1,4), max_log_fold=ifelse(is.na(max_log_fold) | max_log_fold<(-2), -2, max_log_fold), value=pmin(1, abs(max_log_fold))) %>%
  dplyr::mutate(interaction=dplyr::case_when(
    GroupName=="metr" ~ "Previously known", 
    max_interaction_type=="lethal" ~ "Growth Inhibition",
    max_log_fold<0 ~ "Restrained Growth", 
    max_log_fold>0 ~ "Growth Promotion")) %>%
  dplyr::select(source=GroupName, target=org, value, interaction, max_log_fold, side, max_interaction_type, max_pvalue, max_tech.pval)


edges = rbind(
  hits %>% dplyr::mutate(interaction="known") %>% dplyr::select(source, target, value, interaction, side) %>% data.frame(), 
  hits_growth %>% dplyr::select(source, target, value, interaction=max_interaction_type, side) %>% data.frame()) %>%
  dplyr::mutate(drug=ifelse(side=="left", target, source)) %>%
  dplyr::group_by(drug) %>%
  dplyr::do((function(z){
    if(!any(z$side=="left")) { z = rbind(z, data.frame(source="REMOVE1", target=z$drug[1], value=NA, interaction="REMOVE", side="left", drug=z$drug[1])) }
    if(!any(z$side=="right")) { z = rbind(z, data.frame(source=z$drug[1], target="REMOVE2", value=NA, interaction="REMOVE", side="right", drug=z$drug[1])) }
    z
  })(.)) %>% data.frame() %>%
  dplyr::left_join(species.groups %>% dplyr::rename(group.source="group", group.order.source="group.order"), by=c("source"="name")) %>%
  dplyr::left_join(species.groups %>% dplyr::rename(group.target="group", group.order.target="group.order"), by=c("target"="name")) %>%
  dplyr::mutate(group.node=dplyr::case_when(!is.na(group.source)~group.source, !is.na(group.target)~group.target)) %>%
  dplyr::mutate(group.order.node=dplyr::case_when(!is.na(group.order.source)~group.order.source, !is.na(group.order.target)~group.order.target)) %>%
  dplyr::mutate(group.source=NULL, group.target=NULL) %>%
  dplyr::arrange(side, group.order.node, source, target) %>%
  data.frame()

nodes = rbind(edges %>% dplyr::filter(side=="left") %>% dplyr::mutate(position="left") %>% dplyr::select(name=source, position, group=group.node),
              edges %>% dplyr::filter(side=="right") %>% dplyr::mutate(position="right") %>% dplyr::select(name=target, position, group=group.node),
              edges %>% dplyr::mutate(position="middle", group.node="drug") %>% dplyr::left_join(drug_map, by=c("drug"="drug.short")) %>% dplyr::arrange(drug.order) %>% dplyr::select(name=drug, position, group=group.node)) %>% 
  unique() %>%
  dplyr::left_join(drug_map, by=c("name"="drug.short")) %>%
  dplyr::left_join(bug_map, by=c("name"="species.short")) %>%
  dplyr::mutate(id=0:(n()-1), name.long=dplyr::case_when(!is.na(species.long)~species.long, !is.na(drug.long)~drug.long, T~name)) #group=as.character(group), 

edges = edges %>% 
  dplyr::inner_join(nodes %>% dplyr::select(source.id=id, name), by=c("source"="name")) %>% 
  dplyr::inner_join(nodes %>% dplyr::select(target.id=id, name), by=c("target"="name")) %>%
  dplyr::arrange(side, source, target) %>% 
  dplyr::mutate(value=abs(value)) %>%
  dplyr::select(source, source.id, target, target.id, value, group.node)

networkD3::sankeyNetwork(
  Links=edges, 
  Nodes=nodes, 
  Source="source.id", 
  Target="target.id", 
  Value="value", 
  NodeID="name.long",
  NodeGroup="group",
  iterations=0)
