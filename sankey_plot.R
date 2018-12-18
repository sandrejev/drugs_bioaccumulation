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

plot.sankey = function()
{
  drug_map = readr::read_delim("data/drug_map.tsv", "\t")
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")
  
  load("data/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
  plates_with_growth = unique(data.clean %>% dplyr::filter(Status=="GMM" & growth=="growth") %>% .$Plate)
  data.degrad = data.clean %>% dplyr::filter(Status=="sample" & growth=="growth" & dummy==1 & Plate %in% plates_with_growth)
  
  data.get <- ttest.degrad.batch(data.degrad) 
  hits_degrad = data.get[[1]] %>%
    dplyr::inner_join(data.clean %>% dplyr::select(GroupName, Plate, DiffToOwn) %>% unique(), by=c("GroupName", "Plate")) %>%
    dplyr::filter(-0.3>DiffToOwn & Plate!="CramosumA2" & GroupName!="metf") %>%
    dplyr::group_by(GroupName, Species.x) %>%
    dplyr::summarize(nhits=length(GroupName), MeanDiffToOwn=mean(DiffToOwn), side="left") %>%
    data.frame() %>%
    dplyr::mutate(target=as.character(GroupName), source=as.character(Species.x), value=abs(MeanDiffToOwn)) %>%
    dplyr::filter(nhits>1) 
  
  # Only Rosu -> B. vulgatus doesn't appear in thesis
  # From 160109_DataGrowthOwn.R
  hits_growth = read.csv("data/curves.rel_annotation_2016-11-28.tab",sep="\t", header=T) %>%
    dplyr::filter(max_interaction_type!="none" & !grepl("pyri", cond.org)) %>%
    dplyr::mutate(side="right", GroupName=substring(cond.org,1,4), max_log_fold=ifelse(is.na(max_log_fold) | max_log_fold<(-2), -2, max_log_fold), value=pmin(1, abs(max_log_fold))) %>%
    dplyr::mutate(interaction=dplyr::case_when(
      GroupName=="metr" ~ "Previously known", 
      max_interaction_type=="lethal" ~ "Growth Inhibition",
      max_log_fold<0 ~ "Restrained Growth", 
      max_log_fold>0 ~ "Growth Promotion")) %>%
    dplyr::select(source=GroupName, target=org, value, interaction, max_log_fold, side, max_interaction_type, max_pvalue, max_tech.pval)
  
  #
  # Build edges and nodes data.frames
  #
  edges = rbind(
    hits_degrad %>% dplyr::mutate(interaction="known") %>% dplyr::select(source, target, value, interaction, side) %>% data.frame(), 
    hits_growth %>% dplyr::select(source, target, value, interaction=max_interaction_type, side) %>% data.frame()) %>%
    dplyr::mutate(drug=ifelse(side=="left", target, source)) %>%
    dplyr::group_by(drug) %>%
    dplyr::do((function(z){
      if(!any(z$side=="left")) { z = rbind(z, data.frame(source="REMOVE1", target=z$drug[1], value=NA, interaction="REMOVE", side="left", drug=z$drug[1])) }
      if(!any(z$side=="right")) { z = rbind(z, data.frame(source=z$drug[1], target="REMOVE2", value=NA, interaction="REMOVE", side="right", drug=z$drug[1])) }
      z
    })(.)) %>% data.frame() 
  
  nodes = rbind(edges %>% dplyr::filter(side=="left") %>% dplyr::mutate(position="left") %>% dplyr::select(name=source, position), 
                edges %>% dplyr::filter(side=="right") %>% dplyr::mutate(position="right") %>% dplyr::select(name=target, position),
                edges %>% dplyr::mutate(position="middle") %>% dplyr::select(name=drug, position)) %>% unique() %>%
    dplyr::left_join(drug_map, by=c("name"="drug.short")) %>%
    dplyr::left_join(bug_map, by=c("name"="species.short")) %>%
    dplyr::left_join(bug_map, by=c("name"="species.long")) %>%
    dplyr::mutate(
      name.long=dplyr::case_when(!is.na(species.long)~species.long, !is.na(drug.long)~drug.long, T~name),
      name.short=dplyr::case_when(!is.na(species.short)~species.short, T~name),
      group=dplyr::case_when(!is.na(species.group.x)~species.group.x, !is.na(species.group.y)~species.group.y, T~"drug"),
      order=dplyr::case_when(!is.na(drug.order)~drug.order, !is.na(species.group_order.x)~species.group_order.x, !is.na(species.group_order.y)~species.group_order.y)) %>%
    dplyr::arrange(order) %>%
    dplyr::mutate(id=1:n()-1) %>%
    dplyr::select(id, name, name.long, name.short, group, order)
  
  edges = edges %>% 
    dplyr::inner_join(nodes %>% dplyr::select(source.id=id, name), by=c("source"="name")) %>% 
    dplyr::inner_join(nodes %>% dplyr::select(target.id=id, name), by=c("target"="name")) %>%
    dplyr::arrange(side, source, target) %>% 
    dplyr::mutate(value=abs(value)) %>%
    dplyr::select(source, source.id, target, target.id, value)
  
  networkD3::sankeyNetwork(
    Links=edges, 
    Nodes=nodes, 
    Source="source.id", 
    Target="target.id", 
    Value="value", 
    NodeID="name.long",
    NodeGroup="group",
    iterations=0)
}