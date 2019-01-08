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

#hits_uplc = readr::read_delim("data/assay.hits.csv", ";")

plot.sankey = function()
{
  drug_map = readr::read_delim("data/drug_map.tsv", "\t")
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")
  
  #
  # Load depletion experiment (#1) data for LEFT side of the sankey plot
  #
  load("data/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
  plates_with_growth = unique(data.clean %>% dplyr::filter(Status=="GMM" & growth=="growth") %>% .$Plate)
  data.degrad = data.clean %>% dplyr::filter(Status=="sample" & growth=="growth" & dummy==1 & Plate %in% plates_with_growth)
  data.get = ttest.degrad.batch(data.degrad) 
  
  #
  # Calculate edges for LEFT side of the sankey plot and join with UPLC experiment (#2)
  #
  hits_uplc = readr::read_delim("data/degradation_hits.tsv", "\t")
  hits_degrad.all = data.get[[1]] %>%
    dplyr::mutate(Species.x=tidyr::replace_na(Species.x, ""), GroupName=tidyr::replace_na(as.character(GroupName), "")) %>%
    dplyr::left_join(data.clean %>% dplyr::select(GroupName, Plate, DiffToOwn) %>% unique(), by=c("GroupName", "Plate")) %>%
    dplyr::group_by(GroupName, Species.x) %>%
    dplyr::summarize(replicates=length(GroupName), nhits=sum(-0.3>DiffToOwn), MeanDiffToOwn=mean(DiffToOwn[-0.3>DiffToOwn]), side="left") %>%
    data.frame()
  hits_degrad = hits_degrad.all %>%
    dplyr::mutate(target=as.character(GroupName), source=as.character(Species.x), value=abs(MeanDiffToOwn)) %>%
    dplyr::filter(nhits>1) %>% # TODO: Why does it have to be more than one
    dplyr::left_join(hits_uplc %>% dplyr::select(species.long, drug.short, interaction), by=c("source"="species.long", "target"="drug.short")) %>%
    dplyr::inner_join(drug_map %>% dplyr::select(drug.short, drug.known_activity, drug.uplc_excluded), by=c("target"="drug.short")) %>%
    dplyr::mutate(interaction=ifelse(grepl("Depleted", drug.known_activity), "Previously_known", interaction)) %>%
    dplyr::mutate(interaction=ifelse(is.na(interaction), "Untested", interaction)) %>%
    dplyr::mutate(interaction=ifelse(!is.na(drug.uplc_excluded), "Excluded", interaction))
  
  
  #
  # Load growth effect experiment (#0) for RIGHT side of the sankey plot
  #
  hits_growth.all = read.csv("data/curves.rel_annotation_2016-11-28.tab",sep="\t", header=T) %>%
    dplyr::filter(!grepl("pyri", cond.org)) %>% # Pyri is not appearing anywhere in the drug list
    dplyr::mutate(side="right", drug.short=substring(cond.org,1,4), max_log_fold=ifelse(is.na(max_log_fold) | max_log_fold<(-2), -2, max_log_fold), value=pmin(1, abs(max_log_fold))) %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short, drug.known_activity), by=c("drug.short"="drug.short")) %>%
    dplyr::mutate(
      interaction=dplyr::case_when(
        max_interaction_type=="none" ~ "No activity",
        grepl("Inhibits", drug.known_activity) ~ "Previously known",
        max_interaction_type=="lethal" ~ "Growth Inhibition",
        max_log_fold<0 ~ "Growth Inhibition", 
        max_log_fold>0 ~ "Growth Promotion",
        T ~ "Unexpected Results (this shouldn't happen)")) %>%
    dplyr::select(source=drug.short, target=org, value, interaction, max_log_fold, side, max_interaction_type, max_pvalue, max_tech.pval)
  hits_growth = hits_growth.all %>% dplyr::filter(interaction != "No activity")
  
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
  edges = edges.all %>% dplyr::filter(interaction!="No activity" & !grepl("Mix", source) & interaction!="Excluded") %>% data.frame() # Hide everything unrelated to plot
  
  nodes = rbind(edges %>% dplyr::filter(side=="left") %>% dplyr::mutate(position="left") %>% dplyr::select(name=source, position), 
                edges %>% dplyr::filter(side=="right") %>% dplyr::mutate(position="right") %>% dplyr::select(name=target, position),
                edges %>% dplyr::mutate(position="middle") %>% dplyr::select(name=drug, position)) %>% unique() %>%
    dplyr::left_join(drug_map, by=c("name"="drug.short")) %>%
    dplyr::left_join(bug_map, by=c("name"="species.short")) %>%
    dplyr::left_join(bug_map, by=c("name"="species.long")) %>%
    dplyr::mutate(
      name.display=dplyr::case_when(!is.na(species.short)~species.short, !is.na(drug.long)~drug.long, T~name),
      group=dplyr::case_when(!is.na(species.group.x)~species.group.x, !is.na(species.group.y)~species.group.y, T~"drug"),
      order=dplyr::case_when(!is.na(drug.order)~drug.order, !is.na(species.group_order.x)~as.integer(species.group_order.x+100), !is.na(species.group_order.y)~as.integer(species.group_order.y+100))) %>%
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
  colourScale = c(
    "Growth_Inhibition"="#F8523B", "Growth_Promotion"="#6CAA44", 
    "Biotransformation"="#C287C0", "Bioaccumulation"="#4483AA", "Fusobacteria"="#DAC300", 
    "Actinobacteria"="#825BB8", "Firmicutes"="#70A8FF", "Proteobacteria"="#837E79", "Bacteroidetes"="#69BB6F", 
    "Mix"="#FFFFFF", 
    "No_activity"="#FF8A11", "Untested"="#ABABAB", "Excluded"="#ABABAB", "Previously_known"="#AFCFCA",
    "drug"="#27303B", "Remove"="#FFFFFF", "REMOVE"="#FFFFFF"
  )
  colourScaleStr = networkD3::JS(paste0("d3.scaleOrdinal().domain(['", paste(names(colourScale), collapse="', '"), "']).range(['", paste(colourScale, collapse="', '"), "'])"))
  
  #
  # Plot legend for Sankey plot
  #
  legend = cowplot::get_legend(ggplot(data.frame(Group=c("Biotransformation", "Bioaccumulation", "Untested", "Previously_known", "Growth_Inhibition", "Growth_Promotion", "Excluded", "Untested", "No_activity"))) + # names(colourScale)
    geom_bar(aes(x=Group, y=1, fill=Group), alpha=0.5, stat="identity") +
    scale_fill_manual(values=colourScale))
  grid.newpage()
  grid.draw(legend) 
  
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
    dplyr::left_join(hits_degrad.all %>% dplyr::select(Species.x, GroupName, nhits, MeanDiffToOwn), by=c("Species.x", "drug.short"="GroupName")) %>% 
    replace(is.na(.), 0) %>%
    dplyr::mutate(exp1.is_hit=ifelse(nhits>1, "Yes", "No"), exp1.replicates=paste0(nhits, "/", nreplicates))  %>%
    dplyr::select(species.long=Species.x, drug.short, exp1.is_hit, exp1.replicates, exp1.diff=MeanDiffToOwn) %>% 
    data.frame() %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short, drug.long, drug.known_activity, drug.excluded=drug.uplc_excluded), by=c("drug.short")) %>%
    dplyr::left_join(bug_map %>% dplyr::select(species.long, species.short), by=c("species.long")) %>%
    dplyr::left_join(hits_uplc %>% dplyr::select(drug.short, species.long, exp2.padjst_super=padjst.super, exp2.padjst_total=padjst.total, exp2.diff_super=MedianDiff.super, exp2.diff_total=MedianDiff.total, exp2.interaction=interaction), by=c("drug.short", "species.long")) %>%
    dplyr::left_join(hits_growth.all %>% dplyr::select(source, target, growth.effect=interaction, growth.pvalue=max_pvalue, growth.maxod_logfold=value), by=c("drug.short"="source", "species.short"="target")) %>%
    dplyr::filter(!grepl("Mix", species.short)) %>%
    dplyr::mutate(
      drug.known_activity=tidyr::replace_na(drug.known_activity, "No"),
      drug.excluded=tidyr::replace_na(drug.excluded, "No"),
      exp2.interaction=tidyr::replace_na(exp2.interaction, "Not tested"),
      exp1.is_hit=ifelse(drug.excluded=="No", exp1.is_hit, "Excluded"), 
      exp1.diff=round(exp1.diff, 3), 
      exp2.padjst_super=round(exp2.padjst_super, 3), 
      exp2.padjst_total=round(exp2.padjst_total, 3), 
      exp2.diff_super=round(exp2.diff_super, 3), 
      exp2.diff_total=round(exp2.diff_total, 3),
      growth.effect=tidyr::replace_na(growth.effect, "Not tested"),
      growth.pvalue=round(growth.pvalue, 3),
      growth.maxod_logfold=round(growth.maxod_logfold, 3)
    ) %>%
    dplyr::select(species.short, drug.long, drug.excluded, drug.known_activity, 
                  exp1.is_hit, exp2.interaction, growth.effect, exp1.diff, exp1.replicates, 
                  exp2.padjst_super, exp2.padjst_total, exp2.diff_super, exp2.diff_total,
                  growth.pvalue, growth.maxod_logfold) 
  readr::write_tsv(hits_all, "reports/data_sankey.tsv", col_names=T, na="")  
  
  #
  # Summary (second degree)
  #
  for(gr in c("drug.long", "species.short")) {
    hits_all.summary.1 = hits_all %>% 
      dplyr::rename(group=gr) %>%
      dplyr::mutate(exp1.is_hit=paste0("exp1_hit.", exp1.is_hit)) %>%
      dplyr::group_by(group, exp1.is_hit) %>% 
      dplyr::summarise(exp1.count=length(exp1.is_hit)) %>%
      reshape2::dcast(group ~ exp1.is_hit, value.var="exp1.count") %>%
      dplyr::mutate(exp1_total=rowSums(.[,c("exp1_hit.Excluded", "exp1_hit.No", "exp1_hit.Yes")], na.rm=T)) %>%
      data.frame()
    
    hits_all.summary.2 = hits_all %>% 
      dplyr::rename(group=gr) %>%
      dplyr::mutate(exp2.is_hit=paste0("exp2_hit.", exp2.interaction)) %>%
      dplyr::group_by(group, exp2.is_hit) %>% 
      dplyr::summarise(exp2.count=length(exp2.is_hit)) %>%
      reshape2::dcast(group ~ exp2.is_hit, value.var="exp2.count") %>%
      dplyr::mutate(exp2_total=rowSums(.[,c("exp2_hit.Bioaccumulation", "exp2_hit.Biotransformation", "exp2_hit.Excluded", "exp2_hit.No activity")], na.rm=T))%>%
      data.frame()
    
    hits_all.summary.3 = hits_all %>% 
      dplyr::rename(group=gr) %>%
      dplyr::mutate(growth.effect=paste0("growth.", growth.effect)) %>%
      dplyr::group_by(group, growth.effect) %>% 
      dplyr::summarise(growth.count=length(growth.effect)) %>%
      reshape2::dcast(group ~ growth.effect, value.var="growth.count") %>%
      dplyr::mutate(growth_total=rowSums(.[,c("growth.Growth Inhibition", "growth.Growth Promotion", "growth.No activity", "growth.Previously known")], na.rm=T)) %>%
      data.frame()
    
    hits_all.summary = hits_all.summary.1 %>%
      dplyr::left_join(hits_all.summary.2, by="group") %>%
      dplyr::left_join(hits_all.summary.3, by="group") %>% 
      replace(is.na(.), 0)
    hits_all.summary = rbind(hits_all.summary, c("Total", colSums(data.matrix(hits_all.summary[,-1]), na.rm=T)))
    print(head(hits_all.summary))
    
    readr::write_tsv(hits_all.summary, paste0("reports/data_sankey_", gr, ".tsv"), col_names=T, na="")  
  }
  
  
  
  
  
  #
  # Statistics on all edges
  #
  edges.all %>%
    dplyr::filter(interaction!="REMOVE") %>%
    dplyr::group_by(side) %>%
    dplyr::mutate(interactions_sum=length(side)) %>%
    data.frame() %>%
    dplyr::group_by(side, interaction) %>%
    dplyr::summarise(count=length(interaction), proportion=paste0(round(100*count/interactions_sum[1]), "%")) %>%
    data.frame()
  
  
  
  # Experiment #2 summary
  hits_uplc.drugs_summary = hits_uplc %>% 
    dplyr::group_by(drug.long, interaction) %>%
    dplyr::summarise(species=length(species.long)) %>%
    reshape2::dcast(drug.long ~ interaction, value.var="species") %>%
    replace(is.na(.), 0) %>%
    data.frame() %>%
    dplyr::mutate(Total=rowSums(.[,-1]))
  hits_uplc.drugs_summary = rbind(hits_uplc.drugs_summary, c(drug.long="Total", colSums(hits_uplc.drugs_summary[,-1])))
  writeLines("==============================================\n Experiment #2 (transformation type) summary: \n==============================================")
  hits_uplc.drugs_summary
  
  hits_uplc.species_summary = hits_uplc %>% 
    dplyr::group_by(species.short, interaction) %>%
    dplyr::summarise(species=length(drug.long)) %>%
    reshape2::dcast(species.short ~ interaction, value.var="species") %>%
    replace(is.na(.), 0) %>%
    data.frame() %>%
    dplyr::mutate(Total=rowSums(.[,-1]), species.short=ifelse(species.short==0, NA, species.short))
  hits_uplc.species_summary = rbind(hits_uplc.species_summary, c(species.short="Total", colSums(hits_uplc.species_summary[,-1])))
  writeLines("==============================================\n Experiment #2 (transformation type) summary: \n==============================================")
  hits_uplc.species_summary
  
  table(hits_degrad.all$Species.x[hits_degrad.all$replicates==1])
  table(as.character(hits_degrad.all$GroupName)[hits_degrad.all$replicates==1])
  
  
  dim(all_hits %>% dplyr::full_join(hits_growth.all, by=c("drug.short"="source", "species.short"="target")))
  dim(all_hits %>% dplyr::full_join(hits_growth.all %>% dplyr::select(source, target), by=c("species.short"="source", "drug.short"="target")))
  
  setdiff(all_hits$species.short, hits_growth.all$target)
  setdiff(unique(hits_growth.all$target), all_hits$species.short)
  
  dim(all_hits)
  
  View(all_hits %>% dplyr::filter(is.na(exp1.is_hit)))
    
  readr::write_tsv(all_hits, "reports/all_data.tsv", col_names=T, na="")  
  
  #dplyr::left_join(drug_map %>% dplyr::select(drug.short, drug.long, drug.uplc_excluded)) %>%
  dplyr::mutate(f=is.na(drug.uplc_excluded)) %>%
  
  
  hits_degrad.species_summary = hits_degrad.all_counts %>%
    dplyr::group_by(species.long) %>% 
    dplyr::summarise(
      ndrugs.excluded=sum(nreplicates[!f]>0), ndrugs=sum(nreplicates[f]>0),
      nhits.excluded=sum(nhits[!f]>0), nhits=sum(nhits[f]>0),
      nreplicates.excluded=sum(nreplicates[!f]), nreplicates=sum(nreplicates[f])) %>%
    data.frame()
  
  hits_degrad.species_summary = rbind(hits_degrad.species_summary, c(species.long="Total", hits_degrad.all_counts %>%
    dplyr::summarise(
      ndrugs.excluded=sum(nreplicates[!f & !duplicated(drug.long)]>0), ndrugs=sum(nreplicates[f & !duplicated(drug.long)]>0),
      nhits.excluded=sum(nhits[!f]>0), nhits=sum(nhits[f]>0),
      nreplicates.excluded=sum(nreplicates[!f]), nreplicates=sum(nreplicates[f])) %>%
    data.frame()))
    
  
  hits_degrad.all_counts %>%
    dplyr::group_by(drug.long) %>% 
    dplyr::summarise(
      nspecies=paste0(sum(nreplicates[f]>0), " (", sum(nreplicates[!f]>0), ")"), 
      nhits=paste0(sum(nhits[f]>0), " (", sum(nhits[!f]>0), ")"), 
      nreplicates=paste0(sum(nreplicates[f]), " (", sum(nreplicates[!f]), ")"))
    
  pheatmap::pheatmap(d, color=c("#FFFFFF", RColorBrewer::brewer.pal(max(hits_degrad.all$replicates), "Reds")))
    
    )
}


