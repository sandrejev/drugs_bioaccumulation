library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(stringr)
library(httr)
library(xml2)
library(ggrepel)
library(fields)
library(VennDiagram)
library(gridExtra)
library(igraph)
library(KEGGREST)
source("functions.R")

exp8tpp.collect_data = function()
{
  tpp = read.table("data/exp8tpp/tpp.tsv", na.strings="", sep="\t", quote="", stringsAsFactors=F, header=T)
  
  concentrations = c(0, 0.5, 2, 10, 50)
  measurements.regex = "^(signal_sum|rel_fc_protein|norm_rel_fc_protein)_([0-9.]+)_?(unmodified|transformed)?"
  measurements.list = colnames(tpp)[grepl(measurements.regex, colnames(tpp))]
  tpp_long = tpp %>% reshape2::melt(measure.vars=measurements.list)
  vals = stringr::str_match(as.character(tpp_long$variable), measurements.regex); vals[is.na(vals)] = ""
  tpp_long$variable = paste0(vals[,2], ifelse(vals[,4]=="", "", paste0("_", vals[,4])))
  tpp_long$concentration = as.numeric(vals[,3])
  tpp_long = tpp_long %>% 
    dplyr::filter(!grepl("^#",gene_name) & !grepl("^#",protein_id) & qupm>=2) %>%
    dplyr::filter(variable=="norm_rel_fc_protein_unmodified")  %>%
    dplyr::rename(original_protein_id="protein_id") %>%
    dplyr::arrange(organism, preparation, variable, original_protein_id, concentration, temperature) %>%
    dplyr::distinct(organism, preparation, variable, original_protein_id, concentration, temperature, .keep_all=T) %>%
    dplyr::mutate(f=paste(organism, preparation, variable, original_protein_id, concentration, compound_effect)) %>%
    dplyr::mutate(f.lag1=dplyr::lag(f, 1, "NA"), f.lag2=dplyr::lag(f, 2, "NA"), f.lead1=dplyr::lead(f, 1, "NA"), f.lead2=dplyr::lead(f, 2, "NA")) %>%
    dplyr::mutate(f=!grepl("NA$", f) & (f == f.lag1 | f == f.lead1)) %>% #  | f == f.lead2 | f == f.lag2 %>%
    dplyr::mutate(compound_effect.f=ifelse(f, compound_effect, NA)) %>%
    dplyr::group_by(organism, preparation, variable, original_protein_id) %>%
    dplyr::mutate(i=concentration==concentration[1]) %>%
    dplyr::mutate(n_hits.f=sum(!is.na(compound_effect.f[i])), temp_hits.f=paste0(sort(temperature[i & !is.na(compound_effect.f)]), collapse=","),  has_hits.f=paste0(sort(unique(na.omit(compound_effect.f[i]))), collapse=","), has_hits.f=ifelse(has_hits.f=="destabilized,stabilized", "mixed", has_hits.f)) %>%
    dplyr::mutate(n_hits=sum(!is.na(compound_effect[i])), temp_hits=paste0(sort(temperature[i & !is.na(compound_effect)]), collapse=","),  has_hits=paste0(sort(unique(na.omit(compound_effect[i]))), collapse=","), has_hits=ifelse(has_hits=="destabilized,stabilized", "mixed", has_hits)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      compound_effect=tidyr::replace_na(compound_effect, "none"), 
      compound_effect.f=tidyr::replace_na(compound_effect.f, "none")) %>%
    dplyr::select(-i) %>%
    dplyr::select(preparation, organism, temperature, concentration, variable, gene_name, original_protein_id, compound_effect, value, n_hits, temp_hits, has_hits)
  write.table(tpp_long, "data/exp8tpp/tpp_long.tsv", na="", sep="\t", quote=F, row.names=F, col.names=T)
  
  
  # Mapping of single proteins to protein groups which can't be distinguished
  tppsum_long.sep = tpp_long %>% 
    dplyr::mutate(protein_id=original_protein_id) %>% 
    tidyr::separate_rows(protein_id, sep="\\|") %>%
    dplyr::distinct(organism, original_protein_id, protein_id) %>%
    dplyr::group_by(organism, protein_id) %>%
    dplyr::summarize(original_protein_id=strset(unlist(sapply(original_protein_id, strsplit, "\\|")))) %>%
    data.frame()
  write.table(tppsum_long.sep, "data/exp8tpp/tpp_protein_groups.tsv", na="", sep="\t", quote=F, row.names=F, col.names=T)
  
  
  #
  # Sum all the log2(fold-changes) per protein
  #
  tppsum_long = tpp_long %>%
    dplyr::group_by(organism, preparation, variable, original_protein_id, n_hits, temp_hits) %>%
    dplyr::summarise(
      gene_name=paste(unique(na.omit(gene_name)), collapse="|"),
      log2fc=sum(log2(value)), 
      has_hits=paste0(sort(unique(na.omit(has_hits))), collapse=","), has_hits=ifelse(has_hits!="" & !grepl("stabilized", has_hits), "mixed", has_hits)) %>%
    dplyr::group_by(organism, preparation, variable) %>%
    dplyr::mutate(log2fc_q=ecdf(log2fc)(log2fc)) %>%
    dplyr::mutate(log2fc_z=scale(log2fc)) %>%
    dplyr::ungroup()
  write.table(tppsum_long, "data/exp8tpp/tppsum_long.tsv", na="", sep="\t", quote=F, row.names=F, col.names=T)
  
  # create special (Wide) version of hits table and also split undestinguishable proteins to one protein per line
  tppsum_wide = tppsum_long %>%
    dplyr::filter(grepl("stabilized", has_hits)) %>%
    dplyr::inner_join(tppsum_long.sep, by=c("organism", "original_protein_id")) %>%
    reshape2::dcast(protein_id ~ preparation + organism, value.var="has_hits") %>%
    dplyr::rename(
      has_hits.lysate_Csaccharolyticum="lysate_C. saccharolyticum", 
      has_hits.lysate_IAI1="lysate_E. coli IAI1", 
      has_hits.whole_IAI1="whole_cell_E. coli IAI1", 
      has_hits.lysate_ED1a="lysate_E. coli ED1a", 
      has_hits.whole_ED1a="whole_cell_E. coli ED1a")
  write.table(tppsum_wide, "data/exp8tpp/tppsum_wide.tsv", na="", sep="\t", quote=F, row.names=F, col.names=T)
}

exp8tpp.analyze = function()
{
  # UNIPROT mappings
  uniprot2kegg = readr::read_tsv("data/db/uniprot2kegg.tsv", col_names=T, na="")
  gene2enzyme = readr::read_tsv("data/db/KEGG/kegg.gene2enzyme.tsv", col_names=T, na="") %>%
    dplyr::filter(!is.na(KEGG_UNIPROT))
  
  kegg.pathways = readr::read_delim("data/db/KEGG/kegg_pathways.tsv", "\t")
  kegg.gene2enzyme = readr::read_tsv("data/db/KEGG/kegg.gene2enzyme.tsv", col_names=T, na="") %>%
    dplyr::filter(!is.na(KEGG_UNIPROT))
  kegg.pathway2ec = readr::read_delim("data/db/KEGG/kegg.pathway2gene.tsv", "\t") %>%
    dplyr::distinct(KEGG_ORG, KEGG_PATHWAY, KEGG_EC) %>%
    dplyr::mutate(KEGG_EC=gsub("ec:", "", KEGG_EC))
  kegg.pathway2uniprot = readr::read_delim("data/db/KEGG/kegg.pathway2gene.tsv", "\t") %>%
    dplyr::filter(!is.na(KEGG_UNIPROT))
  
  # C. saccharolyticum/E. coli IAI1 blast
  csh2refseq = readr::read_delim("data/db/UNIPROT/csh_proteinids2refseq.list", "\t")
  ecr2refseq = readr::read_delim("data/db/UNIPROT/ecr_proteinids2refseq.list", "\t")
  ecq2refseq = readr::read_delim("data/db/UNIPROT/ecq_proteinids2refseq.list", "\t")
  csh2ecr = read.delim("data/db/UNIPROT/csh_ecr.m8", na="", header=F)
  colnames(csh2ecr) = c("qseqid", "sseqid", "qlen", "slen", "score", "bitscore", "length")
  csh2ecr = csh2ecr %>%
    dplyr::inner_join(csh2refseq %>% dplyr::rename(protein_id.Csaccharolyticum="From"), by=c("qseqid"="To")) %>%
    dplyr::inner_join(ecr2refseq %>% dplyr::rename(protein_id.IAI1="From"), by=c("sseqid"="To"))
  csh2ecr.f = csh2ecr %>% 
    dplyr::filter(bitscore>50 & length>0.8*qlen) %>% 
    dplyr::distinct(protein_id.Csaccharolyticum, protein_id.IAI1, bitscore)
  
  #
  # Read data files
  #
  tpp = readr::read_delim("data/exp8tpp/tpp.tsv", "\t")
  tpp_long = readr::read_delim("data/exp8tpp/tpp_long.tsv", "\t", col_types=c("temp_hits"="c", "has_hits"="c"))
  tppsum_long = readr::read_delim("data/exp8tpp/tppsum_long.tsv", "\t", col_types=c("temp_hits"="c", "has_hits"="c"))
  tppsum_long.sep = readr::read_delim("data/exp8tpp/tpp_protein_groups.tsv", "\t")
  tppsum_wide = readr::read_delim("data/exp8tpp/tppsum_wide.tsv", "\t", col_types=c("has_hits.whole_IAI1"="c", "has_hits.lysate_IAI1"="c", "has_hits.whole_ED1a"="c", "has_hits.lysate_ED1a"="c", "has_hits.lysate_Csaccharolyticum"="c"))
  
  #
  # Display correlation between E. coli strains sum log2 values
  #
  tppsum_long.fit = tppsum_long %>%
    dplyr::group_by(preparation, organism, variable) %>%
    dplyr::summarise(total_hits=sum(!is.na(has_hits)))    
  tppsum_long.ecoli = tppsum_long %>%
    dplyr::inner_join(tppsum_long, by=c("preparation", "variable", "gene_name")) %>%
    dplyr::filter(grepl("ED1a", organism.x) & grepl("IAI1", organism.y) & !is.na(gene_name)) %>% 
    setNames(gsub("\\.x$",".ED1a",names(.))) %>% 
    setNames(gsub("\\.y$",".IAI1",names(.)))
  
  tppsum_long.ecoli.cor = tppsum_long.ecoli %>%
    dplyr::group_by(preparation, variable) %>%
    dplyr::summarise(
      max=max(c(log2fc.ED1a, log2fc.IAI1), na.rm=T), 
      min=min(c(log2fc.ED1a, log2fc.IAI1), na.rm=T), 
      x=min+(max-min)*0.8, 
      y=min+(max-min)*0.2, 
      slope=round(lm(log2fc.ED1a~log2fc.IAI1)$coefficients[2], 2), R=round(cor(log2fc.IAI1, log2fc.ED1a, use="pairwise.complete.obs"), 2))
  
  pdf("reports/exp8tpp_ecoli.pdf", width=10, height=5)
  ggplot(tppsum_long.ecoli, aes(log2fc.IAI1, log2fc.ED1a)) +
    geom_abline(slope=1) +
    geom_point(alpha=0.2, shape=16) +
    geom_text(aes(x=x, y=y, label=paste0("slope: ", slope, "\nR: ", R)), size=6, data=tppsum_long.ecoli.cor) +
    geom_smooth(method="lm", fullrange=T, alpha=0.2) +
    facet_wrap(~preparation, scales="free") +
    labs(x="Sum fold-changes (IAI1)", y="Sum fold-changes (ED1a)")+
    theme_classic() +
    theme(aspect.ratio=1, axis.title=element_text(size=15), axis.text=element_text(size=13),
          legend.title=element_text(size=15), legend.text=element_text(size=13),
          strip.text = element_text(size=20),
          strip.background=element_blank(), 
          panel.grid=element_blank())
  
  
  
  venn.data = with(tppsum_long.ecoli, list(
    IAI1_whole=gene_name[!is.na(has_hits.IAI1) & preparation=="whole_cell"],
    IAI1_lysate=gene_name[!is.na(has_hits.IAI1) & preparation=="lysate"],
    ED1a_whole=gene_name[!is.na(has_hits.ED1a) & preparation=="whole_cell"],
    ED1a_lysate=gene_name[!is.na(has_hits.ED1a) & preparation=="lysate"]
  ))
  venn.colors = c(IAI1_whole="#FDC08699", IAI1_lysate="#386CB0FF", ED1a_whole="#FFD92F99", ED1a_lysate="#6DB5D2")[names(venn.data)]
  grid.newpage()
  grid.draw(VennDiagram::venn.diagram(
    venn.data, 
    fill=venn.colors,
    filename=NULL, force.unique=T))
  dev.off()
  
  tppenrichment.species = tppsum_long %>%
    dplyr::filter(grepl("C. saccharolyticum|IAI1", organism)) %>%
    dplyr::mutate(hit=ifelse(grepl("stabilized", has_hits), "Yes", "No")) %>%
    dplyr::select(organism, preparation, original_protein_id, has_hits)
  
  
  #
  # Analyze TPP hits and overlap
  #
  tppenrichment.species = tppsum_long %>%
    dplyr::filter(grepl("C. saccharolyticum|IAI1", organism)) %>%
    dplyr::mutate(hit=ifelse(grepl("stabilized", has_hits), "Yes", "No")) %>%
    dplyr::select(organism, preparation, original_protein_id, has_hits)
  
  #
  # Overlap between E. coli IAI1 and C. saccharolyticum WM1
  #
  # Cluster proteins based on BLAST
  g = igraph::graph_from_data_frame(csh2ecr.f, directed=F)
  gc = igraph::clusters(g)
  gc = data.frame(gc.names=names(gc$membership), gc.membership=gc$membership)
  
  # Use protein clusters of individual species to cluster cross-species proteins (significant blast hits between proteins of two species)
  tpp_csh2ecr.f = csh2ecr.f %>%
    dplyr::inner_join(gc, by=c("protein_id.Csaccharolyticum"="gc.names")) %>%
    dplyr::inner_join(gc, by=c("protein_id.IAI1"="gc.names")) %>%
    dplyr::mutate(cluster=ifelse(gc.membership.x==gc.membership.y, gc.membership.x, NA_integer_)) %>%
    dplyr::inner_join(tppsum_long.sep %>% dplyr::rename(original_protein_id.Csaccharolyticum="original_protein_id"), by=c("protein_id.Csaccharolyticum"="protein_id")) %>%
    dplyr::inner_join(tppsum_long.sep %>% dplyr::rename(original_protein_id.IAI1="original_protein_id"), by=c("protein_id.IAI1"="protein_id")) %>%
    dplyr::filter(!is.na(cluster)) %>%
    dplyr::distinct(original_protein_id.Csaccharolyticum, original_protein_id.IAI1, cluster) 
  
  # Overlap between E.coli IAI1 and C. saccharolyticum WM1 protein clusters in TPP
  tppenrichment.multispecies = tpp_csh2ecr.f %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(cluster_protein_id.Csaccharolyticum=strset(original_protein_id.Csaccharolyticum), cluster_protein_id.IAI1=strset(original_protein_id.IAI1)) %>%
    dplyr::mutate(protein_id.Csaccharolyticum=cluster_protein_id.Csaccharolyticum, protein_id.IAI1=cluster_protein_id.IAI1) %>%
    tidyr::separate_rows(protein_id.Csaccharolyticum, sep="\\|") %>%
    tidyr::separate_rows(protein_id.IAI1, sep="\\|") %>%
    dplyr::left_join(tppsum_wide %>% dplyr::filter(grepl("stabilized", has_hits.lysate_IAI1)) %>% dplyr::select(protein_id, has_hits.lysate_IAI1), by=c("protein_id.IAI1"="protein_id")) %>%
    dplyr::left_join(tppsum_wide %>% dplyr::filter(grepl("stabilized", has_hits.whole_IAI1)) %>% dplyr::select(protein_id, has_hits.whole_IAI1), by=c("protein_id.IAI1"="protein_id")) %>%
    dplyr::left_join(tppsum_wide %>% dplyr::filter(grepl("stabilized", has_hits.lysate_Csaccharolyticum)) %>% dplyr::select(protein_id, has_hits.lysate_Csaccharolyticum), by=c("protein_id.Csaccharolyticum"="protein_id")) %>%
    dplyr::group_by(cluster, cluster_protein_id.IAI1, cluster_protein_id.Csaccharolyticum) %>%
    dplyr::summarise(
      protein_id.lysate_Csaccharolyticum=strset(protein_id.Csaccharolyticum, grepl("stabilized", has_hits.lysate_Csaccharolyticum)),
      protein_id.lysate_IAI1=strset(protein_id.IAI1, grepl("stabilized", has_hits.lysate_IAI1)),
      protein_id.whole_IAI1=strset(protein_id.IAI1, grepl("stabilized", has_hits.whole_IAI1)),
      has_hits.lysate_Csaccharolyticum=strset(paste0(protein_id.Csaccharolyticum, ":", has_hits.lysate_Csaccharolyticum), grepl("stabilized", has_hits.lysate_Csaccharolyticum)),
      has_hits.lysate_IAI1=strset(paste0(protein_id.IAI1, ":", has_hits.lysate_IAI1), grepl("stabilized", has_hits.lysate_IAI1)),
      has_hits.whole_IAI1=strset(paste0(protein_id.IAI1, ":", has_hits.whole_IAI1), grepl("stabilized", has_hits.whole_IAI1))) %>%
    dplyr::ungroup()
  
  pdf("reports/exp8tpp_csh2ecoli_venn.pdf", width=10, height=5)
  venn.data = with(tppenrichment.multispecies, list(IAI1_whole=cluster[!is.na(has_hits.whole_IAI1)], CSH_lysate=cluster[!is.na(has_hits.lysate_Csaccharolyticum)], IAI1_lysate=cluster[!is.na(protein_id.lysate_IAI1)]))
  venn.colors = c(IAI1_whole="#FDC08699", CSH_lysate="#386CB0FF", IAI1_lysate="#FFD92F99")
  grid.newpage()
  grid.draw(VennDiagram::venn.diagram(
    venn.data, 
    fill=venn.colors,
    filename=NULL, force.unique=T))
  
  tppenrichment.multispecies.export = tppenrichment.multispecies %>% 
    dplyr::mutate(hit=ifelse(!is.na(has_hits.lysate_Csaccharolyticum) & !is.na(has_hits.whole_IAI1) & !is.na(protein_id.lysate_IAI1), "Yes", "No")) %>%
    dplyr::select(cluster_protein_id.Csaccharolyticum, cluster_protein_id.IAI1, hit)
  tppenrichment.species_export = tppenrichment.species %>%
    dplyr::mutate(hit=ifelse(grepl("stabilized", has_hits), "Yes", "No")) %>%
    dplyr::select(organism, preparation, original_protein_id, hit)
  readr::write_delim(tppenrichment.multispecies.export, "reports/exp8tpp_tppenrichment_multispecies.tsv", delim="\t", na="")
  readr::write_delim(tppenrichment.species_export, "reports/exp8tpp_tppenrichment_species.tsv", delim="\t", na="")
  dev.off()
  
  #
  # Export
  #
  tpp_export.a = tpp %>%
    dplyr::mutate(cluster_progein_id=protein_id) %>%
    tidyr::separate_rows(cluster_progein_id, sep="\\|")  
  tpp_export.b = tppenrichment.multispecies.export %>%
    dplyr::mutate(cluster_id=1:n()) %>%
    tidyr::separate_rows(cluster_protein_id.Csaccharolyticum, sep="\\|") %>%
    tidyr::separate_rows(cluster_protein_id.IAI1, sep="\\|") %>%
    dplyr::rename(`C. saccharolyticum`="cluster_protein_id.Csaccharolyticum", `E. coli IAI1`="cluster_protein_id.IAI1") %>%
    reshape2::melt(id.vars=c("cluster_id"), variable.name="organism", value.name="cluster_progein_id") %>%
    dplyr::distinct()
  tpp_export = tpp_export.a %>%
    dplyr::left_join(tpp_export.b, by=c("organism", "cluster_progein_id")) %>%
    dplyr::select(-cluster_progein_id) %>% 
    dplyr::group_by(preparation, organism, temperature, protein_id) %>%
    dplyr::mutate(
      gene_name=paste0(unique(unlist(strsplit(gene_name, "\\|"))), collapse="|"), 
      cluster_id=ifelse(all(is.na(cluster_id)), NA_integer_, unique(na.omit(cluster_id)))
    ) %>%
    dplyr::ungroup() %>%
    unique()
  
  tpp_rossman = readr::read_delim("data/exp8tpp/exp8tpp_rossman_data.tsv", "\t", col_types=c(rossman_family="c"))
  tpp_export.rossman = tpp_export %>%
    dplyr::left_join(tpp_rossman, by=c("organism", "preparation", "protein_id"="original_protein_id")) %>% 
    dplyr::left_join(kegg.gene2enzyme %>% dplyr::filter(!is.na(KEGG_UNIPROT)) %>% dplyr::select(KEGG_UNIPROT, KEGG_EC), by=c("protein_id"="KEGG_UNIPROT")) %>%
    dplyr::select(
      organism, preparation, temperature, protein_id, gene_name, 
      enzyme=KEGG_EC,
      compound_effect, 
      description, 
      qupm, cluster_id, rossman_hit, rossman_family, experiment, 
      signal_sum_0, signal_sum_0.5, signal_sum_2, signal_sum_10, signal_sum_50, 
      norm_rel_fc_protein_0.5_unmodified, norm_rel_fc_protein_2_unmodified, norm_rel_fc_protein_10_unmodified, norm_rel_fc_protein_50_unmodified)
  
  readr::write_delim(tpp_export.rossman, "reports/exp8tpp_export.tsv", delim="\t", na="")

  #
  # GO enrichment
  #
  tpp_enrichment.hits = tpp_long %>%
    dplyr::filter(concentration==50) %>%
    dplyr::group_by(organism, preparation, original_protein_id, variable) %>%
    dplyr::summarise(is_hit=sum(compound_effect==has_hits, na.rm=T)>0)
  tpp_enrichment.go = tppsum_long %>%
    dplyr::mutate(protein_id=original_protein_id) %>%
    tidyr::separate_rows(protein_id, sep="\\|") %>%
    dplyr::inner_join(readr::read_delim("data/db/uniprot2go.tsv", "\t"), by=c("protein_id"="uniprot_id")) %>%
    reshape2::dcast(organism+preparation+original_protein_id~go_term, value.var="go_term", fun.aggregate=length) %>%
    reshape2::melt(id.vars=c("organism", "preparation", "original_protein_id"), variable.name="go_term", value.name="in_go", factorsAsStrings=F) %>%
    dplyr::mutate(in_go=in_go>0)
  tpp_enrichment = tpp_enrichment.go %>%
    dplyr::inner_join(tpp_enrichment.hits, by=c("organism", "preparation", "original_protein_id")) %>%
    dplyr::group_by(organism, preparation, go_term) %>%
    dplyr::do((function(z){
      zz<<-z
      m11 = sum(z$is_hit & z$in_go)
      m12 = sum(z$is_hit & !z$in_go)
      m21 = sum(!z$is_hit & z$in_go)
      m22 = sum(!z$is_hit & !z$in_go)
      m = matrix(c(m11, m12, m21, m22), ncol=2, byrow=T)
      
      z.test = fisher.test(m)
      data.frame(odds=z.test$estimate, pvalue=z.test$p.value, hits_in_go=m11, hits_in_species=sum(z$is_hit), proteins_in_go=sum(z$in_go), proteins_in_species=nrow(z))
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(readr::read_delim("data/db/GO/go_terms.tsv", "\t"), by="go_term") %>%
    dplyr::mutate(is_significant=ifelse(pvalue<=0.1, "Yes", "No")) %>%
    dplyr::select(organism, preparation, go_term, go_description, is_significant, odds, pvalue, hits_in_go, hits_in_species, proteins_in_go, proteins_in_species)
  tpp_enrichment.f = tpp_enrichment %>%
    dplyr::filter(is_significant=="Yes") %>%
    dplyr::inner_join(readr::read_delim("data/db/GO/go_terms.tsv", "\t"), by="go_term") %>%
    dplyr::arrange(organism, preparation, dplyr::desc(abs(log2(odds))), pvalue) %>%
    dplyr::select(-is_significant)
  readr::write_delim(tpp_enrichment.f, "reports/exp8tpp_go_enrichment.tsv", delim="\t", na="")
  
  #
  # TPP / Duloxetine pull-down proteomics
  #
  pdf("reports/exp8tpp_tpp2click_venn.pdf", width=10, height=5)
  click.universe = readr::read_tsv("reports/exp3proteomics_volcano_plot_data.tsv", col_names=T, na="") %>%
    dplyr::mutate(hit.click=is_significant=="Yes") %>%
    dplyr::select(protein_id=proteinId, hit.click) %>% 
    dplyr::distinct(protein_id, hit.click)
  tpp.universe = tppsum_long %>% 
    dplyr::filter(grepl("sacch", organism) & preparation=="lysate") %>%
    dplyr::mutate(hit.tpp=!is.na(has_hits)) %>%
    tidyr::separate_rows(original_protein_id, sep="\\|") %>%
    dplyr::distinct(original_protein_id, hit.tpp)%>%
    dplyr::select(protein_id=original_protein_id, hit.tpp)
  tpp2click.export = click.universe %>%
    dplyr::inner_join(tpp.universe, by="protein_id")
  readr::write_delim(tpp2click.export, "reports/exp8tpp_tpp2click_overlap.tsv", delim="\t", na="")
  tpp2click.venn = with(tpp2click.export, list(tpp=protein_id[hit.click], click=protein_id[hit.tpp]))
  venn.colors = c(tpp="#386CB0FF", click="#38A0B0FF")
  grid.newpage()
  grid.draw(VennDiagram::venn.diagram(
    tpp2click.venn, 
    fill=venn.colors,
    filename=NULL, 
    force.unique=T))
  dev.off()
}

analyze.enrichment = function() {
  tppenrichment.multispecies.export = readr::read_delim("reports/exp8tpp_tppenrichment_multispecies.tsv", delim="\t")
  tppenrichment.species_export = readr::read_delim("reports/exp8tpp_tppenrichment_species.tsv", delim="\t")
  
  rossman_fold.IAI1 = data.frame(str=readLines("data/exp8tpp/enrichment/Ecoli/rosmann_matches.txt"), stringsAsFactors=F) %>%
    dplyr::mutate(organism="E. coli IAI1", rossman_hit="Yes") %>%
    tidyr::extract(str, c("some_id", "protein_id", "rossman_family"), regex='^( *\\d+ *)([^_]+)_.*"id":"([^"]+)".*')
  rossman_fold.CSH = data.frame(str=readLines("data/exp8tpp/enrichment/Csaccharolyticum/all_rosmann_matches.txt"), stringsAsFactors=F) %>%
    dplyr::mutate(organism="C. saccharolyticum", rossman_hit="Yes") %>%
    tidyr::extract(str, c("some_id", "protein_id", "rossman_family"), regex='^( *\\d+ *)([^_]+)_.*"id":"([^"]+)".*')
  rossman_fold = rbind(rossman_fold.IAI1, rossman_fold.CSH) %>%
    dplyr::distinct(organism, protein_id, .keep_all=T)
  
  tppenrichment = tppenrichment.multispecies.export %>%
    dplyr::mutate(original_protein_id=paste0(cluster_protein_id.Csaccharolyticum, "|", cluster_protein_id.IAI1), preparation="lysate", organism="multispecies") %>%
    select(organism, preparation, original_protein_id, hit)
  tppenrichment = rbind(tppenrichment, tppenrichment.species_export)
  tppenrichment.rossman = tppenrichment %>%
    dplyr::mutate(protein_id=original_protein_id) %>%
    tidyr::separate_rows(protein_id, sep="\\|") %>%
    dplyr::left_join(rossman_fold, by=c("organism", "protein_id")) %>%
    dplyr::mutate(rossman_hit=tidyr::replace_na(rossman_hit, "No")) %>%
    dplyr::group_by(organism, preparation, original_protein_id, hit) %>%
    dplyr::summarise(
      rossman_hit=ifelse(any(rossman_hit=="Yes"), "Yes", "No"), 
      rossman_family=paste(unique(rossman_family), collapse="|"), 
      rossman_family=ifelse(rossman_family=="", NA_character_, rossman_family))
  readr::write_delim(tppenrichment.rossman, "reports/exp8tpp_rossman_data.tsv", delim="\t", na="")
  
  
  tppenrichment.rossman_sum = tppenrichment.rossman %>%
    dplyr::group_by(organism, preparation) %>%
    dplyr::do((function(z) {
      m11 = sum(z$hit=="Yes" & z$rossman_hit=="Yes")
      m12 = sum(z$hit=="Yes" & z$rossman_hit=="No")
      m21 = sum(z$hit=="No" & z$rossman_hit=="Yes")
      m22 = sum(z$hit=="No" & z$rossman_hit=="No")
      t = fisher.test(matrix(c(m11,m12,m21,m22), ncol=2, byrow=T))
      
      data.frame(pvalue=t$p.value, odds=t$estimate, hit1_rossman1=m11, hit1_rossman0=m12, hit0_rossman1=m21, hit0_rossman0=m22)
    })(.))
  
  readr::write_delim(tppenrichment.rossman_sum, "reports/exp8tpp_rossman_enrichment.tsv", delim="\t", na="")
  
  # organism           preparation        pvalue    odds   hit1_rossman1 hit1_rossman0 hit0_rossman1 hit0_rossman0
  # C. saccharolyticum lysate              0.631    1.10              42           114           469          1404
  # E. coli IAI1       lysate        0.000000143   0.470              61           309           493          1173
  # E. coli IAI1       whole_cell          0.236   0.844              92           271           386           960
  # multispecies       lysate                  1       0               0             6             0           519
}