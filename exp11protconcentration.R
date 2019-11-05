library(ggplot2)
library(readxl)
library(dplyr)
library(VennDiagram)

exp10bioaccumulation.analyze = function() {
  # top3:      (sum intensity) of the 3 highest intense peptides (based on 3 best flyers) – problematic, as often less than 3 peptides are observed
  # mw:        Molecular weight
  # qssm:      quantified spectrum-sequence matches, which
  #            are the number of spectrum to sequence matches [peptides] in a protein group with
  #            reporter ions used in protein quantification.
  # ssm:       quantified spectrum-sequence matches
  # qupm:      quantified unique peptide matches, which describe
  #            the number of unique peptide sequences in a protein group with reporter ions used
  #            in protein quantification
  # max_score: Maximum Mascot score of all peptides in protein group.
  #
  # https://media.nature.com/original/nature-assets/nprot/journal/v10/n10/extref/nprot.2015.101-S1.pdf
  #
  data_wide = do.call(dplyr::bind_rows, lapply(list.files(path="data/exp11protconcentration/", pattern="*.xlsx$", full.names=T), function(f) {
    readxl::read_xlsx(f, sheet=1) %>% dplyr::mutate(organism=gsub("_.*", "", basename(f)))
  })) %>%
    dplyr::filter(qupm>=2)  %>% 
    dplyr::mutate(
      mw=as.numeric(mw),
      organism=dplyr::case_when(
        grepl("Csaccha", organism)~"C. saccharolyticum", 
        grepl("ED1a", organism)~"E. coli ED1a", 
        grepl("IAI1", organism)~"E. coli IAI1", 
        T~"Other"))
  

  #
  # Calculate long version
  #
  channel_mapping = c("126"="0uM", "127L"="05uM", "127H"="2uM", "128L"="10uM", "128H"="50uM", "129L"="0uM", "129H"="05uM", "130L"="2uM", "130H"="10uM", "131L"="50uM")
  data_long = data_wide %>%
    dplyr::filter(qupm>=2) %>%
    reshape2::melt(measure.vars=colnames(data_wide)[grepl("signal_sum", colnames(data_wide))], variable.name="channel", value.name="signal") %>%
    dplyr::mutate(channel=gsub("signal_sum_", "", channel)) %>%
    dplyr::mutate(concentration=channel_mapping[channel], replicate=ifelse(as.numeric(gsub("[A-Z]+", "", channel))<129, "rep1", "rep2")) %>%
    dplyr::mutate(concentration_r=dplyr::case_when(concentration=="0uM"~0, concentration=="05uM"~0.5, concentration=="2uM"~2, concentration=="10uM"~10, concentration=="50uM"~50)) %>%
    reshape2::melt(measure.vars=colnames(data_wide)[grepl("log2_fc", colnames(data_wide))], variable.name="fc_concentration", value.name="log2_fc.bootstraped") %>%
    dplyr::mutate(fc_replicate=gsub(".*_(rep[^_])", "\\1", fc_concentration), fc_concentration=gsub(".*fc_([^_]+).*", "\\1", fc_concentration)) %>%
    dplyr::filter(concentration==fc_concentration & replicate==fc_replicate) %>%
    dplyr::mutate(signal=ifelse(signal<0, NA_real_, signal)) %>%
    dplyr::mutate(signal=ifelse(grepl("sacch", organism) & concentration=="50uM", NA_real_, signal)) %>%
    dplyr::mutate(log2_fc.bootstraped=ifelse(grepl("sacch", organism) & concentration=="50uM", NA_real_, log2_fc.bootstraped))

  
  #
  # Calculate relative information
  #
  data_long = data_long %>%
    dplyr::inner_join(data_long %>% dplyr::filter(concentration=="0uM") %>% dplyr::select(organism, protein_id, replicate, signal.0=signal), by=c("organism", "protein_id", "replicate")) %>%
    dplyr::group_by(organism, protein_id, concentration, concentration_r) %>%
    dplyr::mutate(
      n=n(),
      pvalue=ifelse(any(is.na(signal) | is.na(signal.0)), NA_real_, t.test(signal, signal.0)$p.value), 
      log2_fc.tmp = mean(signal, na.rm=T)-mean(signal.0, na.rm=T),
      log2_fc=ifelse(all(concentration=="0uM") | all(is.na(signal-signal.0)), NA_real_, log2_fc.tmp),
      signal=signal, signal.0=mean(signal.0, na.rm=T),
      log2_fc.bootstraped_mean=mean(log2_fc.bootstraped, na.rm=T),
      is_significant=ifelse(!is.na(pvalue) & pvalue>0.05, "No","Yes"),
      is_hit=ifelse(is_significant=="Yes" & !is.na(log2_fc.bootstraped) & abs(log2_fc.bootstraped)>=1, "Yes", "No")) %>%
    dplyr::group_by(organism, protein_id) %>%
    dplyr::mutate(is_consistent=ifelse(all(0.5>=log2_fc.bootstraped_mean[concentration!="0uM"], na.rm=T) | all(-0.5<=log2_fc.bootstraped_mean[concentration!="0uM"], na.rm=T), "Yes", "No")) %>%
    data.frame()
  
  ggplot(data_long) +
    geom_histogram(aes(pvalue, fill=concentration)) +
    facet_wrap(concentration~organism)
 
  #
  # Export
  #
  data_export = data_long %>%
    dplyr::mutate(signal=round(signal,2)) %>%
    dplyr::group_by(organism, protein_id) %>%
    dplyr::arrange(organism, protein_id, concentration_r) %>%
    dplyr::mutate(
      is_significant=ifelse(any(is_significant[max(which(!is.na(signal)))]=="Yes"), "Yes", "No"), 
      is_hit=ifelse(any(is_hit[max(which(!is.na(signal)))]=="Yes"), "Yes", "No")) %>%
    dplyr::mutate(
      log2_fc_05uM=mean(log2_fc.bootstraped[concentration=="05uM"], na.rm=T),
      log2_fc_2uM=mean(log2_fc.bootstraped[concentration=="2uM"], na.rm=T),
      log2_fc_10uM=mean(log2_fc.bootstraped[concentration=="10uM"], na.rm=T),
      log2_fc_50uM=mean(log2_fc.bootstraped[concentration=="50uM"], na.rm=T)) %>%
    reshape2::dcast(organism+protein_id+gene_name+protein_fdr+qupm+is_significant+is_hit+is_consistent+log2_fc_05uM+log2_fc_2uM+log2_fc_10uM+log2_fc_50uM~concentration+replicate, value.var="signal") %>%
    dplyr::mutate_all(tidyr::replace_na, replace=NA)
  readr::write_delim(data_export, path="reports/exp11protconcentration_data.tsv", delim="\t", na="", col_names=T)
  
  #
  # Cummulative hits plot
  #
  data_long.sum = data_long %>%
    dplyr::filter(concentration!="0uM" & !is.na(log2_fc.bootstraped)) %>%
    dplyr::group_by(organism, protein_id, gene_name, concentration, pvalue, log2_fc, signal, signal.0, is_significant, is_consistent, is_hit) %>%
    dplyr::summarise(log2_fc.bootstraped=mean(log2_fc.bootstraped, na.rm=T))
  data_long.cumsum = data_long.sum %>%
    tidyr::crossing(data.frame(log2_fc.threshold=seq(0.5, 4, 0.1))) %>%
    dplyr::group_by(organism, concentration, log2_fc.threshold) %>%
    dplyr::summarise(count=sum(abs(log2_fc.bootstraped)>=log2_fc.threshold, na.rm=T))
  
  
  pdf("reports/exp11protconcentration_hits.pdf", width=8, height=5)
  data_long.cumsum_other = data_long.cumsum %>% dplyr::filter(organism!="C. saccharolyticum")
  ggplot(mapping=aes(x=log2_fc.threshold, y=count, color=concentration), data=data_long.cumsum_other) +
    geom_line() +
    ggrepel::geom_text_repel(aes(label=count), data=data_long.cumsum_other %>% dplyr::filter(log2_fc.threshold==1)) +
    facet_wrap(~organism)+ 
    scale_color_manual(values=c("05uM"="#999999", "2uM"="#B27F72", "10uM"="#FF3200", "50uM"="#891B00")) +  
    labs(y="Number of hits", x="log2(FC) threshold") +
    theme_classic(base_size=18) +
    theme(
      strip.text.x = element_text(size=24),
      strip.background = element_blank(),
      legend.justification=c(1,1), 
      legend.position=c(0.95, 0.95), 
      aspect.ratio=1)
  
  data_long.cumsum_csacch = data_long.cumsum %>% dplyr::filter(organism=="C. saccharolyticum")
  ggplot(mapping=aes(x=log2_fc.threshold, y=count, color=concentration), data=data_long.cumsum_csacch) +
    geom_line() +
    ggrepel::geom_text_repel(aes(label=count), data=data_long.cumsum_csacch %>% dplyr::filter(log2_fc.threshold==1)) +
    facet_wrap(~organism)+ 
    scale_color_manual(values=c("05uM"="#999999", "2uM"="#B27F72", "10uM"="#FF3200", "50uM"="#891B00")) +  
    labs(y="Number of hits", x="log2(FC) threshold") +
    theme_classic(base_size=18) +
    theme(
      strip.text.x = element_text(size=24),
      strip.background = element_blank(),
      legend.justification=c(1,1), 
      legend.position=c(0.95, 0.95), 
      aspect.ratio=1)
  dev.off()
  
  
  #
  # hits overlap between concentrations
  #
  grid::grid.newpage()
  pushViewport(plotViewport(layout=grid.layout(1, 3),gp=gpar(cex=2)))
  organisms = unique(data_long.sum$organism)
  for(org in organisms) {
    data_long.sum.f = data_long.sum %>% 
      dplyr::filter(organism==org & concentration!="0uM" & is_hit=="Yes" & is_consistent=="Yes") %>%
      dplyr::group_by(protein_id) %>%
      dplyr::mutate(n_hits=n())
    venn.data = split(data_long.sum.f$protein_id, data_long.sum.f$concentration)
    venn.colors = c("05uM"="#E5C2B9", "2uM"="#B27F72", "10uM"="#FF3200", "50uM"="#891B00")[names(venn.data)]
    venn = VennDiagram::venn.diagram(
      venn.data,
      main=org,
      fill=venn.colors,
      filename=NULL, force.unique=T)
    
    pushViewport(plotViewport(layout.pos.col=which(org==organisms)))
    grid.draw(venn)
    popViewport()
  }
  
  #
  # GO term enrichment
  #
  lu = function(x,f) { length(unique(x$protein_id[f])) }
  go_descriptions = readr::read_delim("data/db/GO/go_terms.tsv", "\t")
  go_terms = readr::read_delim("data/db/uniprot2go.tsv", "\t")

  
  data_long.go = data_long %>% 
    dplyr::mutate(original_protein_id=protein_id) %>%
    tidyr::separate_rows(protein_id, sep="\\|") %>%
    dplyr::inner_join(go_terms, by=c("protein_id"="uniprot_id")) %>%
    reshape2::dcast(organism+concentration+original_protein_id~go_term, value.var="go_term", fun.aggregate=length) %>%
    reshape2::melt(id.vars=c("organism", "concentration", "original_protein_id"), variable.name="go_term", value.name="in_go", factorsAsStrings=F) %>%
    dplyr::mutate(in_go=in_go>0)
  data_long.hits = data_long %>%
    dplyr::mutate(original_protein_id=protein_id) %>%
    dplyr::group_by(organism, concentration, original_protein_id) %>%
    dplyr::summarise(is_hit=sum(is_hit=="Yes", na.rm=T)>0)
  data_enrichment = data_long.go %>%
    dplyr::inner_join(data_long.hits, by=c("organism", "concentration", "original_protein_id")) %>%
    dplyr::group_by(organism, concentration, go_term) %>%
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
    dplyr::select(organism, concentration, go_term, go_description, is_significant, odds, pvalue, hits_in_go, hits_in_species, proteins_in_go, proteins_in_species)
  data_enrichment.f = data_enrichment %>%
    dplyr::filter(is_significant=="Yes") %>%
    dplyr::arrange(organism, concentration, dplyr::desc(abs(log2(odds))), pvalue) %>%
    dplyr::select(-is_significant)
  readr::write_delim(data_long, "reports/exp11protconcentration_go_enrichment.tsv", delim="\t", na="")
}
