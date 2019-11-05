library(dplyr)
library(readr)
library(ggplot2)
library(data.table)
library(FactoMineR)
library(factoextra)
library(VennDiagram)
library(randomForest)
library(KEGGREST)
library(foreach)
library(VennDiagram)
source("functions.R")
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


exp9metconcentration.analyze = function()
{
  #
  # Load KEGG links
  #
  uniprot2kegg = readr::read_tsv("data/db/uniprot2kegg.tsv", col_names=T, na="")
  kegg.gene2cpd = readr::read_delim("data/db/KEGG/species2cpd.tsv", "\t") %>% 
    dplyr::mutate(KEGG_GENE=gsub(".*:", "", KEGG_GENE), KEGG_CPD=gsub("cpd:", "", KEGG_CPD))
  kegg.species2cpd = kegg.gene2cpd %>%
    reshape2::dcast(KEGG_CPD ~ Species, value.var="KEGG_CPD", fun.aggregate=length) %>%
    dplyr::mutate_at(dplyr::vars(-dplyr::matches("KEGG")), list(~.>0))
  kegg.pathways = readr::read_delim("data/db/KEGG/kegg_pathways.tsv", "\t")
  kegg.species2org = readr::read_delim("data/db/KEGG/kegg_species2org.tsv", "\t")
  kegg.pathway2cpd = readr::read_delim("data/db/KEGG/kegg_pathway2cpd.tsv", delim="\t") %>%
    dplyr::inner_join(kegg.species2org, by="KEGG_ORG") %>%
    dplyr::mutate(KEGG_CPD=gsub("cpd:", "", KEGG_CPD)) %>%
    dplyr::distinct(KEGG_ORG, Species_short, KEGG_PATHWAY, KEGG_CPD)
  kegg.uniprot2cpd = kegg.gene2cpd %>% dplyr::inner_join(uniprot2kegg, by=c("KEGG_GENE"="kegg_protein_id")) 
  hmdb2kegg = readr::read_delim("data/db/KEGG/hmdb2kegg.list", "\t")
  
  
  ppm_extended = 100e-6
  ppm = 10e-6
  
  #
  # Calculate expected m/z peaks for roflumilast and duloxetine considering possible adducts and isotopes
  #
  electron_mass = 0.00054858026 
  drug_formulas = readr::read_delim("data/exp9metconcentration/drug_formulas.tsv", "\t", comment="#")
  adducts = readr::read_delim("data/db/adducts.tsv", "\t", comment="#") %>%
    tidyr::separate(adduct_formula, into=c("adduct_formula_add", "adduct_formula_sub"), sep="-")
  drug_peaks = drug_formulas %>%
    tidyr::crossing(adducts) %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z) {
      z.molecule = Rdisop::getMolecule(z$formula)
      if(z$adduct_nmol > 1) {
        for(i in 2:z$adduct_nmol) { z.molecule = Rdisop::addMolecules(z.molecule$formula, z$formula) }
      }
      
      if(!is.na(z$adduct_formula_add) && z$adduct_formula_add!="") {
        z.molecule = Rdisop::addMolecules(z.molecule$formula, z$adduct_formula_add)
      }
      if(!is.na(z$adduct_formula_sub) && z$adduct_formula_sub!="") {
        z.molecule = Rdisop::subMolecules(z.molecule$formula, z$adduct_formula_sub)
      }
      ion_mz = (z.molecule$isotopes[[1]][1,] - z$adduct_charge*electron_mass)/abs(z$adduct_charge)
      cbind(data.frame(z), abundance=z.molecule$isotopes[[1]][2,], mz=ion_mz)
    })(.)) %>%
    dplyr::mutate(lb=mz*(1-ppm), ub=mz*(1+ppm), lb_extended=mz*(1-ppm_extended), ub_extended=mz*(1+ppm_extended))
  
  #
  # Read metabolomics data
  #
  metabolomics.anno = readr::read_delim("data/exp9metconcentration/METABOLOMICS_DATA_ANNOTATION.csv", ",") %>%
    dplyr::mutate(KEGG_ID=gsub("cpd:| ","", KEGG_ID))
  metabolomics.anno.f = metabolomics.anno %>%
    dplyr::filter(Assumed_neutral_mass_shift=="none" & Assumed_ESI_adduct %in% c("[M-2H]", "[M-H]") & Fraction_ion_detected > 80) %>%
    tidyr::separate_rows(HMDB_legacy_ID, sep="\\|") %>% 
    dplyr::left_join(hmdb2kegg, by=c("HMDB_legacy_ID"="HMDB_ID")) %>%
    dplyr::group_by(Ion_mz, Metabolite_name_HMDB, Assumed_ESI_adduct, Assumed_neutral_mass_shift) %>%
    dplyr::arrange(!is.na(cpd)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(KEGG_ID=ifelse(is.na(KEGG_ID), cpd, KEGG_ID)) %>% 
    dplyr::group_by(Ion_mz, Chemical_formula, Assumed_ESI_adduct, Assumed_neutral_mass_shift) %>%
    dplyr::filter(all(is.na(KEGG_ID)) | !is.na(KEGG_ID))
  
  metabolomics.ions = readr::read_delim("data/exp9metconcentration/METABOLOMICS_DATA_IONS.csv", ",")
  metabolomics.meta = readr::read_delim("data/exp9metconcentration/METABOLOMICS_DATA_METADATA.csv", ",", quote='"') %>%
    dplyr::mutate(SPECIES=dplyr::case_when(
      grepl("none", SPECIES) ~ NA_character_,
      grepl("sacch", SPECIES) ~ "C. saccharolyticum",
      grepl("ED1a", SPECIES) ~ "E. coli ED1a",
      grepl("IAI1", SPECIES) ~ "E. coli IAI1",
      grepl("WCFS1", SPECIES) ~ "L. plantarum",
      grepl("IL1403", SPECIES) ~ "L. lactis",
      grepl("salivarius", SPECIES) ~ "S. salivarius")) %>%
    dplyr::mutate(INTERACTION=dplyr::case_when(
      grepl("Duloxetine", DIFFERENTIAL_TREATMENT_OR_CONDITION) & grepl("IAI|sacc|plant|salivarius", SPECIES) ~ "Yes",
      grepl("Roflumilast", DIFFERENTIAL_TREATMENT_OR_CONDITION) & grepl("lactis", SPECIES) ~ "Yes",
      T ~ "No"))
  metabolomics.data = t(readr::read_delim("data/exp9metconcentration/METABOLOMICS_DATA_DATA.csv", ",")[,-1])
  colnames(metabolomics.data) = metabolomics.ions$Ion_mz
  metabolomics_long.data = as.data.frame(metabolomics.data)
  colnames(metabolomics_long.data) = metabolomics.ions$Ion_index
  metabolomics_long.data = metabolomics_long.data %>%
    dplyr::mutate(SAMPLE_INDEX=metabolomics.meta$SAMPLE_INDEX) %>%
    reshape2::melt(variable.name="Ion_index", value.name="Intensity", id.vars="SAMPLE_INDEX") %>%
    dplyr::inner_join(metabolomics.meta, by="SAMPLE_INDEX")
  

  #
  # Select ions
  #
  # Remove all ions that have even SMALL CHANCE of being related tested drugs
  metabolomics.ions.no_drugs2 = metabolomics.ions %>%
    dplyr::filter(!data.table::inrange(Ion_mz, drug_peaks$lb_extended, drug_peaks$ub_extended))
  metabolomics.data.no_drugs2 = metabolomics.data[,metabolomics.ions$Ion_mz %in% metabolomics.ions.no_drugs2$Ion_mz]

  #
  # Calculate T-test 0-concentration vs X-concentration
  #
  th.pval = 0.05
  th.diff = 0.1
  
  #
  # Test one ion difference to zero concentration samples at a time
  #
  metabolomics_long.data.0 = metabolomics_long.data %>% 
    dplyr::filter(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY==0) %>%
    dplyr::select(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, BIOLOGICAL_REPLICATE, Ion_index, Intensity)
  metabolomics.ttest = metabolomics_long.data %>% 
    dplyr::filter(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY>0) %>%
    dplyr::inner_join(metabolomics_long.data.0 ,c("SPECIES", "DIFFERENTIAL_TREATMENT_OR_CONDITION", "BIOLOGICAL_REPLICATE", "Ion_index")) %>%
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY, Ion_index) %>%
    dplyr::do((function(z) {
      zz<<-z
      if(all(z$Intensity.x==z$Intensity.y)) { return(data.frame(p.value=1, fc=1)) }
      z.fc = mean(z$Intensity.x, na.rm=T)/mean(z$Intensity.y, na.rm=T)
      z.rel = (mean(z$Intensity.x, na.rm=T)-mean(z$Intensity.y, na.rm=T))/mean(z$Intensity.y, na.rm=T)
      z.test = t.test(z$Intensity.x, z$Intensity.y)

      data.frame(p.value=z.test$p.value, fc=z.fc, rel=z.rel)
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(fc=ifelse(is.finite(fc), fc, 10), rel=ifelse(is.finite(rel), rel, 1), Ion_index=as.numeric(Ion_index)) 
  # ggplot(metabolomics.ttest) +
  #   geom_histogram(aes(x=p.value, fill=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY)), alpha=0.3, position="identity") +
  #   facet_wrap(~DIFFERENTIAL_TREATMENT_OR_CONDITION) +
  #   labs(fill="Concentration")
  
  
  #
  # Calculate hits (Must only consistenly change into one direction or stay the same)
  #
  metabolomics.hits = metabolomics.ttest %>%
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, Ion_index) %>%
    dplyr::mutate(
      rel.max=max(rel[p.value<=th.pval]), 
      rel.min=min(rel[p.value<=th.pval]), 
      increase=sum(rel[p.value<=th.pval]>=th.diff), 
      decrease=sum(rel[p.value<=th.pval]<=-th.diff),
      is_hit=ifelse(xor(increase>=1, decrease>=1), "Yes", "No")
      ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(pval=p.value, log2_fc=log2(fc)) %>%
    reshape2::melt(measure.vars=c("pval", "log2_fc")) %>%
    reshape2::dcast(SPECIES+DIFFERENTIAL_TREATMENT_OR_CONDITION+Ion_index+rel.max+rel.min+increase+decrease+is_hit ~ variable+DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY, value.var="value")
  
  #
  # Select only duloxetine/C. saccharolyticum hits
  #
  metabolomics.csh_hits = metabolomics.hits %>%
    dplyr::filter(grepl("sacch", SPECIES) & DIFFERENTIAL_TREATMENT_OR_CONDITION=="Duloxetine" & is_hit=="Yes")
  
  # 
  # Annotate duloxetine/C. saccharolyticum hits
  #
  metabolomics.csh_annotated_hits = metabolomics.csh_hits %>%
    dplyr::inner_join(metabolomics.anno.f, by="Ion_index") %>%
    dplyr::filter(!is.na(KEGG_ID)) %>%
    dplyr::mutate(Ion_index=as.numeric(as.character(Ion_index))) %>%
    dplyr::mutate(n_formula=length(unique(Chemical_formula))) %>%
    dplyr::mutate(Display_color=rainbow(n_formula)[match(Chemical_formula, unique(Chemical_formula))]) %>%
    dplyr::mutate(Display_weight=pmax(abs(rel.max), abs(rel.min))) %>%
    dplyr::distinct(KEGG_ID, .keep_all=T) %>%
    data.frame()
  metabolomics.kegg = with(metabolomics.csh_annotated_hits, paste0(KEGG_ID, " ", Display_color, " W", round(Display_weight*200, 2)))
  writeLines(metabolomics.kegg)

  #
  # Export data
  #
  export.kegg = metabolomics.ions %>% 
    dplyr::select(Ion_index, Ion_mz, First_KEGG_ID) %>%
    dplyr::left_join(metabolomics.anno.f %>% dplyr::group_by(Ion_index) %>% dplyr::summarise(KEGG_ID=paste0(na.omit(unique(KEGG_ID)), collapse=";")), by="Ion_index") %>%
    dplyr::mutate(KEGG_ID=ifelse(is.na(KEGG_ID) | KEGG_ID=="", First_KEGG_ID, KEGG_ID)) %>%
    dplyr::mutate(KEGG_ID=gsub("cpd:","", KEGG_ID)) %>%
    dplyr::select(Ion_index, Ion_mz, KEGG_ID)
  export.anno = metabolomics.anno %>%
    dplyr::select(Ion_mz, Metabolite_name=Metabolite_name_HMDB, Chemical_formula, Assumed_ESI_adduct, Assumed_neutral_mass_shift, Number_detected_isotopologues, Isotopologue_pattern_MRE, CV_ion_intensity, Fraction_ion_detected, HMDB_ID, KEGG_ID, Pubchem_ID=Pubchem_compound_ID, CAS_Nr=CAS_Nr_ID, Chebi_ID, FoodDB_ID, Chemspider_ID)
  export.hits = metabolomics.hits
  export.data = metabolomics_long.data %>%
    dplyr::mutate(Media=dplyr::case_when(
      grepl("ACN", GENERAL_TREATMENT_OR_CONDITION) ~ "ACN",
      grepl("GMM", TISSUE_OR_BODY_FLUID_OR_CELL_LINE) ~ "GMM",
      T~"ACN")) %>%
    dplyr::mutate(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY=tidyr::replace_na(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY, 0.0)) %>%
    dplyr::mutate(Ion_index=as.numeric(as.character(Ion_index))) %>%
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, Media, INTERACTION, Ion_index) %>%
    dplyr::mutate(InPCA_DrugSpecies=ifelse(mean(Intensity>0)==1, "Yes", "No")) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(metabolomics.ions.no_drugs2 %>% dplyr::mutate(InPCA_Drug="Yes") %>% dplyr::select(Ion_index, InPCA_Drug), by="Ion_index") %>%
    dplyr::mutate(InPCA_Drug=tidyr::replace_na(InPCA_Drug, "No")) %>%
    dplyr::inner_join(export.kegg, by="Ion_index") %>% 
    dplyr::mutate(str1="Intensity", str2="rep") %>%
    reshape2::dcast(SPECIES+Media+DIFFERENTIAL_TREATMENT_OR_CONDITION+InPCA_DrugSpecies+InPCA_Drug+KEGG_ID+Ion_mz+Ion_index ~ str1+DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY+str2+BIOLOGICAL_REPLICATE, value.var="Intensity") %>%
    dplyr::left_join(export.hits, by=c("SPECIES", "DIFFERENTIAL_TREATMENT_OR_CONDITION", "Ion_index")) %>%
    data.frame() %>%
    dplyr::rename(Species="SPECIES", Drug="DIFFERENTIAL_TREATMENT_OR_CONDITION") %>%
    dplyr::arrange(is.na(Species), Species, is.na(Drug), Drug, is.na(Media), Media, Ion_mz) %>%
    dplyr::select(
      Species, Media, Drug, 
      InPCA_DrugSpecies, InPCA_Drug, 
      KEGG_ID, Ion_mz, is_hit,
      Intensity_0_rep_1, Intensity_0_rep_2, Intensity_0_rep_3, Intensity_0_rep_4, 
      Intensity_12.5_rep_1, Intensity_12.5_rep_2, Intensity_12.5_rep_3, Intensity_12.5_rep_4,
      Intensity_25_rep_1, Intensity_25_rep_2, Intensity_25_rep_3, Intensity_25_rep_4, 
      Intensity_50_rep_1, Intensity_50_rep_2, Intensity_50_rep_3, Intensity_50_rep_4,
      pval_12.5, pval_25, pval_50)
  readr::write_delim(export.data %>% dplyr::select(-is_hit), "reports/exp9metconcentration_data.tsv", delim="\t", na="")
  readr::write_delim(export.anno, "reports/exp9metconcentration_anno.tsv", delim="\t", na="")
  
  
  #
  # Calculate KEGG pathway enrichment
  #
  export_cpd.data = export.data %>%
    dplyr::filter((rowSums(.[grepl("Intensity_", colnames(.))]>1))>0) %>% 
    dplyr::filter(!is.na(KEGG_ID)) %>%
    tidyr::separate_rows(KEGG_ID, sep=";") %>%
    dplyr::inner_join(kegg.pathway2cpd, by=c("Species"="Species_short", "KEGG_ID"="KEGG_CPD"))
  export.enrichment_pathway = export_cpd.data %>%
    reshape2::dcast(Species+Drug+Ion_mz~KEGG_PATHWAY, value.var="KEGG_PATHWAY", fun.aggregate=length) %>%
    reshape2::melt(id.vars=c("Species", "Drug", "Ion_mz"), variable.name="KEGG_PATHWAY", value.name="in_pathway", factorsAsStrings=F) %>%
    dplyr::mutate(in_pathway=in_pathway>0)
  export.enrichment_hits = export_cpd.data %>%
    dplyr::group_by(Species, Drug, Ion_mz) %>%
    dplyr::summarise(is_hit=any(is_hit=="Yes"))
  export.enrichment = export.enrichment_pathway %>%
    dplyr::inner_join(export.enrichment_hits, by=c("Species", "Drug", "Ion_mz")) %>%
    dplyr::group_by(Species, Drug, KEGG_PATHWAY) %>%
    dplyr::do((function(z){
      zz<<-z
      m11 = sum(z$is_hit & z$in_pathway)
      m12 = sum(z$is_hit & !z$in_pathway)
      m21 = sum(!z$is_hit & z$in_pathway)
      m22 = sum(!z$is_hit & !z$in_pathway)
      m = matrix(c(m11, m12, m21, m22), ncol=2, byrow=T)
      
      z.test = fisher.test(m)
      data.frame(odds=z.test$estimate, pvalue=z.test$p.value, hits_in_pathway=m11, hits_in_species=sum(z$is_hit), ions_in_pathway=sum(z$in_pathway), ions_in_species=nrow(z))
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(kegg.pathways, by="KEGG_PATHWAY") %>%
    dplyr::mutate(
      is_significant=dplyr::case_when(
        ions_in_pathway < 2 ~ "No (to few metabolites detected in pathway)",
        pvalue<=0.1 ~ "Yes", 
        T ~ "No")) %>%
    dplyr::arrange(pvalue>0.1, odds<1, pvalue) %>%
    dplyr::select(Species, Drug, KEGG_PATHWAY, KEGG_PATHWAY_NAME, is_significant, odds, pvalue, hits_in_pathway, hits_in_species, ions_in_pathway, ions_in_species)

  readr::write_delim(export.enrichment, "reports/exp9metconcentration_pathway_enrichment.tsv", delim="\t", na="")

  pdf("reports/exp9metconcentration_pathway_enrichment.pdf", width=8, height=8)
  export.enrichment.f = export.enrichment %>%
    dplyr::group_by(Drug, KEGG_PATHWAY_NAME) %>%
    dplyr::filter(sum(hits_in_pathway>0)>1) %>%
    dplyr::ungroup()
  ggplot(export.enrichment.f) +
    geom_boxplot(aes(y=hits_in_pathway/ions_in_pathway, x=KEGG_PATHWAY_NAME)) +
    facet_wrap(~Drug) +
    coord_flip()
  dev.off()
  
  
  #
  # PCA or drug/species data
  #
  pca_n = 3
  metabolomics.pca_species_colors = c("Clostridium saccharolyticum"="#F8766D","Escherichia coli ED1a"="#C49A00","Escherichia coli IAI1"="#53B400","Lactobacillus plantarum WCFS1"="#00C094","Lactococcus lactis IL1403"="#00B6EB","none"="#666666","Streptococcus salivarius"="#FB61D7")
  metabolomics.pca_species_shapes = c("Clostridium saccharolyticum"=21,"Escherichia coli ED1a"=13,"Escherichia coli IAI1"=22,"Lactobacillus plantarum WCFS1"=23,"Lactococcus lactis IL1403"=24,"none"=1,"Streptococcus salivarius"=25)
  metabolomics.pca_species_pch = c("Clostridium saccharolyticum"="\u2776","Escherichia coli ED1a"="\u2777","Escherichia coli IAI1"="\u2778","Lactobacillus plantarum WCFS1"="\u2779","Lactococcus lactis IL1403"="\u277A","none"="\u2B24","Streptococcus salivarius"="\u277B")
  metabolomics.pca = metabolomics.meta %>%
    dplyr::filter(!is.na(SPECIES)) %>%
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, INTERACTION) %>% #  
    dplyr::do((function(z){
      z.data = metabolomics.data.no_drugs2
      z.data = z.data[match(as.character(z$SAMPLE_INDEX), rownames(z.data)),]
      z.data = z.data[,colMeans(z.data>0)==1]
      z.pca = prcomp(z.data)
      z.var = summary(z.pca)$importance[2,]
      cbind(
        DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY=z$DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY, 
        BIOLOGICAL_REPLICATE=z$BIOLOGICAL_REPLICATE,
        data.frame(z.pca$x[,1:pca_n]), 
        data.frame(t(z.var[1:pca_n])) %>% setNames(paste0(names(.), "_var")),
        PC_var=sum(z.var[1:pca_n])
      )
      #data.frame(PCA=z.pca$x[,1], DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY=z$DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY)
    })(.)) %>%
    data.frame() %>% 
    dplyr::mutate(SHAPE=metabolomics.pca_species_shapes[SPECIES], PCH=metabolomics.pca_species_pch[SPECIES]) %>%
    dplyr::filter(SPECIES!="none")
  metabolomics.drug_pca = metabolomics.meta %>%
    dplyr::filter(!is.na(SPECIES) & SPECIES != "none" & DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY %in% c(0, 50)) %>%
    dplyr::group_by(DIFFERENTIAL_TREATMENT_OR_CONDITION) %>% #  
    dplyr::do((function(z){
      #zz<<-z
      z.data = metabolomics.data.no_drugs2
      z.data = z.data[match(as.character(z$SAMPLE_INDEX), rownames(z.data)),]
      z.data = z.data[,colMeans(z.data>0)>=1]
      z.pca = prcomp(z.data)
      z.var = summary(z.pca)$importance[3,]
      cbind(
        z, 
        data.frame(z.pca$x[,1:pca_n]), 
        data.frame(t(z.var[1:pca_n])) %>% setNames(paste0(names(.), "_var")),
        PC_var=sum(z.var[1:pca_n])
      )
    })(.)) %>%
    data.frame() %>% 
    dplyr::mutate(SHAPE=metabolomics.pca_species_shapes[SPECIES], PCH=metabolomics.pca_species_pch[SPECIES])

  #
  # concentration significance
  #
  pca_n = 2
  metabolomics.pca.f = metabolomics.pca %>%
    dplyr::filter(!is.na(SPECIES) & SPECIES != "none" & DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY %in% c(0, 50))
  metabolomics.pca_dist = metabolomics.pca.f %>%
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, INTERACTION, PC_var) %>%
    dplyr::do((function(z){
      zz<<-z
      z.pca = z[,paste0("PC", 1:pca_n), drop=F]
      z.var = z[,paste0("PC", 1:pca_n, "_var"), drop=F]
      z.dist = as.matrix(dist(z.pca*z.var))
      z.dist = z.dist[lower.tri(z.dist)]
      z.same = as.matrix(dist(z$DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY)) > 0
      z.same = ifelse(z.same[lower.tri(z.same)], "different", "same")
      z.ret = data.frame(
        Same_Concentration=z.same,
        Same_Concentration_display=paste0(z.same, " (", ifelse(z$INTERACTION[1]=="Yes", "interaction", "no interaction"), ")"),
        dist=z.dist
      )
      
      # sample_n = 100000
      # z.different = sample(z.ret %>% dplyr::filter(Same_Concentration=="different") %>% .$dist, sample_n, replace=T)
      # z.same = sample(z.ret %>% dplyr::filter(Same_Concentration=="same") %>% .$dist, sample_n, replace=T)
      # z.ret$pvalue = mean(z.different>z.same)
      # z.ret$pvalue = mean(unlist(parallel::mclapply(1:10000, FUN=function(x) {
      #   xx<<-x
      #   sample(z.ret %>% dplyr::filter(Same_Concentration=="different") %>% .$dist, 1) > sample(z.ret %>% dplyr::filter(Same_Concentration=="same") %>% .$dist, 1)
      # })))
      
      z.ret
    })(.)) %>%
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION) %>%
    dplyr::mutate(n=n(), diff=median(dist[Same_Concentration=="different"])-median(dist[Same_Concentration=="same"]), pvalue=wilcox.test(dist[Same_Concentration=="same"], dist[Same_Concentration=="different"])$p.value, pvalue.adj=pvalue*n()) %>%
    dplyr::mutate(
      pvalue_display=paste0("diff: ", sprintf("%0.2f", diff), " (p-value: ", ifelse(pvalue<0.01, sprintf("%1.0e", pvalue), sprintf("%0.2f", pvalue)), ")"),
      pvalue_display.adj=paste0("diff: ", sprintf("%0.2f", diff), " (p-value: ", ifelse(pvalue.adj<0.01, sprintf("%1.0e", pvalue.adj), sprintf("%0.2f", pvalue.adj)), ")"))

  pdf("reports/exp9metconcentration_pca_distance.pdf", width=8, height=8)
  ggplot(metabolomics.pca_dist, aes(x=dist, y=SPECIES)) +
    ggridges::geom_density_ridges(aes(fill=Same_Concentration_display), point_shape=21, alpha=0.8, scale=0.8, jittered_points=T, position=ggridges::position_raincloud(height=0.1)) +
    geom_text(aes(x=9, label=pvalue_display), hjust=1, position=position_nudge(y=0.3), colour="red", size=4, data=metabolomics.pca_dist %>% dplyr::distinct(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, diff, pvalue, pvalue_display)  %>% dplyr::filter(pvalue<=0.05)) +
    scale_fill_manual(values=c("different (interaction)"="#0272A5", "different (no interaction)"="#8BDBF7", "same (interaction)"="#FA3C0D", "same (no interaction)"="#FFBCBA")) +
    guides(fill=guide_legend(ncol=2, title.position="top"), shape=F, color=F) +
    labs(x=paste0("Eucledian distance in PC1-", pca_n, "space\n(median variance: ", round(median(metabolomics.pca_dist$PC_var), 2), ")"), y="", fill="Concentration (0uM, 50uM)") +
    facet_wrap(~DIFFERENTIAL_TREATMENT_OR_CONDITION) +
    theme_classic(base_size=18) +
    theme(axis.title=element_text(size=15), axis.text=element_text(size=13),
                  legend.title=element_text(size=18), legend.text=element_text(size=13), legend.position="bottom",
                  strip.text = element_text(size=20),
                  strip.background=element_blank(), 
                  panel.grid=element_blank())
  dev.off()
  
  
  
  pdf("reports/exp9metconcentration_pca_selected.pdf", width=16, height=8)
  p.all = list()
  for(drug in c("Duloxetine", "Roflumilast")) {
    metabolomics.drug_pca.f = metabolomics.drug_pca %>% dplyr::filter(DIFFERENTIAL_TREATMENT_OR_CONDITION==drug)
    p.all[[drug]] = ggplot(data=metabolomics.drug_pca.f) +
      geom_point(aes(x=PC2, y=PC1, shape=SPECIES, fill=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY), color=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY)), size=5, alpha=0.7) +
      labs(fill="Species", size="Concentration", alpha="Concentration") +
      scale_fill_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#891B00")) +
      scale_color_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#891B00")) +
      scale_shape_manual(name="",values=metabolomics.pca_species_shapes) +
      coord_fixed(ratio=with(metabolomics.drug_pca.f, PC1_var/PC2_var)[1]) +
      labs(
        fill="Concentration", color="Concentration", alpha="Concentration", 
        shape="Species", 
        x=paste0("PC2 (", round(metabolomics.drug_pca.f$PC2_var[1]*100, 1),"%)"), 
        y=paste0("PC1 (", round(metabolomics.drug_pca.f$PC1_var[1]*100, 1),"%)")) +
      theme_classic(base_size=18) +
      theme(strip.text.y=element_text(angle = 180))
  }
  
  # Duloxetine/C. saccharolyticum
  metabolomics.pca_selected = metabolomics.pca %>% dplyr::filter(grepl("sacch", SPECIES) & grepl("Dulox", DIFFERENTIAL_TREATMENT_OR_CONDITION))
  p2 = ggplot(data=metabolomics.pca_selected) +
    geom_point(aes(x=PC2, y=PC1, shape=SPECIES, fill=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY), color=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY)), size=5, alpha=0.7) +
    scale_color_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#891B00")) +
    scale_fill_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#891B00")) +
    scale_shape_manual(name="",values=metabolomics.pca_species_shapes) +
    labs(color="Concentration", 
         x=paste0("PC2 (", round(metabolomics.pca_selected$PC2_var[1]*100, 1),"%)"), 
         y=paste0("PC1 (", round(metabolomics.pca_selected$PC1_var[1]*100, 1),"%)")) +
    guides(size=F, color=F, fill=F, shape=F) +
    coord_fixed(ratio=with(metabolomics.pca_selected, PC1_var/PC2_var)[1]) +
    theme_classic(base_size=18)
  
  # Roflumilast/C. saccharolyticum
  metabolomics.pca_selected = metabolomics.pca %>% dplyr::filter(grepl("sacch", SPECIES) & grepl("Roflumilast", DIFFERENTIAL_TREATMENT_OR_CONDITION))
  p3 = ggplot(data=metabolomics.pca_selected) +
    geom_point(aes(x=PC2, y=PC1, shape=SPECIES, fill=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY), color=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY)), size=5, alpha=0.7) +
    scale_color_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#891B00")) +
    scale_fill_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#891B00")) +
    scale_shape_manual(name="",values=metabolomics.pca_species_shapes) +
    labs(color="Concentration", 
         x=paste0("PC2 (", round(metabolomics.pca_selected$PC2_var[1]*100, 1),"%)"), 
         y=paste0("PC1 (", round(metabolomics.pca_selected$PC1_var[1]*100, 1),"%)")) +
    guides(size=F, color=F, fill=F, shape=F) +
    coord_fixed(ratio=with(metabolomics.pca_selected, PC1_var/PC2_var)[1]) +
    theme_classic(base_size=18)
  
  gridExtra::grid.arrange(p.all[["Duloxetine"]], p2, p3, nrow=1, widths=c(3.4, 1, 1.8))
  dev.off()
  
  pdf("reports/exp9metconcentration_pca.pdf", width=10, height=10)
  ggplot(data=metabolomics.pca) +
    geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=INTERACTION), data=metabolomics.pca %>% dplyr::distinct(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, INTERACTION), alpha=0.1) +
    geom_point(aes(x=PC1, y=PC2, color=factor(SPECIES), size=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY), alpha=factor(DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY))) +
    scale_color_manual(values=metabolomics.pca_species_colors) +
    scale_size_manual(values=c("0"=2, "12.5"=4, "25"=6, "50"=12)) +
    scale_alpha_manual(values=c("0"=1, "12.5"=0.7, "25"=0.5, "50"=0.3)) +
    #scale_color_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#891B00")) +
    scale_fill_manual(values=c(No="#FFFFFF00", Yes="#E41A1C")) +
    labs(color="Species", size="Concentration", alpha="Concentration") +
    scale_y_continuous(position="right") + #, breaks=scales::pretty_breaks(4), limits=c(NA, NA)
    facet_grid(SPECIES~ DIFFERENTIAL_TREATMENT_OR_CONDITION, scales = "free", switch = "y")+
    theme_bw(base_size=18) +
    theme(strip.text.y=element_text(angle = 180))
  dev.off()
  
  #
  # PCA or drug/species data
  #
  pdf("reports/exp9metconcentration_correlated_ions.pdf", width=10, height=10)
  metabolomics.cor = metabolomics.meta %>%
    dplyr::filter(!is.na(SPECIES)) %>%
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, BIOLOGICAL_REPLICATE, INTERACTION) %>% # DIFFERENTIAL_TREATMENT_OR_CONDITION, 
    dplyr::do((function(z){
      z.data = metabolomics.data[match(as.character(z$SAMPLE_INDEX), rownames(metabolomics.data)),]
      z.cor = apply(z.data, 2, function(y) { cor(z$DIFFERENTIAL_TREATMENT_CONCENTRATION_OR_INTENSITY, y, method="spearman") })
      z.cor_sum = sum(abs(z.cor)>0.99, na.rm=T)
      cbind(z[1,,drop=F], CORSUM=z.cor_sum)
    })(.)) %>%
    data.frame() %>% 
    dplyr::filter(SPECIES!="none")  %>% 
    dplyr::group_by(SPECIES, DIFFERENTIAL_TREATMENT_OR_CONDITION, INTERACTION) %>%
    dplyr::summarise(CORSUM.sd=sd(CORSUM), CORSUM=mean(CORSUM))
  
  ggplot(metabolomics.cor) +
    geom_bar(aes(y=CORSUM, x=SPECIES, fill=INTERACTION), stat="identity") +
    facet_wrap(~DIFFERENTIAL_TREATMENT_OR_CONDITION) +
    coord_flip() +
    labs(color="Concentration", x="", y="No. of correlated ions") +
    scale_fill_manual(values=c(No="#666666", Yes="#E41A1C")) +
    theme_bw(base_size=18)
  dev.off()
  
}