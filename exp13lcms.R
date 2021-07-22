library(dplyr)
library(readr)
library(tibble)
library(reshape2)
library(ggplot2)
library(factoextra) 
library(FactoMineR)
library(VennDiagram)
library(ComplexHeatmap)
library(CHNOSZ)
require(XML)
library(Hotelling)
library(shiny)
library(Rdisop)

pca_species_colors = c("Clostridium saccharolyticum"="#F8766D","Escherichia coli ED1a"="#C49A00","Escherichia coli IAI1"="#53B400","Lactobacillus plantarum WCFS1"="#00C094","Lactococcus lactis IL1403"="#00B6EB","none"="#666666","-"="#666666","Streptococcus salivarius"="#FB61D7")
pca_species_shapes = c("Clostridium saccharolyticum"=21,"Escherichia coli ED1a"=13,"Escherichia coli IAI1"=22,"Lactobacillus plantarum WCFS1"=23,"Lactococcus lactis IL1403"=24,"none"=1,"-"=1,"Streptococcus salivarius"=25)
pca_species_pch = c("Clostridium saccharolyticum"="\u2776","Escherichia coli ED1a"="\u2777","Escherichia coli IAI1"="\u2778","Lactobacillus plantarum WCFS1"="\u2779","Lactococcus lactis IL1403"="\u277A","-"="\u2B24","none"="\u2B24","Streptococcus salivarius"="\u277B")
pca_species_pch2 = c("Clostridium saccharolyticum"="\u25CF","Escherichia coli ED1a"="\u2B53","Escherichia coli IAI1"="\u25FC","Lactobacillus plantarum WCFS1"="\u25C6","Lactococcus lactis IL1403"="\u25B2","-"="\u2217","none"="\u2217","Streptococcus salivarius"="\u25BC")

exp13lcms.preprocess_mass = function() {
  
  ppm_extended = 100e-6
  ppm = 10e-6
  
  #
  # Calculate expected m/z peaks for roflumilast and duloxetine considering possible adducts and isotopes
  #
  electron_mass = 0.00054858026 
  
  compounds = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/Table17_test_List.txt", locale=locale(encoding="UTF-8"))
  adducts = readr::read_delim("data/db/adducts.tsv", "\t", comment="#") %>%
    tidyr::separate(adduct_formula, into=c("adduct_formula_add", "adduct_formula_sub"), sep="-") %>%
    dplyr::filter(adduct_name=="[M-H]-")
  
  drug_peaks = compounds %>%
    dplyr::mutate(Formula=gsub(" ", "", Formula)) %>%
    tidyr::crossing(adducts) %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z) {
      zz<<- z
      z.molecule = Rdisop::getMolecule(z$Formula)
      z.exact_mass = z.molecule$isotopes[[1]][1,1]
      
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
      
      cbind(data.frame(z[c("CompoundNo", "Compound", "Mass", "RetentionTime", "Formula", "METLIN_ID")]), exact_ion_mz=ion_mz[1], exact_mass=z.exact_mass)
    })(.)) %>%
    dplyr::mutate(delta_ppm=(Mass-exact_mass)/exact_mass*1e6)
  
  readr::write_tsv(drug_peaks, "C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/Table17_test_List_processed.txt")
  
}

exp13lcms.preprocess_hmdb = function() {
  hmdb_tree2 = xml2::read_xml("data/db/hmdb_metabolites.xml")
  hmdb_nodes_ids = xml2::xml_find_all(hmdb_tree2, "/*[local-name()='hmdb']/*[local-name()='metabolite']/*[local-name()='accession']/text()")
  hmdb_ids = xml2::xml_text(hmdb_nodes_ids)
  hmdb_nodes_formula = xml2::xml_find_all(hmdb_tree2, "/*[local-name()='hmdb']/*[local-name()='metabolite']/*[local-name()='chemical_formula']/text()")
  hmdb_formula = xml2::xml_text(hmdb_nodes_formula)
  hmdb_nodes_smiles = xml2::xml_find_all(hmdb_tree2, "/*[local-name()='hmdb']/*[local-name()='metabolite']/*[local-name()='smiles']/text()")
  hmdb_smiles = xml2::xml_text(hmdb_nodes_smiles)
  hmdb_nodes_inchi = xml2::xml_find_all(hmdb_tree2, "/*[local-name()='hmdb']/*[local-name()='metabolite']/*[local-name()='inchi']/text()")
  hmdb_inchi = xml2::xml_text(hmdb_nodes_inchi)

  # hmdb_nodes_metlin = xml2::xml_find_all(hmdb_tree2, "/*[local-name()='hmdb']/*[local-name()='metabolite']/*[local-name()='metlin_id']/text()")
  # hmdb_metlin = xml2::xml_text(hmdb_nodes_metlin)
  # xml2::write_xml(xml2::xml_find_all(hmdb_tree2, "/*[local-name()='hmdb']/*[local-name()='metabolite']")[[1]], file="test.xml")
  
  hmdb_df = data.frame(HMDB_ID=hmdb_ids, HMDB_FORMULA=hmdb_formula, HMDB_SMILES=hmdb_smiles, HMDB_INCHI=hmdb_inchi, stringsAsFactors=F) %>%
    dplyr::filter(grepl("^([A-Z][a-z]?[0-9]*)+$", HMDB_FORMULA))
  elements = c("C", "H", "D", "F", "Cl", "Br", "Fe", "As", "Mg", "Gd", "N", "Na", "Bi", "B", "I", "O", "P", "S", "Si", "Se", "Sn")
  hmdb_df.formulas =   do.call(rbind, CHNOSZ::makeup(hmdb_df$HMDB_FORMULA, count.zero=T))
  hmdb_df$HMDB_FORMULA2 = apply(hmdb_df.formulas[,intersect(elements, colnames(hmdb_df.formulas))], 1, function(z) paste0(names(z)[z>0], ifelse(z[z>0]>1, z[z>0], ""), collapse=" "))
  
  readr::write_tsv(hmdb_df, path="data/db/hmdb_metabolites.tsv", col_names=T, na="")
}

exp13lcms.preprocess_msms_standards = function() {
  data = readr::read_tsv("data/exp13lcms/150samples/3_annotations_ions.tsv", col_names=T, comment="#", col_types=list("PubChemID"="c"))
  x = data %>% 
    dplyr::mutate(HMP_ID=gsub("HMDB", "HMDB00", HMP_ID)) %>%
    dplyr::filter(!is.na(CAS) | ! is.na(PubChemID) | !is.na(HMP_ID)) %>%
    dplyr::distinct(CAS, PubChemID, HMP_ID, Formula)
  msms = readr::read_tsv("data/exp13lcms/msms_metabolites.tsv", col_names=T) 
  msms.updated = msms %>%
    dplyr::left_join(x %>% dplyr::mutate(Formula.CAS=Formula) %>% dplyr::distinct(CAS, Formula.CAS), by=c("CAS")) %>%
    dplyr::left_join(x %>% dplyr::mutate(Formula.HMP=Formula) %>% dplyr::distinct(HMP_ID, Formula.HMP), by=c("MetaboliteID"="HMP_ID")) %>%
    dplyr::mutate(Formula=dplyr::case_when(
      !is.na(Formula)~Formula,
      !is.na(Formula.CAS)~Formula.CAS,
      !is.na(Formula.HMP)~Formula.HMP,
      T~NA_character_)) %>%
    dplyr::select(MetaboliteID, MetaboliteName, CAS, Formula)
  readr::write_tsv(msms.updated, path="data/exp13lcms/msms_metabolites.tsv", col_names=T)
}

exp13lcms.preprocess_228samples = function() {
  msms_validation = readr::read_tsv("data/exp13lcms/150samples/5_msms_long.tsv") %>%
    dplyr::group_by(Compound) %>%
    dplyr::summarise(MSMS_name=paste0(MSMS_name, ifelse(length(unique(MSMS_MZdelta_rank))>1, paste0(" (", MSMS_MZdelta_rank, ")"), ""), collapse="; "), is_msms_validated=T)
  
  data = readr::read_tsv("data/exp13lcms/228samples/raw/228 samples-kiran-sonja-bugdrug.txt", col_names=T, comment="#", col_types=list("PubChemID"="c"))
  
  annotations_ions = data %>% # Formula, MetlinID, PubChemID, KEGG_ID, HMP_ID, Frequency, ChEBI_ID, CAS, Score, 
    dplyr::select(Compound, Mass, RetentionTime, IonSpecies) %>%
    dplyr::mutate(is_peptide=grepl("\\b[A-Z][a-z][a-z] [A-Z][a-z][a-z] [A-Z][a-z][a-z]\\b", Compound)) %>% # peptide
    dplyr::mutate(is_lipid=grepl("[0-9]{1}:[0-9]{1}", Compound)) %>%
    dplyr::mutate(is_main_fragment=!grepl(" Es.-.*", Compound))  %>%
    dplyr::left_join(msms_validation, by="Compound") %>%
    dplyr::mutate(is_msms_validated=!is.na(is_msms_validated))

  # TODO: sample names to info (Kiran)
  data_long = data %>%
    dplyr::select(Compound, dplyr::starts_with("SD004E")) %>%
    reshape2::melt(id.vars="Compound", variable.name="Sample", value.name="Intensity") %>%
    dplyr::mutate(IntensityTransformation=dplyr::case_when(grepl("Log2", Sample)~"log2", grepl("\\(raw\\)", Sample)~"raw", T~"unknown")) %>%
    dplyr::mutate(SampleLong=Sample, Sample=gsub("SD004E_KP_BugDrug_HILICZ_NP_|: Log2|\\(raw\\)", "", Sample)) %>%
    dplyr::mutate(Intensity=dplyr::case_when(
      IntensityTransformation == "raw"~ifelse(Intensity>0, log2(Intensity), 0),
      IntensityTransformation == "log2"~Intensity,
      T ~ NA_real_
    )) %>%
    tidyr::extract(Sample, into=c("Plate", "Well", "Repeated"), regex="P([0-9])_([A-H][0-9]+)((?:R| \\([0-9]\\))?)") %>%
    dplyr::mutate(Plate=paste("Plate", Plate), Repeated=gsub("[^0-9]", "", gsub("R", 10, Repeated))) %>%
    dplyr::mutate(Intensity=as.numeric(Intensity)) %>%
    dplyr::select(-IntensityTransformation)

  
  annotations_samples = readr::read_tsv("data/exp13lcms/228samples/annotations_samples.tsv", col_names=T)
  data_long.repeated = data_long %>%
    dplyr::inner_join(annotations_samples %>% dplyr::select(SampleName, Plate, Well, Species, Drug, Medium, Replicate, ConcentrationFinal), by=c("Plate", "Well")) %>% 
    dplyr::group_by(Species, Drug, Replicate, ConcentrationFinal, Plate, Compound) %>%
    dplyr::arrange(dplyr::desc(Repeated)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  options(digits=6)
  x = data.frame(data.frame(data_long.repeated %>%
    dplyr::mutate(Intensity=ifelse(is.na(Intensity), "", ifelse(Intensity==0, "0.0", as.character(Intensity)))) %>%
    reshape2::dcast(Compound ~ SampleLong, value.var="Intensity") %>%
    dplyr::inner_join(data %>% dplyr::select(-dplyr::starts_with("SD004E")), by="Compound") %>%
    dplyr::mutate(Mass=ifelse(is.na(Mass), "", ifelse(Mass==0, "0.0", as.character(Mass)))) %>%
    dplyr::mutate(RetentionTime=ifelse(is.na(RetentionTime), "", ifelse(RetentionTime==0, "0.0", as.character(RetentionTime)))) %>%
    t(), stringsAsFactors=F) %>% tibble::rownames_to_column("SampleLong") %>%
    dplyr::left_join(data_long.repeated %>% dplyr::distinct(SampleLong, SampleName, Drug, ConcentrationFinal, Plate, Well, Replicate, Medium), by="SampleLong") %>%
    dplyr::select(SampleLong, SampleName, Drug, ConcentrationFinal, Plate, Well, Replicate, Medium, grep("^(SampleLong|Plate|Well|Replicate|Drug)$", colnames(.), invert=T)) %>%
    t(), stringsAsFactors=F)
  x = x[c(1:8, match(data$Compound, x[,1])),] # Order in the same way as original data
  readr::write_tsv(x, path="data/exp13lcms/228samples/raw/228 samples-kiran-sonja-bugdrug_filtered_with_sample_annotation.txt", col_names=F, na="")
  
  readr::write_tsv(data_long.repeated %>% dplyr::select(-Species, -Drug, -Replicate, -Medium, -Repeated, -ConcentrationFinal, -SampleName), path="data/exp13lcms/228samples/data_long.tsv", col_names=T)
  readr::write_tsv(annotations_ions, path="data/exp13lcms/228samples/annotations_ions.tsv", col_names=T)
}

exp13lcms.preprocess_150samples = function() {
  data = readr::read_tsv("data/exp13lcms/150samples/raw/4_SD004D_150Samples_modified_with_5000counts_updated_130330.txt", col_names=T, comment="#", col_types=list("PubChemID"="c"))
  drugbank = readr::read_csv("data/db/drugbank/drugbank_vocabulary.csv", col_names=T, comment="#") %>%
    dplyr::group_by(CAS) %>%
    dplyr::summarise(DrugBankID=paste0(`DrugBank ID`, collapse=","))
  
  kegg2compounds = readr::read_tsv("data/db/KEGG/kegg_compounds.tsv", col_names=T) %>%
    tidyr::separate_rows(KEGG_METABOLITE_NAME, sep=";") %>%
    dplyr::group_by(KEGG_CPD, KEGG_FORMULA) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  kegg2species = readr::read_tsv("data/db/KEGG/kegg_species2cpd.tsv", col_names=T) %>%
    dplyr::mutate(KEGG_CPD=gsub("cpd:", "", KEGG_CPD)) %>%
    dplyr::filter(grepl("sacch", Species) & !is.na(KEGG_FORMULA) & !grepl("X|R|n\\b", KEGG_FORMULA)) %>%
    dplyr::inner_join(kegg2compounds %>% dplyr::distinct(KEGG_CPD, KEGG_METABOLITE_NAME), by="KEGG_CPD") %>%
    dplyr::group_by(KEGG_FORMULA) %>%
    dplyr::summarise(is_saccharolyticum=T, kegg.names=paste(unique(KEGG_METABOLITE_NAME), collapse=";"), is_saccharolyticum.metabolite=substr(kegg.names, 1, pmin(50, nchar(kegg.names)))) %>%
    dplyr::select(-kegg.names)
  elements = c("C", "H", "D", "F", "Cl", "Br", "Fe", "As", "Mg", "Gd", "N", "Na", "Bi", "B", "I", "O", "P", "S", "Si", "Se", "Sn")
  kegg2species.formulas = do.call(rbind, CHNOSZ::makeup(kegg2species$KEGG_FORMULA, count.zero=T))
  kegg2species$Formula = apply(kegg2species.formulas[,intersect(elements, colnames(kegg2species.formulas))], 1, function(z) paste0(names(z)[z>0], ifelse(z[z>0]>1, z[z>0], ""), collapse=" "))
  kegg2species = kegg2species %>% dplyr::filter(!is.na(Formula) & Formula!="") %>% dplyr::distinct(Formula, is_saccharolyticum, is_saccharolyticum.metabolite)
  
  msms = readr::read_tsv("data/exp13lcms/msms_metabolites.tsv", col_names=T) %>%
    dplyr::mutate(is_msms=T) %>%
    dplyr::filter(!is.na(Formula)) %>%
    dplyr::distinct(Formula, is_msms)
  
  
  # The new curated list preservers only metabolites that map with PPM < 20
  msms_validation.new = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp13lcms/msms_data_curated.tsv")
  msms_validation = readr::read_tsv("data/exp13lcms/150samples/5_msms_long.tsv") %>%
    dplyr::group_by(Compound) %>%
    dplyr::inner_join(msms_validation.new %>% dplyr::distinct(Compound, `Mass (Da)`, `mass error (ppm)`), by=c("MSMS_name"="Compound")) %>%
    dplyr::summarise(MSMS_name=paste0(unique(MSMS_name), collapse="; "), is_msms_validated=T)
  
  
  annotations_ions = data %>%
    dplyr::select(Compound, CompoundName, Mass, RetentionTime, IonSpecies, Annotations, Formula, MetlinID, PubChemID, KEGG_ID, HMP_ID, Frequency, ChEBI_ID, CAS, Score) %>%
    dplyr::left_join(drugbank, by="CAS") %>%
    dplyr::left_join(kegg2species, by="Formula") %>%
    dplyr::left_join(msms, by="Formula") %>%
    dplyr::left_join(msms_validation, by=c("Compound")) %>%
    dplyr::mutate(is_saccharolyticum=!is.na(is_saccharolyticum)) %>%
    dplyr::mutate(is_msms=!is.na(is_msms)) %>%
    dplyr::mutate(is_msms_validated=!is.na(is_msms_validated)) %>%
    dplyr::mutate(is_peptide=grepl("\\b[A-Z][a-z][a-z] [A-Z][a-z][a-z] [A-Z][a-z][a-z]\\b", Compound)) %>% # peptide
    dplyr::mutate(is_lipid=grepl("[0-9]{1}:[0-9]{1}", Compound)) %>%
    dplyr::mutate(is_main_fragment=!grepl(" Es.-.*", Compound)) %>%
    dplyr::mutate(is_mapped_formula=!is.na(Formula)) %>%
    dplyr::mutate(is_mapped_metabolite=is_mapped_formula & !grepl("\\b((?:C|Cl|H|N|O|P|S)(?:\\b|[0-9]{1,2})(?: |\\b)){3}", Compound)) 

  data_long = data %>%
    dplyr::select(Compound, dplyr::starts_with("SD004D")) %>%
    reshape2::melt(id.vars="Compound", variable.name="Sample", value.name="Intensity") %>%
    dplyr::mutate(IntensityTransformation=dplyr::case_when(grepl("Log2", Sample)~"log2", grepl("\\(raw\\)", Sample)~"raw", T~"unknown")) %>%
    dplyr::mutate(Sample=gsub("SD004D_KP_HILICZ_NP_BD_Sonja_|: Log2|\\(raw\\)", "", Sample)) %>%
    dplyr::mutate(Intensity=dplyr::case_when(
      IntensityTransformation == "raw"~ifelse(Intensity>0, log2(Intensity), 0),
      IntensityTransformation == "log2"~Intensity,
      T ~ NA_real_
    )) %>%
    tidyr::extract(Sample, into=c("Plate", "Well", "Repeated"), regex="P([12])_([A-H][0-9]+)((?:R| \\([0-9]\\))?)") %>%
    dplyr::mutate(Plate=paste("Plate", Plate), Repeated=gsub("[^0-9]", "", gsub("R", 10, Repeated))) %>%
    dplyr::mutate(Intensity=as.numeric(Intensity)) %>%
    dplyr::select(-IntensityTransformation)
  
  annotations_samples = readr::read_tsv("data/exp13lcms/150samples/annotations_samples.tsv", col_names=T)
  data_long.repeated = data_long %>%
    dplyr::inner_join(annotations_samples %>% dplyr::select(Plate, Well, Species, Drug, Medium, Replicate, ConcentrationFinal, Time), by=c("Plate", "Well")) %>% 
    dplyr::group_by(Species, Drug, Replicate, ConcentrationFinal, Plate, Time, Compound) %>%
    dplyr::arrange(dplyr::desc(Repeated)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-Species, -Drug, -Replicate, -Medium, -Repeated, -ConcentrationFinal, -Time)
  
  readr::write_tsv(data_long.repeated, path="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp13lcms/150samples/data_long.tsv", col_names=T)
  readr::write_tsv(annotations_ions, path="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp13lcms/150samples/annotations_ions.tsv", col_names=T)
}

exp13lcms.228pca = function() 
{
  annotations_samples = readr::read_tsv("data/exp13lcms/228samples/annotations_samples.tsv", col_names=T)
  annotations_ions = readr::read_tsv("data/exp13lcms/228samples/annotations_ions.tsv", col_names=T)
  data_long_raw = readr::read_tsv("data/exp13lcms/228samples/data_long.tsv", col_names=T)
  
  data_long = data_long_raw %>%
    dplyr::inner_join(annotations_samples, by=c("Plate", "Well")) %>%
    dplyr::inner_join(annotations_ions, by="Compound") %>%
    dplyr::group_by(Plate, Compound, Medium, Drug, Species, ConcentrationFinal) %>%
    dplyr::mutate(Intensity.has4replicates=sum(Intensity>0)>=4) %>%
    dplyr::group_by(Compound, Medium, Drug) %>%
    dplyr::mutate(Intensity.has4replicates_all=any(Intensity.has4replicates)) %>%
    dplyr::group_by(Compound, Medium, Drug, Species) %>%
    dplyr::mutate(Intensity.has4replicates_species=any(Intensity.has4replicates)) %>%
    dplyr::ungroup()

  #
  # Find which features are related to drug ions
  #
  data_long.min = min(data_long$Intensity[data_long$Intensity>0])
  data_long.0 = data_long %>% dplyr::filter(ConcentrationFinal==0)
  data_long.diff = data_long %>%
    dplyr::inner_join(data_long.0 %>% dplyr::select(Medium, Drug, Species, Replicate, Compound, Intensity.0=Intensity), by=c("Medium", "Species", "Drug", "Replicate", "Compound")) %>%
    dplyr::group_by(Medium, Species,Drug, Compound, ConcentrationFinal) %>%
    dplyr::do((function(z) {
      zz<<-z
      f1 = z$Intensity>=4
      f0 = z$Intensity.0>=4
      
      if(any(is.na(f1)) | any(is.na(f0))) {
        zz<<-z
        asddsa()
        
      }
      if(sum(f1)<3 & sum(f0)<3) {
        return(data.frame(IntensityPvalue=1, IntensityDiff=0, IntensityCount=sum(f1), Intensity0Count=sum(f0)))
      }
      
      if(sum(f0)>=3) {
        int0 = z$Intensity.0[f0]
      } else {
        int0 = ifelse(f0, z$Intensity.0, data_long.min)
      }
      
      if(sum(f1)>=3) {
        int1 = z$Intensity[f1]
      } else {
        int1 = ifelse(f1, z$Intensity, data_long.min)
      }
      z.test = t.test(int1, int0)
      data.frame(IntensityPvalue=z.test$p.value, IntensityDiff=mean(int1-int0), IntensityCount=sum(f1), Intensity0Count=sum(f0))
    })(.)) %>%
    dplyr::ungroup()  %>%
    dplyr::mutate(is_significant=IntensityPvalue<=0.05 & abs(IntensityDiff)>=log2(2))
  data_long = data_long[,setdiff(colnames(data_long), c("IntensityCount", "IntensityDiff", "IntensityPvalue", "Intensity0Count", "is_significant"))] %>%
    dplyr::left_join(data_long.diff, by=c("Medium", "Species", "Drug", "Compound", "ConcentrationFinal")) %>%
    dplyr::mutate(is_significant=!is.na(is_significant) & is_significant) 
  
  #
  # Filter data for all species and per species PCA
  #
  data_long.samples = data_long %>%
    dplyr::group_by(Compound, Medium, Drug) %>%
    dplyr::filter(!any(is_significant & Species=="-")) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Species!="-")

  # 
  #   PCA
  #   
  data_long.pca_all = data_long.samples %>%
    dplyr::filter(Intensity.has4replicates_all & ConcentrationFinal %in% c(0,50)) %>%
    dplyr::group_by(Medium, Drug) %>%
    dplyr::do((function(z){
      zz<<-z
      z.wide = z %>% reshape2::dcast(Plate + Medium + Drug + Species + ConcentrationFinal + Replicate ~ Compound, value.var="Intensity")
      z.mat = as.matrix(z.wide %>% dplyr::select(-(Plate:Replicate)))
      z.pca = prcomp(z.mat)
      z.eigs = as.data.frame(t(z.pca$sdev^2/sum(z.pca$sdev^2)))
      z.pca_norm = t(t(z.pca$x)*as.numeric(z.eigs))
      colnames(z.eigs) = paste0("PC", 1:ncol(z.eigs), ".var")
      colnames(z.pca_norm) = paste0("PC", 1:ncol(z.pca$x), ".norm")
      z.out = cbind(z.wide %>% dplyr::select(Plate:Replicate), z.pca$x[,1:2], z.eigs[1:2], z.pca_norm[,1:2])
      z.out
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SHAPE=pca_species_shapes[Species], PCH=pca_species_pch2[Species]) %>%
    dplyr::mutate(Variance=paste0(round(PC1.var*100),"% / ", round(PC2.var*100),"%"))

  data_long.pca_species = data_long.samples %>%
    dplyr::filter(Intensity.has4replicates_species) %>%
    dplyr::group_by(Medium, Drug, Species) %>%
    dplyr::do((function(z){
      zz<<-z
      z.wide = z %>% reshape2::dcast(Plate + Medium + Drug + Species + ConcentrationFinal + Replicate ~ Compound, value.var="Intensity")
      z.mat = as.matrix(z.wide %>% dplyr::select(-(Plate:Replicate)))
      z.pca = prcomp(z.mat)
      z.eigs = as.data.frame(t(z.pca$sdev^2/sum(z.pca$sdev^2)))
      z.pca_norm = t(t(z.pca$x)*as.numeric(z.eigs))
      colnames(z.eigs) = paste0("PC", 1:ncol(z.eigs), ".var")
      colnames(z.pca_norm) = paste0("PC", 1:ncol(z.pca$x), ".norm")
      z.out = cbind(z.wide %>% dplyr::select(Plate:Replicate), z.pca$x[,1:2], z.eigs[1:2], z.pca_norm[,1:2])
      z.out
    })(.)) %>%
    dplyr::mutate(SHAPE=pca_species_shapes[Species], PCH=pca_species_pch2[Species]) %>%
    dplyr::mutate(Condition=paste0(Drug, " (", Species, ")"), Variance=paste0(round(PC1.var*100),"% / ", round(PC2.var*100),"%"))

  
  cairo_pdf("reports/exp13lcms_228pca.pdf", width=8, height=8,onefile=T)
  for(drug in unique(data_long.pca_all$Drug)) {
    pca_species_pch2.x = pca_species_pch2[names(pca_species_pch2)!="-"]
    data_long.pca_all.f = data_long.pca_all %>% dplyr::filter(Drug==drug & ConcentrationFinal %in% c(0,50) & Species!="-")
    data_long.pca_all.f.zoom = data_long.pca_all.f %>% 
      dplyr::mutate(SpeciesGroup=dplyr::case_when(grepl("coli", Species)~"E. coli", grepl("Lacto", Species)~"LAC", T~"Other")) %>% 
      dplyr::group_by(SpeciesGroup) %>%
      dplyr::do((function(z) data.frame(PC1=range(z$PC1), PC2=range(z$PC2)))(.))
    
    p1 = ggplot(data_long.pca_all.f) +
      geom_text(aes(PC2, PC1, color=ConcentrationFinal, label=PCH, alpha=ConcentrationFinal), size=6, family="Lucida Sans Unicode") +
      labs(color="Concentration", title=paste0(drug, ":\n", paste0(pca_species_pch2.x[1:3], " ",  names(pca_species_pch2.x[1:3]), collapse=" "), "\n", paste0(pca_species_pch2.x[4:7], " ",  names(pca_species_pch2.x[4:7]), collapse=" "))) +
      scale_color_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#CD1C20")) +
      scale_alpha_manual(values=c("0"=0.7, "12.5"=0.7, "25"=0.7, "50"=0.9)) +
      coord_fixed(ratio=with(data_long.pca_all.f, PC1.var/PC2.var)[1]) +
      labs(
        fill="Concentration", color="Concentration", alpha="Concentration", 
        x=paste0("PC2 (", round(data_long.pca_all.f$PC2.var[1]*100, 1),"%)"), 
        y=paste0("PC1 (", round(data_long.pca_all.f$PC1.var[1]*100, 1),"%)")) +
      theme_classic(base_size=18) +
      theme(strip.text.y=element_text(angle = 180), legend.text=element_text(family="Lucida Sans Unicode"), plot.title=element_text(size=10))
    p2 = with(data_long.pca_all.f.zoom %>% dplyr::filter(SpeciesGroup=="E. coli"), p1 + coord_cartesian(ylim=PC1, xlim=PC2) + labs(title="E. coli"))
    p3 = with(data_long.pca_all.f.zoom %>% dplyr::filter(SpeciesGroup=="LAC"), p1 + coord_cartesian(ylim=PC1, xlim=PC2) + labs(title="LAC"))
    print(p1) 
    print(p2)
    print(p3)
  }
  
  data_long.pca_species.f = data_long.pca_species %>% dplyr::filter(grepl("sacch", Species) & Drug=="Duloxetine")
  p = ggplot(data_long.pca_species.f) +
    geom_text(aes(PC2, PC1, color=ConcentrationFinal, label=PCH, alpha=ConcentrationFinal), size=4, family="Lucida Sans Unicode", show.legend=F) +
    labs(color="Concentration") +
    scale_color_manual(values=c("0"="#999999", "12.5"="#B27F72", "25"="#FF3200", "50"="#CD1C20")) +
    scale_alpha_manual(values=c("0"=0.7, "12.5"=0.7, "25"=0.7, "50"=0.9)) +
    coord_fixed(ratio=with(data_long.pca_species.f, PC1.var/PC2.var)[1]) +
    labs(
      title="Clostridium saccharolyticum",
      fill="Concentration", color="Concentration", alpha="Concentration", 
      x=paste0("PC2 (", round(data_long.pca_species.f$PC2.var[1]*100, 1),"%)"), 
      y=paste0("PC1 (", round(data_long.pca_species.f$PC1.var[1]*100, 1),"%)")) +
    theme_classic(base_size=18) +
    theme(strip.text.y=element_text(angle = 180))
  print(p)
  dev.off()
  
  data_long.export = data_long %>%
    dplyr::filter(Intensity.has4replicates_all & ConcentrationFinal %in% c(0,50)) %>%
    dplyr::mutate(UniqueID=paste(Mass, RetentionTime)) %>%
    dplyr::mutate(Species=ifelse(Species=="-", "Control", Species)) %>%
    dplyr::mutate(Species=factor(Species, unique(Species))) %>%
    dplyr::filter(Drug=="Duloxetine") %>%
    dplyr::arrange(Species, ConcentrationFinal, Replicate) %>%
    dplyr::mutate(Sample=paste0(ConcentrationFinal, "uM (", Species, ", ", Replicate, ")"), Sample=factor(Sample, unique(Sample))) %>%
    reshape2::dcast(UniqueID ~ Sample, value.var="Intensity")
  readr::write_tsv(data_long.export, path="reports/exp13lcms_228pca.txt")    
  
  
  data_long.pca_all.test = data_long.pca_all %>%
    reshape2::melt(measure.vars=c("PC1", "PC2")) %>%
    reshape2::dcast(Species+Drug+Replicate+ConcentrationFinal~variable, value.var="value") %>%
    dplyr::group_by(Species, Drug) %>%
    dplyr::do((function(z){
      zz <<- z
      z.0 = z %>% dplyr::filter(ConcentrationFinal==0) %>% dplyr::select(PC1, PC2)
      z.50 = z %>% dplyr::filter(ConcentrationFinal==50) %>% dplyr::select(PC1, PC2)
      z.dist = z %>% 
        dplyr::group_by(ConcentrationFinal) %>% dplyr::summarise_all(mean) %>%
        dplyr::select(PC1, PC2) %>%
        dist() %>% as.vector()
      z.test = Hotelling::hotelling.test(z.0, z.50)
      data.frame(pval=z.test$pval, dist=z.dist)
    })(.))
  readr::write_tsv(data_long.pca_all.test, path="reports/exp13lcms_228pca_test.tsv", na="")
}

validate.msms_retentions = function() {
  msms_validation = readr::read_tsv("data/exp13lcms/150samples/5_msms_long.tsv") %>%
    dplyr::group_by(Compound) %>%
    dplyr::summarise(MSMS_name=paste0(MSMS_name, collapse="; "), MSMS_RetentionTime=`Retention Time`[1], MSMS_Mass=Mass[1]) %>%
    dplyr::ungroup() %>%
    dplyr::rename(MSMS_Compound="Compound")
  hits = readr::read_tsv("reports/exp13lcms_candidates_lcms_all.tsv") %>%
    dplyr::distinct(Compound, Mass, RetentionTime)
  all = annotations_ions %>%
    dplyr::distinct(Compound, Mass, RetentionTime)
  
  step=0.05; step_robustness=1.2; rt=seq(1.2, 19, by=step)
  x = data.frame(MSMS_RetentionTime.1=rt) %>%
    tidyr::crossing(msms_validation) %>%
    dplyr::group_by(MSMS_RetentionTime.1) %>%
    dplyr::summarise(n.msms=as.numeric(any(step*step_robustness>=abs(MSMS_RetentionTime.1-MSMS_RetentionTime))))
  y = data.frame(RetentionTime.1=rt) %>%
    tidyr::crossing(hits) %>%
    dplyr::group_by(RetentionTime.1) %>%
    dplyr::summarise(n.hits=as.numeric(any(step*step_robustness>=abs(RetentionTime.1-RetentionTime))))
  xy = x %>% dplyr::inner_join(y, by=c("MSMS_RetentionTime.1"="RetentionTime.1"))
  xy.ccf = ccf(x$n.msms, y$n.hits, na.action=na.omit)
  table(xy$n.msms, xy$n.hits)
  fisher.test(table(xy$n.msms, xy$n.hits))
  xy.ccf$lag[which.max(xy.ccf$acf)] * diff(x$MSMS_RetentionTime.1)[1]
  
  
  ggplot(xy %>% reshape2::melt(measure.vars=c("n.hits", "n.msms"))) +
    geom_bar(aes(x=MSMS_RetentionTime.1, y=value, fill=variable), alpha=0.5, stat="identity") +
    scale_x_log10()
  
}

exp13lcms.150analyze = function() 
{
  annotations_samples = readr::read_tsv("data/exp13lcms/150samples/annotations_samples.tsv", col_names=T)
  annotations_ions = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp13lcms/150samples/annotations_ions.tsv", col_names=T)
  data_long_raw = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp13lcms/150samples/data_long.tsv", col_names=T)
  
  data_long = data_long_raw %>%
    dplyr::inner_join(annotations_samples, by=c("Plate", "Well")) %>%
    dplyr::inner_join(annotations_ions, by="Compound") %>%
    dplyr::filter(SampleType!="-" & Time=="2 days")  %>%  # Plate=="Plate 1
    dplyr::group_by(Plate, Compound, Medium, Drug, Species, ConcentrationFinal) %>%
    dplyr::mutate(Intensity.has4replicates=sum(Intensity>0)>=4) %>%
    # dplyr::group_by(Compound, Medium, Drug, Species) %>%
    # dplyr::filter(any(Intensity.has4replicates)) %>%
    dplyr::ungroup()  %>%
    dplyr::mutate(UniqueID=paste(Mass, RetentionTime))
  
  #
  # Differential analysis
  #
  intensity_min = min(data_long$Intensity[data_long$Intensity>0])
  data_long.diff = data_long %>% 
    dplyr::filter(dplyr::between(as.numeric(ConcentrationFinal), 10, 70)) %>% 
    dplyr::inner_join(data_long %>%  dplyr::filter(ConcentrationFinal==0) %>% dplyr::select(Species, Medium, Drug, Replicate, Compound, Intensity.0=Intensity, Well.0=Well), by=c("Species", "Medium", "Drug", "Replicate", "Compound")) %>%
    # bind_cols(group_i = group_indices(., Species, Medium, Drug, Compound, ConcentrationFinal)) %>%
    # dplyr::filter(group_i<10) %>%
    dplyr::group_by(Species, Medium, Drug, Compound, UniqueID, ConcentrationFinal, Intensity.has4replicates) %>%
    dplyr::do((function(z) {
      
      z.f = z %>% dplyr::filter(Intensity>0 & Intensity.0>0)
      if(nrow(z.f) >= 2) {
        z.test = t.test(z.f$Intensity, z.f$Intensity.0)
        return(data.frame(pval=z.test$p.value, diff=mean(z.f$Intensity)-mean(z.f$Intensity.0)))
      }
      z.f = z %>% dplyr::filter(Intensity>0 | Intensity.0>0)
      if(nrow(z.f) >= 2) {
        z.intensity = pmax(jitter(rep(intensity_min,nrow(z.f))), z.f$Intensity)
        z.intensity0 = pmax(jitter(rep(intensity_min,nrow(z.f))), z.f$Intensity.0)
        z.test = t.test(z.intensity, z.intensity0)
        return(data.frame(pval=ifelse(nrow(z.f)>2, z.test$p.value, 1), diff=mean(z.intensity)-mean(z.intensity0)))
      }
      
      return(data.frame(pval=1, diff=0))
    })(.)) %>%
    dplyr::mutate(is_significant=pval<0.01 & abs(diff)>=3) %>%
    dplyr::ungroup()
  
  data_long.diff2 = data_long.diff %>%
    dplyr::arrange(dplyr::desc(ConcentrationFinal)) %>%
    dplyr::group_by(Species, Medium, Drug, UniqueID) %>%
    dplyr::mutate(
      total_significant=sum(is_significant), 
      total_significant_consecutive=with(rle(is_significant), ifelse(values[1], lengths[1], 0)), 
      total_concentrations=n()) %>%
    dplyr::group_by(Medium, Drug, UniqueID) %>%
    dplyr::arrange(dplyr::desc(ConcentrationFinal)) %>%
    dplyr::mutate(
      total_significant_sacch=sum(is_significant[grepl("sacch", Species)]),
      total_significant_consecutive_sacch=with(rle(is_significant[grepl("sacch", Species)]), ifelse(values[1], lengths[1], 0)),
      total_significant_none=sum(is_significant[Species=="no bug"])) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(annotations_ions, by="Compound")

  data_long = data_long[,setdiff(colnames(data_long), c("is_significant", "pval", "diff"))] %>%
    dplyr::left_join(data_long.diff2 %>% dplyr::distinct(Medium, Species, Drug, Compound, ConcentrationFinal, is_significant, pval, diff), by=c("Medium", "Species", "Drug", "Compound", "ConcentrationFinal")) %>%
    dplyr::mutate(is_significant=!is.na(is_significant) & is_significant)
  
  
  
  # 
  #   PCA
  #   
  d = data_long %>% 
    dplyr::group_by(UniqueID) %>%
    dplyr::mutate(is_any_significant=ifelse(any(is_significant), "Yes", "No"), is_significant=ifelse(is_significant, "Yes", "No"))  %>%
    dplyr::ungroup()
  
  #
  # Wide version of the data
  # 
  data_long.diff2.diff = data_long.diff2 %>%
    dplyr::mutate(ConcentrationFinal=paste0("Effect size (", ConcentrationFinal, "uM)")) %>%
    reshape2::dcast(Medium + Drug + Species + UniqueID ~ ConcentrationFinal, value.var="pval")
  data_long.diff2.pval = data_long.diff2 %>%
    dplyr::mutate(ConcentrationFinal=paste0("p-value (", ConcentrationFinal, "uM)")) %>%
    reshape2::dcast(Medium + Drug + Species + UniqueID ~ ConcentrationFinal, value.var="diff")

  d.f = d %>%
    dplyr::filter(grepl("sacch", Species)) %>%
    dplyr::mutate(Sample=paste0("Intensity (", ConcentrationFinal, "uM, ", Replicate, ")")) %>%
    dplyr::mutate(SampleLong=stringr::str_glue("SD004D_KP_HILICZ_NP_BD_Sonja_P{plate}_{well}(raw)", plate=gsub("Plate ", "", Plate), well=Well))
  d.wide = d.f %>%
    reshape2::dcast(Medium + Drug + Species + MSMS_name + UniqueID + is_peptide + is_mapped_metabolite + is_lipid + is_main_fragment + is_any_significant ~ Sample , value.var="Intensity") %>%
    dplyr::inner_join(data_long.diff2.pval, by=c("Medium", "Drug", "Species", "UniqueID")) %>%
    dplyr::inner_join(data_long.diff2.diff, by=c("Medium", "Drug", "Species", "UniqueID"))
  sample2sample_raw.df = d.f %>% dplyr::distinct(Sample, SampleLong)
  sample2sample_raw = sample2sample_raw.df$SampleLong
  names(sample2sample_raw) = sample2sample_raw.df$Sample

  #
  # Subset of data for PCA (All) and for PCA (MS/MS)
  #
  d1 = d %>%
    dplyr::mutate(subset="MS/MS validated only") %>%
    dplyr::filter(grepl("sacch", Species) & is_msms_validated)
  d2 = d %>% 
    # dplyr::filter(!is_peptide & is_mapped_metabolite & !is_lipid & is_main_fragment) %>% 
    dplyr::mutate(subset="mapped metabolites only") %>%  
    dplyr::group_by(Compound, Medium, Drug, Species) %>%
    dplyr::filter(any(Intensity.has4replicates)) %>%
    dplyr::group_by(Compound, Medium, Drug) %>%
    dplyr::filter(!any(is_significant=="Yes" & Species=="-"))

  library(IRanges)
  selected_metabolites = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp13lcms/selected_metabolites.tsv") %>%
    dplyr::select(CompoundName, MetlinID, mz)
  selected_metabolites_ranges = IRanges::IRanges(start=round(selected_metabolites$mz*1e6), end=round(selected_metabolites$mz*1e6))
  d_ranges = IRanges::IRanges(start=d$Mass*1e6, end=d$Mass*1e6)
  d3.index = as.data.frame(IRanges::distanceToNearest(d_ranges, selected_metabolites_ranges))
  d3 = d[d3.index$distance/1e6 < 5e-6,] %>%
    dplyr::left_join(selected_metabolites, by="MetlinID") %>%
    # dplyr::filter(!is_peptide & is_mapped_metabolite & !is_lipid & is_main_fragment) %>%
    dplyr::mutate(subset="predefined list") %>%
    dplyr::group_by(Compound, Medium, Drug, Species) %>%
    dplyr::filter(any(Intensity.has4replicates)) %>%
    dplyr::group_by(Compound, Medium, Drug) %>%
    dplyr::filter(!any(is_significant=="Yes" & Species=="-"))
  
  #
  # Export PCA data
  #
  d1.wide = d.wide %>% 
    dplyr::inner_join(d1 %>% dplyr::distinct(Species, Drug, UniqueID), by=c("Drug", "Species", "UniqueID")) %>%
    dplyr::select(UniqueID, MSMS_name, Significant=is_any_significant, dplyr::matches("Intensity|p-value|Effect size"))
  readr::write_tsv(d1.wide, path="reports/exp13lcms_150profileSsacch_msms.txt")

  d2.wide = d.wide %>% 
    dplyr::inner_join(d2 %>% dplyr::distinct(Species, Drug, UniqueID), by=c("Drug", "Species", "UniqueID")) %>% 
    dplyr::mutate(MappedToPeptide=ifelse(is_peptide, "Yes", "No"), MappedToMetabolite=ifelse(is_mapped_metabolite, "Yes", "No"), MappedToLipid=ifelse(is_lipid, "Yes", "No"), MappedToMainFragment=ifelse(is_main_fragment, "Yes", "No")) %>%
    dplyr::select(UniqueID, dplyr::matches("Intensity|p-value|Effect size"))
  d2.wide = rbind(c("", sample2sample_raw), d2.wide)
  readr::write_tsv(d2.wide, path="reports/exp13lcms_150profileSsacch_all_extra.txt")
  
  d3.wide = d.wide %>% 
    dplyr::inner_join(d3 %>% dplyr::distinct(Species, Drug, UniqueID), by=c("Drug", "Species", "UniqueID")) %>%
    dplyr::select(UniqueID, MSMS_name, Significant=is_any_significant, dplyr::matches("Intensity|p-value|Effect size"))
  readr::write_tsv(d3.wide, path="reports/exp13lcms_150selected.txt")
  
  d12 = dplyr::bind_rows(
    d1 %>% dplyr::distinct(Species, Drug, UniqueID), 
    d2 %>% dplyr::distinct(Species, Drug, UniqueID), 
    d3 %>% dplyr::distinct(Species, Drug, UniqueID)) %>% dplyr::distinct(Species, Drug, UniqueID)
  d12.wide = d.wide %>% 
    dplyr::inner_join(d12, by=c("Drug", "Species", "UniqueID")) %>%
    dplyr::select(UniqueID, MSMS_name, Significant=is_any_significant, dplyr::matches("Intensity|p-value|Effect size"))
  readr::write_tsv(d12.wide, path="reports/exp13lcms_150profileSsacch.txt", na="")
  
  #
  # Plot raw MS/MS values
  #
  d1  %>%
    dplyr::filter(ConcentrationFinal!=70) %>%
    dplyr::group_by(UniqueID) %>%
    dplyr::mutate(Compound=paste0(ifelse(is_any_significant=="Yes", "+", ""), ifelse(any(Intensity.has4replicates=="Yes"), "!", ""), ifelse(any(is_any_significant=="Yes" | Intensity.has4replicates=="Yes"), " ", ""), MSMS_name)) %>% #, "\n", Compound
    dplyr::ungroup() %>%
    bind_cols(group_col = floor(group_indices(., UniqueID) %/% 10)+1) %>% 
    ggplot() +
      geom_boxplot(aes(x=Compound, y=Intensity, fill=ConcentrationFinal, color=is_significant)) +
      scale_color_manual(values=c("No"="#333333", "Yes"="#FF0000")) +
      coord_flip() +
      facet_wrap(~group_col, scales="free", nrow=3)
  
  
  data_long.pca_species = dplyr::bind_rows(d1, d2, d3) %>%
    dplyr::group_by(subset, Medium, Drug, Species) %>%
    dplyr::do((function(z){
      zz <<- z
      z.wide = z %>%
        dplyr::mutate(Row=paste0(ConcentrationFinal, " (", Replicate, ")")) %>%
        reshape2::dcast(Plate + Medium + Drug + Species + ConcentrationFinal + Replicate + Row ~ UniqueID, value.var="Intensity") %>%
        tibble::column_to_rownames("Row")
      z.mat = as.matrix(z.wide %>% dplyr::select(-(Plate:Replicate)))
      z.pca = prcomp(z.mat)
      z.eigs = as.data.frame(t(z.pca$sdev^2/sum(z.pca$sdev^2)))
      colnames(z.eigs) = paste0("PC", 1:ncol(z.eigs), ".var")
      z.out = cbind(z.wide %>% dplyr::select(Plate:Replicate), z.pca$x[,1:2], z.eigs[1:2])
      
      z.out
    })(.)) %>%
    dplyr::mutate(subset_variance=paste0(subset, "(", round(PC1.var*100), "%, ", round(PC2.var*100), "%)")) %>%
    dplyr::mutate(SHAPE=pca_species_shapes[Species], PCH=pca_species_pch2[Species])
  
  cairo_pdf("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/reports/exp13lcms_150pcaSsacch.pdf", width=30, height=10, onefile=T)
  for(s in unique(data_long.pca_species$subset)) {
    data_long.pca_species.f = data_long.pca_species %>% dplyr::filter(grepl("sacch", Species) & Drug=="Duloxetine" & subset==s) #  & subset=="mapped metabolites only"
    p = ggplot(data_long.pca_species.f) +
      geom_text(aes(PC2, PC1, color=ConcentrationFinal, label=PCH, alpha=ConcentrationFinal), size=4, family="Lucida Sans Unicode", show.legend=F) +
      scale_color_manual(values=c("0"="#999999", "10"="#B27F72", "20"="#D85839", "30"="#FF3200", "40"="#E62710", "50"="#CD1C20", "70"="#780A0D")) +
      scale_alpha_manual(values=c("0"=0.7, "10"=0.7, "20"=0.7, "30"=0.7, "40"=0.7, "50"=0.9, "70"=0.9)) +
      coord_fixed(ratio=with(data_long.pca_species.f, PC1.var/PC2.var)[1]) +
      labs(
        color="Concentration",
        title=paste0("C. saccharolyticum (", s, ")"),
        fill="Concentration", color="Concentration", alpha="Concentration", 
        x=paste0("PC2 (", round(data_long.pca_species.f$PC2.var[1]*100, 1),"%)"), 
        y=paste0("PC1 (", round(data_long.pca_species.f$PC1.var[1]*100, 1),"%)")) +
      theme_classic(base_size=18) +
      theme(strip.text.y=element_text(angle = 180))
    print(p)
  }

  dev.off()
}

exp13lcms.150heatmap = function()
{
  
  #
  # Plots hits statistics
  #
  pdf("reports/exp13lcms_msms_validated.pdf", width=24, height=15)
  data_long.significant_all = data_long.diff2 %>%
    dplyr::filter(total_significant>=3 & total_significant_none==0)
  data_long.significant_all_sacch = data_long.significant_all %>%
    dplyr::filter(total_significant_sacch>=3)
  data_long.significant = data_long.significant_all %>%
    # dplyr::filter(!is_peptide & is_mapped_metabolite & !is_lipid & is_main_fragment) %>%
    dplyr::filter(!is_peptide & is_mapped_metabolite & !is_lipid) #  & is_main_fragment
  data_long.significant_sacch = data_long.significant %>%
    dplyr::filter(total_significant_sacch>=3) #  & is_msms is_saccharolyticum &
  data_long.msms = data_long.diff2 %>% dplyr::filter(is_msms_validated)
  
  
  #
  # Hist overview
  #
  cap = 300
  data_long.significant %>% 
    dplyr::mutate(Species=paste0(Species, " (", total_concentrations, " concentrations)")) %>%
    dplyr::group_by(Species, Drug) %>% 
    dplyr::summarise(
      count2=length(unique(Compound[total_significant>=2])),
      count3=length(unique(Compound[total_significant>=3])), 
      count4=length(unique(Compound[total_significant>=4])), 
      count5=length(unique(Compound[total_significant>=5])), 
      count2_cap=pmin(count2, cap), count3_cap=pmin(count3, cap), count4_cap=pmin(count4, cap), count5_cap=pmin(count5, cap), 
      total_concentrations=unique(total_concentrations)) %>%
    ggplot() +
    geom_bar(aes(x=Species, y=count2, fill=Drug), stat="identity", position=position_dodge2(preserve="single"), alpha=0.3) +
    geom_bar(aes(x=Species, y=count3, fill=Drug), stat="identity", position=position_dodge2(preserve="single"), alpha=0.3) +
    geom_bar(aes(x=Species, y=count4, fill=Drug), stat="identity", position=position_dodge2(preserve="single"), alpha=0.3) +
    geom_bar(aes(x=Species, y=count5, fill=Drug), stat="identity", position=position_dodge2(preserve="single")) +
    geom_text(aes(x=Species, y=count2_cap*0.95, color=Drug, label=paste0(count2, "\n", count3, "\n", count4, "\n", count5)), stat="identity",  color="black", position=position_dodge2(preserve="single", width=1), hjust=1) +
    labs(y="Number of significant compounds (in 2/3/4/5 concentrations)", x="") +
    coord_flip(ylim=c(0, cap))
  
  
  #
  # Export hits to MetaboLight
  #
  hits.kegg = unique(na.omit(data_long.significant_sacch$KEGG_ID))
  unichem.chembl2hmdb = unichem.mapping(1, 18) %>% dplyr::rename(hmdb="dest", chembl.hmdb="src") %>% dplyr::mutate(hmdb=gsub("HMDB0*", "HMDB", hmdb))
  unichem.chembl2kegg = unichem.mapping(1, 6) %>% dplyr::rename(kegg="dest", chembl.kegg="src")
  unichem.chembl2drugbank = unichem.mapping(1, 2) %>% dplyr::rename(drugbank="dest", chembl.drugbank="src")
  hmdb_df = readr::read_tsv("data/db/hmdb_metabolites.tsv", na="") %>% dplyr::mutate(HMDB_ID=gsub("HMDB0*", "HMDB", HMDB_ID))
  drugbank_df = readr::read_csv("data/db/drugbank/structure links.csv", na="") %>% setNames(., paste0(names(.), ".drugbank"))
  
  data_long.export_all_metabolight = data_long.significant_all_sacch %>%
    dplyr::mutate(MetlinID=as.character(MetlinID)) %>% 
    dplyr::distinct(Compound, .keep_all=T) %>%
    dplyr::left_join(unichem.chembl2hmdb, by=c("HMP_ID"="hmdb")) %>%
    dplyr::left_join(unichem.chembl2kegg, by=c("KEGG_ID"="kegg")) %>%
    dplyr::left_join(unichem.chembl2drugbank, by=c("DrugBankID"="drugbank")) %>%
    dplyr::left_join(hmdb_df, by=c("HMP_ID"="HMDB_ID")) %>%
    dplyr::left_join(drugbank_df, by=c("DrugBankID"="DrugBank ID.drugbank")) %>%
    dplyr::mutate(
      database_identifier=dplyr::case_when(
        !is.na(chembl.drugbank)~chembl.drugbank, 
        !is.na(chembl.kegg)~chembl.kegg,
        !is.na(chembl.hmdb)~chembl.hmdb,
        grepl("Hippuric acid", MSMS_name) ~ "CHEMBL461",
        T~NA_character_))
  chembl_df = chembl.molecules(unique(na.omit(data_long.export_all_metabolight$database_identifier)))
  metaboanalyst_metlin = metaboanalyst.mapper(unique(na.omit(data_long.export_all_metabolight$MetlinID)), "metlin") %>% setNames(paste0(names(.), ".metaboanalyst"))
  data_long.export_all_metabolight = data_long.export_all_metabolight %>%
    dplyr::left_join(metaboanalyst_metlin, by=c("MetlinID"="metlin_id.metaboanalyst")) %>%
    dplyr::left_join(chembl_df, by=c("database_identifier"="chembl_id")) %>%
    dplyr::mutate(
      smiles=dplyr::case_when(
        !is.na(chembl_smiles)~chembl_smiles,
        !is.na(HMDB_SMILES)~HMDB_SMILES, 
        !is.na(smiles.metaboanalyst)~smiles.metaboanalyst,
        T~NA_character_)) %>%
    dplyr::mutate(
      inchi=dplyr::case_when(
        !is.na(chembl_inchi)~chembl_inchi,
        !is.na(HMDB_INCHI)~HMDB_INCHI,
        T~NA_character_))  %>%    
    dplyr::mutate(
      chemical_formula=gsub(" ", "", Formula),
      metabolite_identification=CompoundName, #	Compound/Annotation name 
      mass_to_charge=Mass, #	Mass of the compound
      fragmentation=NA_character_, #	we dont need this
      modifications=IonSpecies, #	whether it is [M+H]+ or [M+H]- charge 	its generally 1
      retention_time=RetentionTime,
      taxid="610130", #	microbe classification , Please see attached excel sheet for eg.
      species=Species, #	Microbe species name
      database="metlin", #	what is this? Metlin metabolite or lipid
      database_version="", #	metlin database version name
      reliability="", #	i dont know, maybe we can skip it
      uri="", #	what is this?
      search_engine="metlin", #	metlin database
      search_engine_score=Score, #	database score
      smallmolecule_abundance_sub="", #	I dont know, we can skip it
      smallmolecule_abundance_stdev_sub="", #	I dont know, we can skip it
      smallmolecule_abundance_std_error_sub="" 	# I dont know, we can skip it)
    ) %>% dplyr::distinct(database_identifier, chemical_formula, smiles, inchi, metabolite_identification, mass_to_charge, fragmentation, modifications, retention_time, taxid, species, database, database_version, reliability, uri, search_engine, search_engine_score, smallmolecule_abundance_sub, smallmolecule_abundance_stdev_sub, smallmolecule_abundance_std_error_sub)
  readr::write_tsv(data_long.export_all_metabolight, path="reports/exp13lcms_candidates_lcms_all_metabolight.tsv  ", col_names=T, na="")
  
  
  #
  # Heatmap
  #
  data_long.heatmap = data_long.significant_sacch %>%
    reshape2::dcast(Compound ~ Species + ConcentrationFinal, value.var="diff") %>%
    replace(is.na(.), 0) %>%
    tibble::column_to_rownames("Compound")
  data_long.heatmap_sign = data_long.significant_sacch %>%
    reshape2::dcast(Compound ~ Species + ConcentrationFinal, value.var="is_significant") %>%
    replace(is.na(.), F) %>%
    tibble::column_to_rownames("Compound")
  data_long.heatmap_adducts = data_long.significant_sacch %>% 
    tidyr::separate_rows(IonSpecies, sep=" ") %>%
    dplyr::distinct(Compound, IonSpecies) %>%
    reshape2::dcast(Compound ~ IonSpecies, value.var="IonSpecies", fun.aggregate=function(x) ifelse(length(x)>0, "Yes", "No")) %>%
    tibble::column_to_rownames("Compound") %>%
    as.matrix()
  data_long.heatmap_mass = data_long.significant_sacch %>%
    dplyr::mutate(RetentionTime=round(RetentionTime, 3)) %>%
    dplyr::mutate(
      Significant=total_significant_sacch,
      IsValidated=ifelse(is_msms_validated, "Yes", "No"),
      MSMS_name=ifelse(is.na(MSMS_name), "", MSMS_name),
      InDrugbank=ifelse(!is.na(DrugBankID), "Yes", "No"), 
      IsSaccharolyticum=ifelse(is_saccharolyticum, "Yes", "No"),
      KEGG_METNAME=ifelse(is.na(is_saccharolyticum.metabolite), "", is_saccharolyticum.metabolite),
      MainFragment=ifelse(is_main_fragment, "Yes", "No"), 
      MSMS=ifelse(is_msms, "Yes", "No"),
      Formula=gsub(" ", "", Formula))  %>%
    dplyr::distinct(Compound, .keep_all=T)%>%
    tibble::column_to_rownames("Compound")
  data_long.heatmap_abs_summary = data_long %>%
    dplyr::filter(ConcentrationFinal <= 50) %>%
    dplyr::inner_join(data_long.significant_sacch %>% dplyr::distinct(Compound, Medium, Drug, Species), by=c("Compound", "Medium", "Drug", "Species")) %>%
    dplyr::group_by(Species, Medium, Drug, Compound, ConcentrationFinal) %>%
    dplyr::summarise(IntensityN=sum(Intensity>0, na.rm=T), Intensity=ifelse(any(Intensity>0), median(Intensity[Intensity>0], na.rm=T), 0)) %>%
    dplyr::group_by(Medium, Drug, Compound) %>%
    dplyr::summarise(
      IntensityN.0=ifelse(any(grepl("sacch", Species) & ConcentrationFinal==0), IntensityN[grepl("sacch", Species) & ConcentrationFinal==0], 0),
      IntensityN.max=ifelse(any(grepl("sacch", Species) & ConcentrationFinal!=0), max(IntensityN[grepl("sacch", Species) & ConcentrationFinal!=0]), 0),
      Intensity.0=ifelse(any(grepl("sacch", Species) & ConcentrationFinal==0), Intensity[grepl("sacch", Species) & ConcentrationFinal==0], 0), 
      Intensity.max=ifelse(any(ConcentrationFinal!=0), Intensity[grepl("sacch", Species)][which.max(abs(Intensity[grepl("sacch", Species)]-Intensity[grepl("sacch", Species) & ConcentrationFinal==0]))], 0)) %>%
    dplyr::ungroup() 
  data_long.heatmap_absN = data_long.heatmap_abs_summary %>%
    dplyr::distinct(Compound, IntensityN.0, IntensityN.max) %>%
    dplyr::mutate(N=paste0(IntensityN.0, " / ", IntensityN.max)) %>%
    tibble::column_to_rownames("Compound") 
  data_long.heatmap_abs = data_long.heatmap_abs_summary %>%
    dplyr::distinct(Compound, Intensity.0, Intensity.max) %>%
    tibble::column_to_rownames("Compound") %>%
    as.matrix()
  heatmap.empty = data_long.heatmap_absN$N;  names(heatmap.empty) = rownames(data_long.heatmap_absN)
  
  
  colors.yesno = c("Yes"="#1F78B4FF", "No"="#A6CEE3")
  column_labels = paste0(paste(rep(" ", 20), collapse=""), gsub(".*_", "", colnames(data_long.heatmap)))
  row_labels = rownames(data_long.heatmap) #paste(gsub(" Es.-.*", "", rownames(data_long.heatmap))) 
  ComplexHeatmap::Heatmap(data_long.heatmap, 
                          cluster_columns=F, 
                          column_labels=column_labels, 
                          row_labels=row_labels, 
                          column_split=gsub("^([A-Z])[^ ]+", "\\1.", gsub("_.*", "", colnames(data_long.heatmap))), 
                          col=circlize::colorRamp2(seq(-10, 10, length = 11), rev(RColorBrewer::brewer.pal(11, "Spectral"))),
                          row_names_max_width = max_text_width(rownames(data_long.heatmap)),
                          # cluster_rows=data_long.heatmap_rowdend,
                          right_annotation=ComplexHeatmap::rowAnnotation(
                            Significant=ComplexHeatmap::anno_simple(data_long.heatmap_mass[,"Significant", drop=F], gp=grid::gpar(border=T, col="#000000")),
                            Adduct=data_long.heatmap_adducts, 
                            # Score=ComplexHeatmap::anno_barplot(data_long.heatmap_mass[,"Score", drop=F]),
                            IsValidated=ComplexHeatmap::anno_simple(data_long.heatmap_mass[,"IsValidated", drop=F], col=colors.yesno, gp=grid::gpar(border=T, col="#000000")),
                            MainFragment=ComplexHeatmap::anno_simple(data_long.heatmap_mass[,"MainFragment", drop=F], col=colors.yesno, gp=grid::gpar(border=T, col="#000000")),
                            Saccharolyticum=ComplexHeatmap::anno_simple(data_long.heatmap_mass[,"IsSaccharolyticum", drop=F], col=colors.yesno, gp=grid::gpar(border=T, col="#000000")),
                            IntensityNonZero=ComplexHeatmap::anno_text(heatmap.empty, width=max_text_width(heatmap.empty)),
                            Intensity=ComplexHeatmap::anno_points(cbind(data_long.heatmap_abs, intensity_min), gp=grid::gpar(col=rep(c(1,2,3), nrow(data_long.heatmap_abs)))), 
                            Mass=ComplexHeatmap::anno_text(data_long.heatmap_mass$Mass, width=max_text_width(data_long.heatmap_mass$Mass)),  
                            Retention=ComplexHeatmap::anno_text(data_long.heatmap_mass[,"RetentionTime"], width=max_text_width(data_long.heatmap_mass$RetentionTime*1.5)),  
                            Formula=ComplexHeatmap::anno_text(data_long.heatmap_mass[,"Formula"], width=max_text_width(data_long.heatmap_mass$Formula)),
                            # KEGG_METNAME=ComplexHeatmap::anno_text(data_long.heatmap_mass[,"KEGG_METNAME"], width=max_text_width(data_long.heatmap_mass$KEGG_METNAME)),
                            MSMS_name=ComplexHeatmap::anno_text(data_long.heatmap_mass[,"MSMS_name"], width=max_text_width(data_long.heatmap_mass$MSMS_name)),
                            gp=grid::gpar(border=T, col="#FFFFFF"),
                            col=list(Adduct=colors.yesno),
                            gap=unit(5, "points")
                          ),
                          rect_gp=grid::gpar(border=T, col="#FFFFFF"),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.rect(x=x, y=y, width=width, height=height, gp=gpar(col=ifelse(data_long.heatmap_sign[i,j], "#FFD700", "#FFFFFF00"), fill=fill))
                          })
  
  
  pdf("reports/exp12lcms_venn.pdf", width=10, height=5)
  data_long.f.pcavars = data_long.f %>%
    dplyr::group_by(Medium, Drug, Species) %>%
    dplyr::do((function(z){
      zz <<- z
      z.wide = z %>%
        reshape2::dcast(Plate + Medium + Drug + Species + ConcentrationFinal + Replicate ~ Compound, value.var="Intensity")
      z.pca = prcomp(z.wide %>% dplyr::select(-(Plate:Replicate)))
      z.pca_pov = z.pca$sdev^2/sum(z.pca$sdev^2)
      z.pca_pov.df = as.data.frame(z.pca_pov) %>% tibble::rownames_to_column("PC") %>% dplyr::mutate(PC=paste0("PC",1:n())) %>% dplyr::rename(ExplainedVariance="z.pca_pov")
      z.pca_vars = z.pca$rotation^2/rowSums(z.pca$rotation^2)
      z.pca_vars.rel = t(t(z.pca_vars)*z.pca_pov)
      
      z.out = z.pca_vars %>% 
        reshape2::melt(varnames=c("Compound", "PC"), value.name="Contribution") %>%
        dplyr::inner_join(z.pca_pov.df, by="PC") %>%
        dplyr::arrange(dplyr::desc(Contribution)) %>%
        dplyr::group_by(PC) %>%
        dplyr::mutate(ContributionRank=1:n(), ContributionRel=ExplainedVariance*Contribution) %>%
        dplyr::ungroup()
      z.out = cbind(z.out, z %>% dplyr::distinct(Plate, Medium, Drug, Species))
      
      z.out
    })(.)) %>%
    dplyr::ungroup()
  
  data_long.f.pcavars.f = data_long.f.pcavars %>%
    dplyr::group_by(Compound, Drug, Species) %>%
    dplyr::arrange(dplyr::desc(ContributionRel)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() 
  data_long.f.pcavars.f_sum = data_long.f.pcavars.f %>%
    dplyr::group_by(Compound, Medium, Drug, Species) %>%
    dplyr::summarize(ContributionRelMax=max(ContributionRel)) %>%
    dplyr::group_by(Medium, Drug, Species) %>%
    dplyr::arrange(dplyr::desc(ContributionRelMax)) %>%
    dplyr::mutate(CompoundRank=1:n()) %>%
    dplyr::filter(CompoundRank<=100) %>%
    dplyr::ungroup()
  data_long.f.pcavars.f.ggplot = data_long.f.pcavars.f %>%
    dplyr::mutate(PC=as.character(PC)) %>% 
    dplyr::filter(ContributionRel>0.01) %>%
    dplyr::inner_join(data_long.f.pcavars.f_sum %>% dplyr::distinct(Compound), by="Compound") %>%
    dplyr::inner_join(annotations_ions, by="Compound")
  x = data_long %>%
    dplyr::inner_join(data_long.f.pcavars.f_sum %>% dplyr::distinct(Compound), by="Compound") %>%
    dplyr::anti_join(data_long.f0 %>% dplyr::distinct(Compound), by="Compound") %>%
    dplyr::inner_join(annotations_ions, by="Compound") %>%
    dplyr::mutate(rowname=paste0(Drug, " ", ConcentrationFinal, "uM ", Species, " (", Plate, " ", Replicate, ")"), colname=paste0(Mass, "    ", Compound, " ", `Ion Species`)) %>%
    reshape2::dcast(colname~rowname, value.var="Intensity") %>%
    tibble::column_to_rownames("colname")
  pheatmap::pheatmap(x, cluster_cols=F, gaps_col=seq(4, ncol(x), 4), fontsize_col=6)
  
  
  ggplot(data_long.f.pcavars.f.ggplot) +
    geom_bar(aes(x=paste0(Compound, `Ion Species`, " (", Mass, ")"), y=ContributionRel, fill=PC), stat="identity", position="dodge") +
    coord_flip() +
    facet_wrap(Species~Drug, nrow=1)
}
