dir.create("data/exp2metabolomics/metabolights", showWarnings=F)
dir.create("data/exp1depletion/metabolights", showWarnings=F)
library(dplyr)
library(readr)
library(ggplot2)
source("functions.R")

exp12screen.export_metabolights = function()
{
  #
  # Load Species and drug data
  #
  drug_map = readr::read_delim("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/drug_map.tsv", "\t")
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")

  bioaccumulation_processed_map = readr::read_delim("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp2metabolomics/data.depletionmodeassay_long_new.csv", ",") %>%
    # dplyr::filter(!is.na(Drugs)) %>%
    dplyr::distinct(Extraction, Sample.Set.Name, SampleName, Channel.Description, Vial, RetentionTime=RT, drug.short2=Drugs, Replicate, Ctrl, .keep_all=T) %>% 
    dplyr::left_join(drug_map %>% dplyr::select(drug.short2, drug.long, drug.chebi, drug.formula, drug.smiles, drug.inchi, drug.uplc_uv), by="drug.short2") %>% 
    dplyr::left_join(bug_map %>% dplyr::select(species.code, species.short, species.long, species.taxid), by=c("Bugs"="species.code")) %>% 
    dplyr::mutate(Status=dplyr::case_when(Ctrl=="ctrl"~"ctrl", Ctrl=="smpl"~"sample", Ctrl=="zero"~"GMM")) %>%
    dplyr::distinct(Extraction, Sample.Set.Name, SampleName, Channel.Description, Status, Replicate, Vial, RetentionTime, drug.long, drug.chebi, drug.formula, drug.smiles, drug.inchi, species.short, species.long, species.taxid, drug.uplc_uv, .keep_all=T) %>% 
    group_by(SampleName, Sample.Set.Name, Channel.Description, Vial, Replicate) %>%
    dplyr::mutate(ChannelNo=1:n()) %>%
    dplyr::ungroup()

  bioaccumulation_plate_map = readr::read_tsv("data/exp2metabolomics/raw_annotations/uplc_interactions_mapping.tsv") %>% 
    dplyr::mutate(GroupName=gsub("^[0-9]+_|_martina", "", Acq.Method.Set)) %>%
    group_by(AssayName, SampleName, Sample.Set.Name, Channel.Description, Vial) %>%
    dplyr::mutate(ChannelNo=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(AssayName, SampleName, Sample.Set.Name, Plate=Plate_Martina, GroupName, Channel.Description, Vial, ChannelNo, Channel.Id) 
  
  
  bioaccumulation_export_samples = bioaccumulation_plate_map %>% 
    dplyr::inner_join(bioaccumulation_processed_map, by=c("SampleName", "Sample.Set.Name", "Channel.Description", "Vial", "ChannelNo")) %>%
    dplyr::group_by(Extraction, Sample.Set.Name, Bugs, GroupName, Status) %>%
    dplyr::mutate(WellReplicate=match(Vial, unique(Vial))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-ChannelNo) %>%
    dplyr::filter(Bugs!="med") %>%
    dplyr::select(AssayName, Extraction, Sample.Set.Name, Plate, SampleName, GroupName, RetentionTime, Status, WellReplicate, Replicate, Vial, drug.long, drug.chebi, drug.formula, drug.smiles, drug.inchi, drug.uplc_uv, species.short, species.long, species.taxid, Channel.Description, Channel.Id, Area, Area.caffee)  %>%
    dplyr::filter(!grepl("^med_dot", SampleName)) %>%
    dplyr::distinct(AssayName, Plate, Extraction, GroupName, Status, WellReplicate, .keep_all=T) %>%
    dplyr::group_by(GroupName) %>%
    dplyr::mutate(GroupName.long=na.omit(drug.long)[1], GroupName.chebi=na.omit(drug.chebi)[1]) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Plate) %>%
    dplyr::mutate(Plate.long=na.omit(species.long)[1], Plate.taxid=na.omit(species.taxid)[1]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SampleName.long=paste0("Bioactivity-", Plate, "-", Extraction, "-", GroupName, "-", Status, "-", WellReplicate)) %>%
    dplyr::mutate(report=paste0(SampleName.long,"-channel", Channel.Id, ".arw"))
  
  
  #
  # Collect all the raw files in one folder
  #
  bioaccumulation_input_dir.raw = "/g/patil/Share/EmpowerProjects/Martina"
  bioaccumulation_input_dir.arw = "data/exp2metabolomics/arw"
  bioaccumulation_output_dir.arw = "data/exp2metabolomics/metabolights/arw_filtered/Bioactivity"
  bioaccumulation_output_dir.raw = "data/exp2metabolomics/metabolights/raw_filtered/Bioactivity"
  bioaccumulation_output.collect_script = "data/exp2metabolomics/metabolights/filter_files.sh"
  
  writeLines(sprintf("Duplicated channels: %i", sum(duplicated(unique(bioaccumulation_export_samples$Channel.Id)))))
  
  # Check missing report files
  bioaccumulation_input_dir.report_files = list.files(bioaccumulation_input_dir.arw)
  bioaccumulation_input_dir.report_files.expected = sprintf("report_raw%s.arw", bioaccumulation_export_samples$Channel.Id)
  writeLines(sprintf("Missing report files: %i", length(setdiff(bioaccumulation_input_dir.report_files.expected, bioaccumulation_input_dir.report_files))))
  
  if(F) {
    # Collect all reports into single file
    
    i.run = 0
    reports = data.frame()
    reports.100 = data.frame()
    for(i in unique(bioaccumulation_export_samples$Channel.Id)) {
      if(i %in% reports$Channel.Id) next
      i.file = sprintf(sprintf("%s/report_raw%i.arw", bioaccumulation_input_dir.arw, i))
      if(file.exists(i.file)) {
        i.run = i.run+1
        report.i = readr::read_tsv(i.file, skip=2, col_names=c("RetentionTime", "Intensity"), col_types=readr::cols()) %>% dplyr::mutate(Channel.Id=i)
        reports.100 = rbind(reports.100, report.i)
        if(i.run %% 100==0) {
          print(sprintf("%i/%i (added: %i)", which(i==unique(bioaccumulation_export_samples$Channel.Id)), length(unique(bioaccumulation_export_samples$Channel.Id)), i.run))
          reports = rbind(reports, reports.100)
          reports.100 = data.frame()
        }
      } else {
        writeLines(sprintf("File doesn't exist: %s", i.file))
      }
    }
    reports = rbind(reports, reports.100)
    reports_distinct = reports %>% dplyr::distinct(Channel.Id, RetentionTime, .keep_all=T) %>% dplyr::select(Channel.Id, ElutionTime=RetentionTime, Intensity)
    reports.wide = reports_distinct %>%
      reshape2::dcast(Channel.Id ~ RetentionTime, value.var="Intensity")
    readr::write_tsv(reports_distinct, path="data/exp2metabolomics/uplc_raw.tsv", na="")
  }
  
  
  cmd = list()
  if(!file.exists(bioaccumulation_output_dir.arw)) cmd[["create_arw"]] = sprintf("mkdir %s", bioaccumulation_output_dir.arw)
  cmd[["copy_reports"]] = sprintf("cp %s/%s %s/%s", bioaccumulation_input_dir.arw, bioaccumulation_input_dir.report_files.expected, bioaccumulation_output_dir.arw, bioaccumulation_export_samples$report)
  # cmd[["copy_eic"]] = sprintf("data/exp2metabolomics/data.eic.csv %s/raw.zip", bioaccumulation_output_dir)
  writeLines(unlist(cmd), con=bioaccumulation_output.collect_script)
  
  
  
  ######################################################################################################################################
  depletion_plate_map = readr::read_tsv("data/exp2metabolomics/raw_annotations/uplc_plata_mapping.tsv")
  depletion_channel_map = readr::read_tsv("data/exp2metabolomics/raw_annotations/uplc_interactions_mapping.tsv")  %>%
    dplyr::filter(Plate_Martina!="EcoliiAi1A2_one") %>%
    dplyr::mutate(Plate_Martina=ifelse(Plate_Martina == "BlonguminfantisA_one", "BlonguminfantisA_EcoliiAi1A2", Plate_Martina)) # Data for this is wrong
  
  #
  # Load mapping between raw files and final data
  #
  growth_file2plate = readr::read_csv("data/exp0growth/160106_FileToPlate.csv")
  growth_curves_annotations = readr::read_tsv("data/exp0growth/raw_annotations/curves_annotation_2016-11-28.tab") %>%
    dplyr::inner_join(growth_file2plate, by="File")
  
  
  #
  # Prepare annotations 
  #
  # Replicate - Technical (measurement replicate) of the same well
  # Well - Replicate
  # Key (Plate, GroupName, Status, Well, Replicate) - Unique key
  load("data/exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
  depletion_processed_map = data.clean %>% 
    dplyr::filter(!grepl("Mix", Plate)) %>%
    dplyr::mutate(WellReplicate=ifelse(Status=="sample", Well, 1)) %>%
    dplyr::mutate(MeasurementReplicate=ifelse(Status!="sample", Well, 1)) %>%
    dplyr::rename(species.long="Species.x") %>% 
    dplyr::left_join(bug_map %>% dplyr::select(species.long, species.taxid), by=c("species.long"))
  
  #
  # Prepare mappings
  #
  depletion_channel_map_final = depletion_channel_map %>% 
    dplyr::filter(AssayName=="Screen") %>% 
    dplyr::filter(!grepl("Spect", Channel.Description)) %>%
    # dplyr::filter(!grepl("down|purging", Acq.Method.Set)) %>%
    dplyr::group_by(AssayName, Sample.Set.Name, SampleName, Vial) %>%
    dplyr::mutate(ChannelNo=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Plate_Martina_Number=as.numeric(gsub(":.*", "", Vial)), Well=gsub(".*:", "", Vial)) %>% 
    dplyr::inner_join(depletion_plate_map, by=c("Plate_Martina", "Plate_Martina_Number")) %>%
    dplyr::mutate(drug.something=gsub("^[^_]+_", "", Acq.Method.Set)) %>%
    dplyr::mutate(drug.something=dplyr::case_when(drug.something=="eyet"~"ezet", drug.something=="rosuva"~"rosu", drug.something=="donez"~"done", T~drug.something)) %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short, drug.long1=drug.long, drug.chebi1=drug.chebi, drug.formula1=drug.formula, drug.smiles1=drug.smiles, drug.inchi1=drug.inchi, drug.uplc_uv1=drug.uplc_uv, drug.uplc_retention1=drug.uplc_retention), by=c("drug.something"="drug.short")) %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short2, drug.long2=drug.long, drug.chebi2=drug.chebi, drug.formula2=drug.formula, drug.smiles2=drug.smiles, drug.inchi2=drug.inchi, drug.uplc_uv2=drug.uplc_uv, drug.uplc_retention2=drug.uplc_retention), by=c("drug.something"="drug.short2")) %>%
    dplyr::mutate(drug.long=dplyr::case_when(!is.na(drug.long2)~drug.long2, !is.na(drug.long1)~drug.long1, T~NA_character_)) %>%
    dplyr::mutate(drug.chebi=dplyr::case_when(!is.na(drug.chebi2)~as.character(drug.chebi2), !is.na(drug.chebi1)~as.character(drug.chebi1), T~NA_character_)) %>%
    dplyr::mutate(drug.formula=dplyr::case_when(!is.na(drug.formula2)~as.character(drug.formula2), !is.na(drug.formula1)~as.character(drug.formula1), T~NA_character_)) %>%
    dplyr::mutate(drug.smiles=dplyr::case_when(!is.na(drug.smiles2)~as.character(drug.smiles2), !is.na(drug.smiles1)~as.character(drug.smiles1), T~NA_character_)) %>%
    dplyr::mutate(drug.inchi=dplyr::case_when(!is.na(drug.inchi2)~as.character(drug.inchi2), !is.na(drug.inchi1)~as.character(drug.inchi1), T~NA_character_)) %>%
    dplyr::mutate(drug.uplc_uv=dplyr::case_when(!is.na(drug.uplc_uv2)~as.character(drug.uplc_uv2), !is.na(drug.uplc_uv1)~as.character(drug.uplc_uv1), T~NA_character_)) %>%
    dplyr::mutate(drug.uplc_retention=dplyr::case_when(!is.na(drug.uplc_retention2)~as.character(drug.uplc_retention2), !is.na(drug.uplc_retention1)~as.character(drug.uplc_retention1), T~NA_character_)) %>%
    dplyr::select(-drug.long1, -drug.long2) %>%
    dplyr::mutate(Status=dplyr::case_when(
      grepl("ctrl", SampleName)~"ctrl",
      grepl("GMM", SampleName)~"GMM",
      grepl("down|purging", Acq.Method.Set)~"technical",
      T~"sample")) %>%
    dplyr::mutate(GroupName=ifelse(Status=="technical", NA_character_, substr(Acq.Method.Set, 1, 4))) %>%
    dplyr::mutate(drug.long=ifelse(Status=="GMM", NA_character_, drug.long)) %>%
    dplyr::group_by(Plate, GroupName, Status) %>%
    dplyr::mutate(WellReplicate=match(Vial, unique(Vial))) %>%
    dplyr::group_by(Plate, GroupName, Status, Vial) %>%
    dplyr::mutate(MeasurementReplicate=match(SampleName, unique(SampleName))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Extraction="Supernatant") %>%
    dplyr::select(AssayName, Sample.Set.Name, Extraction, Plate_Martina, SampleName, GroupName, Plate, Vial, WellReplicate, MeasurementReplicate, Status, drug.long, drug.chebi, drug.formula, drug.smiles, drug.inchi, drug.uplc_uv, drug.uplc_retention, Channel.Description, Channel.Id, ChannelNo) 
  
  # Pyri is excluded from drugs
  depletion_export_samples = depletion_processed_map %>%
    dplyr::left_join(depletion_channel_map_final %>% dplyr::mutate(has_data=T), by=c("Plate", "GroupName", "Status", "MeasurementReplicate", "WellReplicate")) %>%
    dplyr::select(AssayName, Extraction, Sample.Set.Name, Plate, SampleName, GroupName, Status, WellReplicate, MeasurementReplicate, Vial, drug.long, drug.formula, drug.chebi, drug.smiles, drug.inchi, drug.uplc_uv, drug.uplc_retention, species.long, species.taxid, Channel.Description, Channel.Id, DrugRatio, has_data) %>%
    dplyr::group_by(GroupName) %>%
    dplyr::mutate(GroupName.long=na.omit(drug.long)[1], GroupName.chebi=na.omit(drug.chebi)[1]) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Plate) %>%
    dplyr::mutate(Plate.long=na.omit(species.long)[1], Plate.taxid=na.omit(species.taxid)[1]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SampleName.long=paste0("Depletion-", Plate, "-", Extraction, "-", GroupName, "-", Status, "-", WellReplicate)) %>%
    dplyr::mutate(report=paste0(SampleName.long,"-channel", Channel.Id, ".arw"))
  stopifnot(any(!is.na(depletion_export_samples$has_data)))
  
  depletion_input_dir.arw = "data/exp1depletion/arw"
  depletion_input_dir.raw = "data/exp1depletion/raw" 
  depletion_output_dir.raw = "data/exp2metabolomics/metabolights/raw_filtered/Depletion"
  depletion_output_dir.arw = "data/exp2metabolomics/metabolights/arw_filtered/Depletion"
  depletion_output.collect_script = "data/exp1depletion/metabolights/filter_files.sh"
  
  writeLines(sprintf("Duplicated channels: %i", sum(duplicated(unique(depletion_export_samples$Channel.Id)))))
  
  # Check missing report files
  depletion_input_dir.report_files = list.files(depletion_input_dir.arw)
  depletion_input_dir.report_files.expected = sprintf("d%s.arw", depletion_export_samples$Channel.Id)
  writeLines(sprintf("Missing report files: %i", length(setdiff(depletion_input_dir.report_files.expected, depletion_input_dir.report_files))))
  
  #
  # Collect all reports into single file
  #
  if(F) {
    i.run = 0
    reports = data.frame()
    reports.100 = data.frame()
    for(i in unique(depletion_export_samples$Channel.Id)) {
      if(i %in% unique(reports$Channel.Id)) next
      i.file = sprintf(sprintf("%s/d%i.arw", depletion_input_dir.arw, i))
      if(file.exists(i.file)) {
        i.run = i.run+1
        print(sprintf("%i/%i (added: %i)", which(i==unique(depletion_export_samples$Channel.Id)), length(unique(depletion_export_samples$Channel.Id)), i.run))
        report.i = readr::read_tsv(i.file, skip=2, col_names=c("RetentionTime", "Intensity"), col_types=readr::cols())
        report.i_meta = readr::read_tsv(i.file, n_max=1, locale=readr::locale(asciify=T, encoding="Windows-1252"), col_types=readr::cols()) %>%
          data.frame() %>%
          dplyr::select(Channel.Description, Channel.Id, Sample.Set.Name, Sample.Set.Method, SampleName, Instrument.Method.Name, Vial)
        reports.100 = rbind(reports.100, cbind(report.i_meta, report.i))
        if(i.run %% 100==0) {
          reports = rbind(reports, reports.100)
          reports.100 = data.frame()
        }
      }
    }
    reports = rbind(reports, reports.100)
    reports_distinct = reports %>% dplyr::distinct(Channel.Id, RetentionTime, .keep_all=T) %>% dplyr::select(Channel.Id, ElutionTime=RetentionTime, Intensity)
    readr::write_tsv(reports_distinct, path="data/exp1depletion/uplc_raw.tsv", na="")
  }
  
  #
  # Copy related files to a separate directory
  #
  cmd = list()
  if(!file.exists(depletion_output_dir.arw)) cmd[["create_arw"]] = sprintf("mkdir %s", depletion_output_dir.arw)
  cmd[["copy_reports"]] = sprintf("cp %s/%s %s/%s", depletion_input_dir.arw, depletion_input_dir.report_files.expected, depletion_output_dir.arw, gsub("^d", "report_bioaccumulation", depletion_export_samples$report))
  # cmd[["copy_eic"]] = sprintf("data/exp2metabolomics/data.eic.csv %s/raw.zip", depletion_output_dir)
  writeLines(unlist(cmd), con=depletion_output.collect_script)
  
  #
  # Prepare metabolights samples
  #
  
  # TODO: seems ready
  # Samples (depletion)
  metabolights_export_samples.depletion = depletion_export_samples %>%
    dplyr::distinct(AssayName, Plate, Extraction, GroupName, Status, WellReplicate, .keep_all=T) %>%
    dplyr::mutate(`Sample Name`=SampleName.long) %>%
    dplyr::mutate(`Source Name`="Cultivated bacteria")  %>%
    dplyr::mutate(`Protocol REF[Sample collection]`="Sample collection") %>%
    # dplyr::mutate(`Protocol REF[Extraction]`="Supernatant") %>%
    dplyr::mutate(
      `Characteristics[Organism]`=dplyr::case_when(
        species.long == "Standard" | Status %in% c("ctrl") ~ "experimental blank",
        T ~ species.long),
      `Term Source REF[Organism]`=dplyr::case_when(
        species.long == "Standard" | Status %in% c("ctrl") ~ "MTBLS",
        T ~ "NCBITaxon"),
      `Term Accession Number[Organism]`=dplyr::case_when(
        species.long == "Standard" | Status %in% c("ctrl") ~ "http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218",
        T ~ paste0("http://purl.obolibrary.org/obo/NCBITaxon_", species.taxid))) %>%
    dplyr::mutate(
      `Characteristics[Organism part]`=dplyr::case_when(
        species.long == "Standard" | Status %in% c("ctrl") ~ "experimental control",
        T ~ "whole organism"),
      `Term Source REF[Organism part]`=dplyr::case_when(
        species.long == "Standard" | Status %in% c("ctrl") ~ "MTBLS",
        T ~ "CARO"),
      `Term Accession Number[Organism part]`=dplyr::case_when(
        species.long == "Standard" | Status %in% c("ctrl") ~ "http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218",
        T ~ "http://purl.obolibrary.org/obo/CARO_0000064")) %>%
    dplyr::mutate(
      `Characteristics[Variant]`="normal variability",
      `Term Source REF[Variant]`="PATO",
      `Term Accession Number[Variant]`="http://purl.obolibrary.org/obo/PATO_0045074") %>%
    dplyr::mutate(
      `Characteristics[Sample type]`=dplyr::case_when(
        species.long == "Standard" ~ "standard sample",
        Status=="ctrl" ~ "control sample (no bacteria)",
        Status=="GMM" ~ "GMM media reference (no drug)",
        Status=="sample" ~ "sample"),
      `Term Source REF[Sample type]`="AFO",
      `Term Accession Number[Sample type]`=dplyr::case_when(
        species.long == "Standard" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000255",
        Status=="ctrl" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000253",
        Status=="GMM" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000258",
        Status=="sample" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000171")) %>%
    dplyr::mutate(
      `Factor[Extraction]`=dplyr::case_when(
        species.long == "Standard" ~ "standard sample",
        T ~ "lysate"),
      `Term Source REF[Extraction factor]`=dplyr::case_when(
        species.long == "Standard" ~ "AFO",
        T ~ "OBI"),
      `Term Accession Number[Extraction factor]`=dplyr::case_when(
        species.long == "Standard" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000255",
        T ~ "http://purl.obolibrary.org/obo/OBI_1000036")
    ) %>%
    dplyr::mutate(
      `Factor[Sample type]`=dplyr::case_when(
        species.long == "Standard" ~ "standard sample",
        Status=="ctrl" ~ "control sample (no bacteria)",
        Status=="GMM" ~ "GMM media reference (no drug)",
        Status=="sample" ~ "sample"),
      `Term Source REF[Sample type factor]`="AFO",
      `Term Accession Number[Sample type factor]`=dplyr::case_when(
        species.long == "Standard" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000255",
        Status=="ctrl" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000253",
        Status=="GMM" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000258",
        Status=="sample" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000171")
    ) %>%
    dplyr::mutate(
      `Factor[Drug]`=dplyr::case_when(
        species.long == "Standard" ~ "no drug",
        T ~ GroupName.long),
      `Term Source REF[Drug factor]`=dplyr::case_when(
        species.long == "Standard" ~ "MTBLS",
        T ~ "CHEBI"),
      `Term Accession Number[Drug factor]`=dplyr::case_when(
        species.long == "Standard" ~ "http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218", 
        T ~ paste0("http://purl.obolibrary.org/obo/CHEBI_", GroupName.chebi))) %>%
    dplyr::mutate(
      `Factor[Organism]`=Plate.long,
      `Term Source REF[Organism factor]`=dplyr::case_when(
        species.long == "Standard"~"AFO",
        T ~ "NCBITaxon"),
      `Term Accession Number[Organism factor]`=dplyr::case_when(
        species.long == "Standard" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000255",
        T ~ paste0("http://purl.obolibrary.org/obo/NCBITaxon_", species.taxid))) %>%
    dplyr::mutate(
      `Factor[Screen]`="Depletion assay",
      `Term Source REF[Screen factor]`="NCIT",
      `Term Accession Number[Screen factor]`="http://purl.obolibrary.org/obo/NCIT_C48261") %>%
    dplyr::select(dplyr::matches("Sample Name|Source Name|Protocol REF|Source Name|Characteristics|Factor|Term Source REF|Term Accession Number")) %>%
    data.frame(check.names=F)

  
  # Samples (bioaccumulation)
  metabolights_export_samples.bioaccumulation = bioaccumulation_export_samples %>%
    dplyr::distinct(AssayName, Plate, Extraction, GroupName, Status, WellReplicate, .keep_all=T) %>%
    dplyr::mutate(`Sample Name`=SampleName.long) %>%
    dplyr::mutate(`Source Name`="Cultivated bacteria")  %>%
    # dplyr::mutate(`Protocol REF[Extraction]`=ifelse(Extraction == "Supernatant", "Supernatant", "Lysed cells with media")) %>%
    dplyr::mutate(`Protocol REF[Sample collection]`="Sample collection") %>%
    dplyr::mutate(
      `Characteristics[Organism]`=dplyr::case_when(
        Status %in% c("ctrl") ~ "experimental blank",
        T ~ species.long),
      `Term Source REF[Organism]`=dplyr::case_when(
        Status %in% c("ctrl") ~ "MTBLS",
        T ~ "NCBITaxon"),
      `Term Accession Number[Organism]`=dplyr::case_when(
        Status %in% c("ctrl") ~ "http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218",
        T ~ paste0("http://purl.obolibrary.org/obo/NCBITaxon_", species.taxid))) %>%
    dplyr::mutate(
      `Characteristics[Organism part]`=dplyr::case_when(
        Status %in% c("ctrl") ~ "experimental control",
        T ~ "whole organism"),
      `Term Source REF[Organism part]`=dplyr::case_when(
        Status %in% c("ctrl") ~ "MTBLS",
        T ~ "CARO"),
      `Term Accession Number[Organism part]`=dplyr::case_when(
        Status %in% c("ctrl") ~ "http://www.ebi.ac.uk/metabolights/ontology/MTBLS_000218",
        T ~ "http://purl.obolibrary.org/obo/CARO_0000064")) %>%
    dplyr::mutate(
      `Characteristics[Variant]`="normal variability",
      `Term Source REF[Variant]`="PATO",
      `Term Accession Number[Variant]`="http://purl.obolibrary.org/obo/PATO_0045074") %>%
    dplyr::mutate(
      `Characteristics[Sample type]`=dplyr::case_when(
        Status=="ctrl" ~ "control sample (no bacteria)",
        Status=="GMM" ~ "GMM media reference (no drug)",
        Status=="sample" ~ "sample"),
      `Term Source REF[Sample type]`="AFO",
      `Term Accession Number[Sample type]`=dplyr::case_when(
        Status=="ctrl" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000253",
        Status=="GMM" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000258",
        Status=="sample" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000171")
    ) %>%
    dplyr::mutate(
      `Factor[Extraction]`=dplyr::case_when(
        Extraction == "Supernatant" ~ "supernatant",
        T ~ "lysate"),
      `Term Source REF[Extraction factor]`="OBI",
      `Term Accession Number[Extraction factor]`=dplyr::case_when(
        Extraction == "Supernatant" ~ "http://purl.obolibrary.org/obo/OBI_0000034",
        T ~ "http://purl.obolibrary.org/obo/OBI_1000036")
    ) %>%
    dplyr::mutate(
      `Factor[Sample type]`=dplyr::case_when(
        Status=="ctrl" ~ "control sample (no bacteria)",
        Status=="GMM" ~ "GMM media reference (no drug)",
        Status=="sample" ~ "sample"),
      `Term Source REF[Sample type factor]`="AFO",
      `Term Accession Number[Sample type factor]`=dplyr::case_when(
        Status=="ctrl" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000253",
        Status=="GMM" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000258",
        Status=="sample" ~ "http://purl.allotrope.org/ontologies/role#AFRL_0000171")
    ) %>%
    dplyr::mutate(
      `Factor[Drug]`= GroupName.long,
      `Term Source REF[Drug factor]` = "CHEBI",
      `Term Accession Number[Drug factor]`=paste0("http://purl.obolibrary.org/obo/CHEBI_", GroupName.chebi)) %>%
    dplyr::mutate(
      `Factor[Organism]`=Plate.long,
      `Term Source REF[Organism factor]`="NCBITaxon",
      `Term Accession Number[Organism factor]`=paste0("http://purl.obolibrary.org/obo/NCBITaxon_", species.taxid)) %>%
    dplyr::mutate(
      `Factor[Screen]`="Bioactivity assay",
      `Term Source REF[Screen factor]`="NCIT",
      `Term Accession Number[Screen factor]`="http://purl.obolibrary.org/obo/NCIT_C48261") %>%
    dplyr::select(dplyr::matches("Sample Name|Source Name|Protocol REF|Source Name|Characteristics|Factor|Term Source REF|Term Accession Number")) %>%
    data.frame(check.names=F)
  
  metabolights_export_samples.all = rbind(metabolights_export_samples.depletion, metabolights_export_samples.bioaccumulation)
  metabolights_export_samples.final = metabolights_export_samples.all %>% setNames(gsub("(Term Source REF|Term Accession Number|Protocol REF)\\[.*", "\\1", names(.)))
  readr::write_tsv(metabolights_export_samples.final, path="C:/Users/sandrejev/Desktop/UPLC2/december3/s_MTBLS1264.txt") # path="data/exp2metabolomics/metabolights/s_MTBLS1264.txt",  na="")
  
  
  #
  # Prepare Metabolights Assays
  #
  metabolights_export_assays.depletion = depletion_export_samples %>%
    # dplyr::distinct(AssayName, Plate, Extraction, GroupName, Status, WellReplicate, .keep_all=T) %>%
    dplyr::mutate(`Sample Name`=SampleName.long)  %>%
    dplyr::mutate(`Protocol REF[Extraction]`="Supernatant") %>%
    dplyr::mutate(`Parameter Value[Post Extraction]`="") %>%
    dplyr::mutate(`Parameter Value[Derivatization]`="") %>%
    dplyr::mutate(`Extract Name`="Supernatant") %>%
    dplyr::mutate(`Protocol REF[Chromatography]`="UPLC") %>%
    dplyr::mutate(`Parameter Value[Chromatography Instrument]`="Vanquish UHPLC system") %>%
    dplyr::mutate(`Term Source REF[Chromatography Instrument]`="OBI") %>%
    dplyr::mutate(`Term Accession Number[Chromatography Instrument]`="http://purl.obolibrary.org/obo/OBI_0000485") %>%
    dplyr::mutate(`Parameter Value[Autosampler model]`="") %>%
    dplyr::mutate(`Parameter Value[Column model]`="ACQUITY UPLC CSH C18 column (Waters, Part number 186005297; 130Å, 1.7 µm, 2.1 mm X 100 mm)") %>%
    dplyr::mutate(`Parameter Value[Column type]`="reverse phase") %>%
    dplyr::mutate(`Parameter Value[Guard column]`="") %>%
    dplyr::mutate(`Parameter Value[Detector]`="ACQUITY-UPLC-PDA-Detektor", `Term Source REF[Detector]`="MS", `Term Accession Number[Detector]`="http://purl.obolibrary.org/obo/MS_1000818") %>%
    dplyr::mutate(`Parameter Value[Signal range]`="190 to 500 nm") %>%
    dplyr::mutate(`Parameter Value[Resolution]`="1.2 nm") %>%
    dplyr::mutate(`Labeled Extract Name`="", `Label`="", `Term Source REF[Label]`="", `Term Accession Number[Label]`="") %>%
    dplyr::mutate(`Raw Spectral Data File`=paste0("Depletion/", SampleName.long,"-channel", Channel.Id, ".arw")) %>% 
    dplyr::mutate(`Protocol REF[Normalization]`="") %>%
    dplyr::mutate(`Normalization Name`="") %>%
    dplyr::mutate(`Derived Spectral Data File`="") %>%
    dplyr::mutate(`Protocol REF[Data Transformation]`="Data Transformation") %>%
    dplyr::mutate(`MS Assay Name`=`Sample Name`) %>%
    dplyr::mutate(`Data Transformation Name`="Metabolite identification") %>%
    dplyr::mutate(`Metabolite Assignment File`="m_MTBLS1264_LC-DAD___metabolite_profiling_v2_maf.tsv") %>%
    dplyr::select(`Sample Name`:`Metabolite Assignment File`) %>%
    data.frame(check.names=F)
  
  metabolights_export_assays.bioaccumulation = bioaccumulation_export_samples %>%
    # dplyr::distinct(AssayName, Plate, Extraction, GroupName, Status, WellReplicate, .keep_all=T) %>%
    dplyr::mutate(`Sample Name`=SampleName.long)  %>%
    dplyr::mutate(`Protocol REF[Extraction]`=ifelse(Extraction == "Supernatant", "Supernatant", "Lysed cells with media")) %>%
    dplyr::mutate(`Parameter Value[Post Extraction]`="") %>%
    dplyr::mutate(`Parameter Value[Derivatization]`="") %>%
    dplyr::mutate(`Extract Name`=ifelse(Extraction == "Supernatant", "Supernatant", "Lysed cells with media")) %>%
    dplyr::mutate(`Protocol REF[Chromatography]`="UPLC") %>%
    dplyr::mutate(`Parameter Value[Chromatography Instrument]`="Vanquish UHPLC system", `Term Source REF[Chromatography Instrument]`="OBI", `Term Accession Number[Chromatography Instrument]`="http://purl.obolibrary.org/obo/OBI_0000485") %>%
    dplyr::mutate(`Parameter Value[Autosampler model]`="") %>%
    dplyr::mutate(`Parameter Value[Column model]`="ACQUITY UPLC CSH C18 column (Waters, Part number 186005297; 130?, 1.7 µm, 2.1 mm X 100 mm)") %>%
    dplyr::mutate(`Parameter Value[Column type]`="reverse phase") %>%
    dplyr::mutate(`Parameter Value[Guard column]`="") %>%
    dplyr::mutate(`Parameter Value[Detector]`="ACQUITY-UPLC-PDA-Detektor", `Term Source REF[Detector]`="MS", `Term Accession Number[Detector]`="http://purl.obolibrary.org/obo/MS_1000818") %>%
    dplyr::mutate(`Parameter Value[Signal range]`="190 to 500 nm") %>%
    dplyr::mutate(`Parameter Value[Resolution]`="1.2 mm") %>%
    dplyr::mutate(`Labeled Extract Name`="", `Label`="", `Term Source REF[Label]`="", `Term Accession Number[Label]`="") %>%
    dplyr::mutate(`Raw Spectral Data File`=paste0("Bioactivity/", SampleName.long,"-channel", Channel.Id, ".arw")) %>% 
    dplyr::mutate(`Protocol REF[Normalization]`="") %>%
    dplyr::mutate(`Normalization Name`="") %>%
    dplyr::mutate(`Derived Spectral Data File`="") %>%
    dplyr::mutate(`Protocol REF[Data Transformation]`="Data Transformation") %>%
    dplyr::mutate(`Data Transformation Name`="Metabolite identification") %>%
    dplyr::mutate(`Metabolite Assignment File`="m_MTBLS1264_LC-DAD___metabolite_profiling_v2_maf.tsv") %>%
    dplyr::select(`Sample Name`:`Metabolite Assignment File`) %>%
    data.frame(check.names=F)
  
  
  
  # table(basename(metabolights_export_assays.depletion$`Raw Spectral Data File`) %in% list.files("data/exp1depletion/metabolights/arw_filtered"))
  # table(basename(metabolights_export_assays.bioaccumulation$`Raw Spectral Data File`) %in% list.files("data/exp2metabolomics/metabolights/arw_filtered"))
  # 
  
  # readr::write_tsv(metabolights_export_assays.depletion %>% setNames(gsub("(Term Source REF|Term Accession Number|Protocol REF)\\[.*", "\\1", names(.))), path="data/exp1depletion/metabolights/a_MTBLS1264_LC-DAD___metabolite_profiling.txt",  na="") # 
  # readr::write_tsv(metabolights_export_assays.bioaccumulation %>% setNames(gsub("(Term Source REF|Term Accession Number|Protocol REF)\\[.*", "\\1", names(.))), path="data/exp2metabolomics/metabolights/a_MTBLS1264_LC-DAD___metabolite_profiling.txt",  na="") # path="data/exp2metabolomics/metabolights/a_MTBLS1264_LC-DAD___metabolite_profiling.txt",  na="\"\"")
  readr::write_tsv(metabolights_export_assays.depletion %>% setNames(gsub("(Term Source REF|Term Accession Number|Protocol REF)\\[.*", "\\1", names(.))), path="C:/Users/sandrejev/Desktop/UPLC2/december3/a_MTBLS1264_LC-DAD___metabolite_profiling.txt",  na="\"\"") # 
  readr::write_tsv(metabolights_export_assays.bioaccumulation %>% setNames(gsub("(Term Source REF|Term Accession Number|Protocol REF)\\[.*", "\\1", names(.))), path="C:/Users/sandrejev/Desktop/UPLC2/december3/a_MTBLS1264_LC-DAD___metabolite_profiling-3.txt",  na="\"\"") # path="data/exp2metabolomics/metabolights/a_MTBLS1264_LC-DAD___metabolite_profiling.txt",  na="\"\"")
  
  
  #
  # MAF (Metabolites)
  #
  
  drugs_annotations = depletion_export_samples %>% 
    dplyr::filter(!is.na(drug.long)) %>%
    dplyr::distinct(drug.long, drug.chebi, drug.formula, drug.smiles, drug.inchi, drug.uplc_retention, GroupName) 
  drugs_annotations = rbind(drugs_annotations, data.frame(drug.long="Cafffeine", drug.chebi="27732", drug.formula="C8H10N4O2", drug.smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C", drug.inchi="InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3", drug.uplc_uv=230, GroupName="caffeine"))
  maf_drugs_annotations = drugs_annotations %>%
    dplyr::mutate(
      M1_database_identifier=drug.chebi, 
      M1_chemical_formula=drug.formula, 
      M1_smiles=drug.smiles, 
      M1_inchi=drug.inchi, 
      M1_metabolite_identification=drug.long,
      M1_mass_to_charge="",
      M1_fragmentation="",
      M1_modifications="",
      M1_charge="",
      M1_retention_time=drug.uplc_retention,
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
    dplyr::select(dplyr::matches("M1_"), GroupName) %>%
    setNames(gsub("^M1_", "", colnames(.)))
  
  
  #
  # MAF (depletion)
  #
  depletion_metabolites = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp1depletion/peaks_exported.tsv")    
  maf_depletion = depletion_metabolites %>%
    dplyr::mutate(SampleName.long=paste0("Depletion-", Plate, "-", Extraction, "-", GroupName, "-", Status, "-", WellReplicate)) %>%
    reshape2::melt(measure.vars=c("Area", "Area_IS"), value.name="Area", variable.name="metabolite_short") %>%
    dplyr::mutate(GroupName=ifelse(metabolite_short=="Area", GroupName, "caffeine")) %>%
    dplyr::group_by(SampleName.long, GroupName) %>%
    dplyr::summarise(Area=mean(Area, na.rm=T))
  
  maf_depletion_values.wide = maf_depletion %>%
    reshape2::dcast(GroupName~SampleName.long, value.var="Area")
  maf_bioaccumulation_values.wide = bioaccumulation_export_samples %>%
    reshape2::melt(measure.vars=c("Area", "Area.caffee"), value.name="Area") %>%
    dplyr::mutate(GroupName=ifelse(variable=="Area.caffee", "caffeine", GroupName)) %>%
    reshape2::dcast(GroupName ~ SampleName.long, value.var="Area")
  
  maf_final = maf_drugs_annotations %>%
    dplyr::left_join(maf_depletion_values.wide, by=c("GroupName")) %>%
    dplyr::left_join(maf_bioaccumulation_values.wide, by=c("GroupName")) %>%
    dplyr::select(-GroupName)
  
  readr::write_tsv(maf_final, "C:/Users/sandrejev/Desktop/UPLC2/december3/m_MTBLS1264_LC-DAD___metabolite_profiling_v2_maf.tsv", na = "")

  
}