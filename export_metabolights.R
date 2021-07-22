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

exp1depletion.export_metabolights = function()
{
  #
  # Load Species and drug data
  #
  drug_map = readr::read_delim("data/drug_map.tsv", "\t")
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")
  
  load("data/exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
  # plates_with_growth = unique(data.clean %>% dplyr::filter(Status=="GMM" & growth=="growth") %>% .$Plate)
  # data.degrad = data.clean %>% dplyr::filter(Status=="sample" & growth=="growth" & dummy==1 & Plate %in% plates_with_growth)
  uplc_interaction_map = readr::read_tsv("data/exp2metabolomics/raw_annotations/uplc_interactions_mapping.tsv")
  plate_mapping = readr::read_tsv("data/exp2metabolomics/raw_annotations/uplc_plata_mapping.tsv")
  
  #
  # Load mapping between raw files and final data
  #
  growth_file2plate = readr::read_csv("data/exp0growth/160106_FileToPlate.csv")
  growth_curves_annotations = readr::read_tsv("data/exp0growth/raw_annotations/curves_annotation_2016-11-28.tab") %>%
    dplyr::inner_join(growth_file2plate, by="File")
  
  #
  # Load one mapping between raw files and final data that was missing in the main mapping (Plate_Martina: Buniformis66A2_Buniformis65A2)
  #
  uplc_interaction_map_cols = paste0("^(", paste(c("AssayName", "SampleName", "Sample.Set.Name", "Vial", "Plate_Martina", "Acq.Method.Set", "Channel.Description", "Channel.Id"), collapse="|"), ")$")
  uplc_interaction_map_Buniformis66A2_Buniformis65A2 = readr::read_tsv("data/exp2metabolomics/raw_annotations/report157342.txt") %>%
    dplyr::mutate(AssayName="Screen", Plate_Martina="Buniformis66A2_Buniformis65A2") %>%
    setNames(gsub(" ", ".", colnames(.)))
  uplc_interaction_map_extended = dplyr::bind_rows(
    uplc_interaction_map %>% dplyr::filter(Plate_Martina!="Buniformis66A2_Buniformis65A2"),
    uplc_interaction_map_Buniformis66A2_Buniformis65A2) %>% 
    dplyr::select(dplyr::matches(uplc_interaction_map_cols))
  
  #
  # Prepare annotations 
  #
  annotations = data.clean %>%
    dplyr::distinct(GroupName, Plate, Batch, drug, Status, Well, Replicate, Species.x) %>%
    dplyr::group_by(GroupName, Plate, drug, Status, Species.x) %>% 
    dplyr::mutate(Replicate=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::rename(species.long="Species.x") 
  
  #
  # Prepare mappings
  #
  channel_mappings = uplc_interaction_map_extended %>% 
    dplyr::filter(AssayName=="Screen") %>% 
    dplyr::filter(!grepl("Spect", Channel.Description)) %>%
    # dplyr::filter(!grepl("down|purging", Acq.Method.Set)) %>%
    dplyr::group_by(AssayName, SampleName, Sample.Set.Name, Vial) %>%
    dplyr::mutate(ChannelNo=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Plate_Martina_Number=as.numeric(gsub(":.*", "", Vial)), Well=gsub(".*:", "", Vial)) %>% 
    dplyr::inner_join(plate_mapping, by=c("Plate_Martina", "Plate_Martina_Number")) %>%
    dplyr::mutate(drug.something=gsub("^[^_]+_", "", Acq.Method.Set)) %>%
    dplyr::mutate(drug.something=dplyr::case_when(drug.something=="eyet"~"ezet", drug.something=="rosuva"~"rosu", drug.something=="donez"~"done", T~drug.something)) %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short, drug.long1=drug.long), by=c("drug.something"="drug.short")) %>%
    dplyr::left_join(drug_map %>% dplyr::select(drug.short2, drug.long2=drug.long), by=c("drug.something"="drug.short2")) %>%
    dplyr::mutate(drug.long=dplyr::case_when(!is.na(drug.long2)~drug.long2, !is.na(drug.long1)~drug.long1, T~NA_character_)) %>%
    dplyr::select(-drug.long1, -drug.long2) %>%
    dplyr::mutate(Status=dplyr::case_when(
      grepl("ctrl", SampleName)~"ctrl",
      grepl("GMM", SampleName)~"GMM",
      grepl("down|purging", Acq.Method.Set)~"technical",
      T~"sample")) %>%
    dplyr::mutate(GroupName=ifelse(Status=="technical", NA_character_, substr(Acq.Method.Set, 1, 4))) %>%
    dplyr::mutate(drug.long=ifelse(Status=="GMM", NA_character_, drug.long)) %>%
    dplyr::group_by(Plate, GroupName, Status) %>%
    dplyr::mutate(Replicate=match(paste(SampleName, Vial), unique(paste(SampleName, Vial)))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Extraction="Supernatant") %>%
    dplyr::select(AssayName, Sample.Set.Name, Extraction, Plate_Martina, SampleName, GroupName, Plate, Vial, Replicate, Status, drug.long, Channel.Description, Channel.Id, ChannelNo)  # Replicate, mz, RetentionTime, species.short
  
  # Pyri is excluded from drugs
  output = annotations %>%
    dplyr::left_join(channel_mappings %>% dplyr::mutate(has_data=T), by=c("Plate", "GroupName", "Status", "Replicate")) %>%
    # dplyr::mutate(species.long=dplyr::case_when(Status=="ctrl"~NA_character_, species=="Standard"~NA_character_, T~species.long)) %>%
    dplyr::select(AssayName, Extraction, Sample.Set.Name, Plate, SampleName, GroupName, Status, Replicate, Vial, drug.long, species.long, Channel.Description, Channel.Id, has_data) 
  stopifnot(any(!is.na(output$has_data)))
  readr::write_tsv(output %>% dplyr::select(-has_data), path="data/exp1depletion/uplc_channels_annotation.tsv")
  
  if(F) {
    input_dir.report = "/g/patil/Share/EmpowerProjects/DepletionModeAssay_Rawraw/Raw_readable"
    input_dir.raw = "/g/patil/Share/EmpowerProjects/Martina" 
    # input_dir.report = "/g/scb2/patil/andrejev/UPLC_depletion/reports"
    # input_dir.raw = "/g/scb2/patil/andrejev/UPLC_depletion/raw" 
    output_dir = "/g/scb2/patil/andrejev/UPLC_depletion"
    
    input_dir.report = "data/exp1depletion/arw"
    input_dir.raw = "data/exp1depletion/raw" 
    output_dir = "data/exp1depletion/raw_filtered"
    
    # Check missing report files
    input_dir.report_files = list.files(input_dir.report)
    input_dir.report_files.expected = sprintf("d%s.arw", xy$Channel.Id)
    writeLines(sprintf("Missing report files: %i", length(setdiff(input_dir.report_files.expected, input_dir.report_files))))
    
    input_dir.raw_files = list.files(input_dir.raw)
    input_dir.raw_files.expected = sprintf("d%s.dat", xy$Channel.Id)
    writeLines(sprintf("Missing RAW files: %i", length(setdiff(input_dir.raw_files.expected, input_dir.raw_files))))
    
    
    #
    # Collect all reports into single file
    #
    i.run = 0
    reports = data.frame()
    reports.100 = data.frame()
    for(i in output$Channel.Id) {
      if(i %in% unique(reports$Channel.Id)) next
      i.file = sprintf(sprintf("%s/d%i.arw", input_dir.report, i))
      if(file.exists(i.file)) {
        i.run = i.run+1
        print(sprintf("%i/%i (added: %i)", which(i==unique(output$Channel.Id)), length(unique(output$Channel.Id)), i.run))
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
    
    
    
    #
    # Copy related files to a separate directory
    #
    cmd = list()
    if(!file.exists(output_dir)) cmd[["create_dir"]] = sprintf("mkdir %s", output_dir)
    if(!file.exists(sprintf("%s/raw", output_dir))) cmd[["create_dir/raw"]] = sprintf("mkdir %s/raw", output_dir)
    # if(!file.exists(sprintf("%s/reports", output_dir))) cmd[["create_dir/raw"]] = sprintf("mkdir %s/reports", output_dir)
    # cmd[["copy_reports"]] = sprintf("cp %s/report_raw%i.arw %s/reports", input_dir.report, xy$Channel.Id, output_dir)
    # cmd[["zip_reports"]] = sprintf("zip -r %s/reports.zip %s/reports", output_dir, output_dir)
    # cmd[["copy_eic"]] = sprintf("data/exp2metabolomics/data.eic.csv %s/raw.zip", output_dir)
    cmd[["copy_raw"]] = sprintf("cp %s/d%i.dat %s/raw/", input_dir.raw, xy$Channel.Id, output_dir)
    cmd[["zip_raw"]] = sprintf("zip -r %s/raw.zip %s/raw", output_dir, output_dir)
    writeLines(unlist(cmd), con=sprintf("%s/collect_uplc_depletion.sh", output_dir))
  }
}

exp2metabolomics.export_metabolights = function()
{
  drug_map = readr::read_delim("data/drug_map.tsv", "\t")
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")
  
  uplc_interaction = readr::read_delim("data/exp2metabolomics/data.depletionmodeassay_long.csv", ",") %>%
    # dplyr::filter(!is.na(Drugs)) %>%
    dplyr::mutate(species.code2=substr(SampleName, 1, 2), species.code3=substr(SampleName, 1, 3)) %>%
    dplyr::distinct(Extraction, Sample.Set.Name, SampleName, Channel.Description, Vial, mz=X., RetentionTime=RT, drug.short2=Drugs, Replicate, species.code2, species.code3, Ctrl) %>% 
    dplyr::left_join(drug_map %>% dplyr::select(drug.short2, drug.long), by="drug.short2") %>% 
    dplyr::left_join(bug_map %>% dplyr::select(species.code, species.short2=species.short), by=c("species.code2"="species.code")) %>% 
    dplyr::left_join(bug_map %>% dplyr::select(species.code, species.short3=species.short), by=c("species.code3"="species.code")) %>%
    dplyr::mutate(species.short=dplyr::case_when(!is.na(species.short2)~species.short2, !is.na(species.short3)~species.short3, T~NA_character_)) %>%
    dplyr::mutate(Status=dplyr::case_when(Ctrl=="ctrl"~"ctrl", Ctrl=="smpl"~"sample", Ctrl=="zero"~"GMM")) %>%
    dplyr::distinct(Extraction, Sample.Set.Name, SampleName, Channel.Description, Status, Replicate, Vial, mz, RetentionTime, drug.long, species.short)
  
  uplc_interaction_map = readr::read_tsv("data/exp2metabolomics/raw_annotations/uplc_interactions_mapping.tsv")
  channel_mappings = uplc_interaction_map %>% 
    dplyr::mutate(GroupName=gsub("^[0-9]+_|_martina", "", Acq.Method.Set)) %>%
    group_by(AssayName, SampleName, Sample.Set.Name, Channel.Description, Vial) %>%
    dplyr::mutate(ChannelNo=1:n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(AssayName, SampleName, Sample.Set.Name, Plate=Plate_Martina, GroupName, Channel.Description, Vial, ChannelNo, Channel.Id)
  
  annotations = uplc_interaction %>% 
    group_by(SampleName, Sample.Set.Name, Channel.Description, Vial, Replicate) %>%
    dplyr::mutate(ChannelNo=1:n()) %>%
    dplyr::ungroup()
  
  output = channel_mappings %>% 
    dplyr::inner_join(annotations, by=c("SampleName", "Sample.Set.Name", "Channel.Description", "Vial", "ChannelNo")) %>%
    dplyr::select(-ChannelNo) %>%
    dplyr::select(AssayName, Extraction, Sample.Set.Name, Plate, SampleName, GroupName, Status, Replicate, Vial, drug.long, species.short, Channel.Description, Channel.Id) 
  readr::write_tsv(output, path="data/exp2metabolomics/uplc_channels_annotation.tsv")
  
  
  #
  # Collect all the raw files in one folder
  #
  if(F) {
    input_dir.report = "/g/patil/Share/EmpowerProjects/DepletionModeAssay_Rawraw/Raw_readable"
    input_dir.raw = "/g/patil/Share/EmpowerProjects/Martina"
    output_dir = "/g/scb2/patil/andrejev/UPLC"
    input_dir.report = "data/exp2metabolomics/raw_reports/"
    # input_dir.raw = "/g/patil/Share/EmpowerProjects/Martina" 
    # output_dir = "/g/scb2/patil/andrejev/UPLC"
    
    
    #
    # Collect all reports into single file
    #
    i.run = 0
    reports = data.frame()
    reports.100 = data.frame()
    for(i in unique(output$Channel.Id)) {
      if(i %in% reports$Channel.Id) next
      i.file = sprintf(sprintf("%s/report_raw%i.arw", input_dir.report, i))
      if(file.exists(i.file)) {
        i.run = i.run+1
        report.i = readr::read_tsv(i.file, skip=2, col_names=c("RetentionTime", "Intensity"), col_types=readr::cols()) %>% dplyr::mutate(Channel.Id=i)
        reports.100 = rbind(reports.100, report.i)
        if(i.run %% 100==0) {
          print(sprintf("%i/%i (added: %i)", which(i==unique(output$Channel.Id)), length(unique(output$Channel.Id)), i.run))
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
    readr::write_tsv(reports_distinct, path="data/exp2metabolomics//uplc_raw.tsv", na="")
    
    
    cmd = list()
    if(!file.exists(output_dir)) cmd[["create_dir"]] = sprintf("mkdir %s", output_dir)
    if(!file.exists(sprintf("%s/raw", output_dir))) cmd[["create_dir/raw"]] = sprintf("mkdir %s/raw", output_dir)
    if(!file.exists(sprintf("%s/reports", output_dir))) cmd[["create_dir/raw"]] = sprintf("mkdir %s/reports", output_dir)
    cmd[["copy_reports"]] = sprintf("cp %s/report_raw%i.arw %s/reports", input_dir.report, xy$Channel.Id, output_dir)
    cmd[["copy_raw"]] = sprintf("cp %s/d%i.dat %s/raw/", input_dir.raw, xy$Channel.Id, output_dir)
    cmd[["zip_reports"]] = sprintf("zip -r %s/reports.zip %s/reports", output_dir, output_dir)
    cmd[["zip_raw"]] = sprintf("zip -r %s/raw.zip %s/raw", output_dir, output_dir)
    cmd[["copy_eic"]] = sprintf("data/exp2metabolomics/data.eic.csv %s/raw.zip", output_dir)
    writeLines(unlist(cmd), con="collect_uplc.sh")
  }
}
