library(dplyr)
library(readr)
library(CHNOSZ)
library(Rdisop)


exp2metabolomics_extend_supplementary = function() {
  setwd("F:/home/andrejev/Workspace/drugs_bioaccumulation")
  
  kegg2compounds = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/db/KEGG/kegg_compounds.tsv", col_names=T) %>%
    tidyr::separate_rows(KEGG_METABOLITE_NAME, sep=";") %>%
    dplyr::distinct(KEGG_CPD, KEGG_FORMULA) %>%
    dplyr::mutate(KEGG_CPD=paste0("cpd:", KEGG_CPD), KEGG_FORMULA=as.character(KEGG_FORMULA))
  
  adducts = readr::read_delim("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/db/adducts.tsv", "\t", comment="#") %>%
    tidyr::separate(adduct_formula, into=c("adduct_formula_add", "adduct_formula_sub"), sep="-") %>%
    dplyr::filter(adduct_name %in% c("[M-H]-", "[M+ACN+H]+", "[M+H]+"))
  
  # save(r1, file="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r1.rda")
  # save(r2, file="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r2.rda")
  # save(r3, file="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r3.rda")
  # save(r4, file="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r4.rda")
  # save(r5, file="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r5.rda")
  # for(i in 1:5) {
  #   load(stringr::str_glue("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r{i}.rda", i=i))
  # }
  # r = rbind(r1,r2,r3,r4,r5) %>% dplyr::ungroup()
  # save(r, file="C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r_x.rda")
  # load("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r.rda")
  
  # 44826
  
  electron_mass = 0.00054858026 
  r1 = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp13lcms/Supplementary_Tables_14_4Sergej.txt") %>%
    dplyr::mutate(KEGG_1_ID=`KEGG ID`) %>%
    tidyr::separate_rows(KEGG_1_ID, sep=",") %>%
    dplyr::left_join(kegg2compounds %>% dplyr::mutate(KEGG_IN_DB=T) %>% dplyr::select(KEGG_CPD, KEGG_IN_DB), by=c("KEGG_1_ID"="KEGG_CPD")) %>%
    dplyr::arrange(dplyr::desc(KEGG_IN_DB)) %>%
    dplyr::group_by(Peak, Preparation, Bacteria) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(kegg2compounds, by=c("KEGG_1_ID"="KEGG_CPD")) %>%
    tidyr::crossing(adducts) %>%
    # dplyr::slice(00001:10000) %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z){
      zz<<-z
      z = data.frame(z, stringsAsFactors=F, check.names=F)
      
      # if(z$Preparation=="Extracellular" & z$Bacteria=="Bacteroides uniformis" & z$Peak=="111.1/15") {
      #   asdasdsadad
      # }
      
      # z = data.frame(
      #   KEGG_FORMULA="C3H6O3S", adduct_formula_sub="H", check.names = F, stringsAsFactors = F
      # )
      # 
      if(!is.na(z$KEGG_FORMULA)) {
        print(paste(z$KEGG_FORMULA, z$adduct_formula_sub, z$adduct_formula_ad))
        element_matrix = dplyr::bind_rows(
          as.data.frame(t(CHNOSZ::makeup(z$KEGG_FORMULA, count.zero=T))),
          as.data.frame(t(CHNOSZ::makeup(z$adduct_formula_sub, count.zero=T)))
        ) %>% replace(is.na(.), 0)
        element_matrix.diff = min(element_matrix[1,] - element_matrix[2,])
        element_matrix.sum = sum(element_matrix[1,] - element_matrix[2,])
        if(element_matrix.diff >= 0 & element_matrix.sum > 0) {
          z.molecule = Rdisop::getMolecule(z$KEGG_FORMULA)
          z.exact_mass = z.molecule$isotopes[[1]][1,1]
          
          if(z$adduct_nmol > 1) {
            for(i in 2:z$adduct_nmol) { z.molecule = Rdisop::addMolecules(z.molecule$formula, z$KEGG_FORMULA) }
          }
          
          if(!is.na(z$adduct_formula_add) && z$adduct_formula_add!="") {
            z.molecule = Rdisop::addMolecules(z.molecule$formula, z$adduct_formula_add)
          }
          if(!is.na(z$adduct_formula_sub) && z$adduct_formula_sub!="") {
            z.molecule = Rdisop::subMolecules(z.molecule$formula, z$adduct_formula_sub)
          }
          ion_mz = (z.molecule$isotopes[[1]][1,] - z$adduct_charge*electron_mass)/abs(z$adduct_charge)
          z$exact_mass = z.molecule$exactmass
        } else {
          z$exact_mass = NA_real_
        }
      } else {
        z$exact_mass = NA_real_
      }
      z
    })(.)) 
  
  load("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/r_x.rda")
  rr = r %>%
    dplyr::ungroup() %>%
    dplyr::mutate(delta_ppm=(`m/z (median)`-electron_mass-exact_mass)/exact_mass*1e6) %>%
    dplyr::arrange(abs(delta_ppm)) %>%
    dplyr::group_by(Preparation, Bacteria, Peak) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  rr.out = rr %>% 
    dplyr::select(-dplyr::matches("^adduct_"), adduct_name, -KEGG_IN_DB, KEGG_1_ID)
  
  readr::write_tsv(rr.out, "C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp3metabolomics_supplementary_table.tsv")
}


exp2metabolomics_preprocess = function() {
  bioaccumulation_processed_map = readr::read_delim("data/exp2metabolomics/data.depletionmodeassay_long.csv", ",") %>%
    dplyr::mutate(order=1:n())
  
  bioaccumulation_raw = data.frame()
  for(f in list.files("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp2metabolomics/raw", full.names=T, pattern="*.csv")) {
    d = readr::read_csv(f) %>%
      dplyr::select(Name, RT, Area, Height, SampleName, Condition, Extraction, Method, Replicate, Vial, Sample.Set.Name=`Sample Set Name`, Channel.Description=`Channel Description`)
    bioaccumulation_raw = rbind(bioaccumulation_raw, d)
  }
  
  bioaccumulation_raw.caffee = bioaccumulation_raw %>%
    dplyr::filter(grepl("_caffe", Name)) %>%
    dplyr::group_by(Sample.Set.Name, SampleName) %>%
    dplyr::summarise(Area.caffee=mean(Area, na.rm=T))
  
  
  bioaccumulation_raw.processed = bioaccumulation_raw %>%
    dplyr::filter(grepl("_peak", Name)) %>%
    dplyr::select(Sample.Set.Name, SampleName, RT, Area.peak=Area) %>%
    dplyr::inner_join(bioaccumulation_raw.caffee, by=c("Sample.Set.Name", "SampleName"))
  
  bioaccumulation_processed_map.new = bioaccumulation_processed_map %>%
    dplyr::inner_join(bioaccumulation_raw.processed, by=c("Sample.Set.Name", "SampleName", "RT")) %>%
    dplyr::group_by(Sample.Set.Name, SampleName, RT, order) %>%
    dplyr::arrange(abs(Area-Area.peak)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Area.caffee=tidyr::replace_na(Area.caffee, 0)) %>%
    dplyr::arrange(order) %>%
    dplyr::select(-Area.peak, -order)
  
  write.table(bioaccumulation_processed_map.new, "C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/exp2metabolomics/data.depletionmodeassay_long_new.csv", sep=",", quote=T, qmethod="double", row.names=F)
  
  for(s in colnames(bioaccumulation_processed_map)) {
    writeLines(paste0(s, ": ", sum(is.na(bioaccumulation_processed_map[[s]]) & is.na(bioaccumulation_processed_map.new[[s]]) | bioaccumulation_processed_map[[s]]==bioaccumulation_processed_map.new[[s]], na.rm=T), "/", nrow(bioaccumulation_processed_map)))
  }
}