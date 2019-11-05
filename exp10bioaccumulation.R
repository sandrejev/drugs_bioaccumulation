library(dplyr)
library(readr)
library(ggplot2)
library(readxl)
library(reshape2)
library(tidyverse)
library(stringr)
library(httr)
library(jsonlite)
library(xml2)
library(dendextend)
library(Rtsne)

exp10bioaccumulation.collect_data = function() {
  data = readr::read_tsv("data/exp10bioaccumulation/data.tsv", na="NF")
  data_long = data %>%
    #dplyr::filter(!(Batch=="Batch 3" & grepl("IAI1", Species_name))) %>%
    dplyr::mutate(Row=gsub("[0-9]", "", Well), Collumn=as.numeric(gsub("[A-H]", "", Well)), Collumn=factor(Collumn, sort(unique(Collumn)))) %>%
    dplyr::group_by(Batch, Plate) %>%
    data.frame() %>% 
    reshape2::melt(measure.vars=c("Duloxetine", "Montelukast", "Roflumilast", "Rosiglitazone"), variable.name="Drug_measured", value.name="Intensity") %>%
    dplyr::mutate(Intensity=dplyr::case_when(Concentration==0 & is.na(Intensity)~0, is.na(Intensity)~NA_real_, T~Intensity)) %>%
    dplyr::mutate(None_IS=1e6) 
  data_long = data_long %>%
    reshape2::melt(measure.vars=colnames(data_long)[grepl("_IS", colnames(data_long))], variable.name="Standard", value.name="StandardIntensity") %>%
    dplyr::mutate(IntensityNorm=Intensity/StandardIntensity) %>%
    dplyr::select(Batch, Plate, Well, No, Extraction, Species_name, Incubated, Drug, Concentration, Replicate, Drug_measured, Intensity, IntensityNorm, Standard, StandardIntensity)
  readr::write_delim(data_long, "data/exp10bioaccumulation/data_long.tsv", delim="\t", na="")
}


exp10bioaccumulation.analyze = function() {
  # Are there "DMSO only" samples?
  data_long = readr::read_tsv("data/exp10bioaccumulation/data_long.tsv") %>%
    dplyr::filter(!is.na(IntensityNorm) & Incubated=="Yes") %>%
    dplyr::mutate(BatchPlate=paste0("Plate: ", gsub("Batch ", "", Batch), "/", gsub("Plate ", "", Plate))) 
                  
  #
  # Select only drug concentration correction samples
  #
  data_long.0 = data_long %>% 
    dplyr::filter(Drug==Drug_measured & Species_name=="no bug" & Concentration > 0) %>% 
    data.frame()
  
  #
  # Plot performance of different internal standards
  #
  data_long.0.cor = data_long.0 %>%
    dplyr::group_by(Drug, Standard) %>%
    dplyr::mutate(y=0.01+0.2*max(IntensityNorm)**as.numeric(gsub("Plate ", "", Plate)), x=0) %>%
    dplyr::group_by(Batch, Plate, BatchPlate, Extraction, Drug, Species_name, Standard, x, y) %>% 
    dplyr::summarise(n=length(Plate), R=cor(IntensityNorm, Concentration, use="pairwise.complete.obs"), Slope=lm(Concentration~IntensityNorm+0)$coefficients[1])
  
  pdf("reports/exp10bioaccumulation_nobug_calibration.pdf", width=10, height=14)
  ggplot(data_long.0, aes(x=IntensityNorm, y=Concentration)) +
    geom_smooth(aes(color=BatchPlate), method="lm", formula=y~x, size=1, alpha=.2) + 
    geom_point(aes(group=Replicate, color=BatchPlate)) + 
    geom_line(aes(group=Replicate, color=BatchPlate)) + 
    #geom_text(aes(x=x, y=y, label=paste0(gsub("Plate ", "P", Plate), " R= ", round(R, 2), " (Slope:", round(Slope), ")")), data=data_long.0.cor, hjust=0) +
    facet_grid(Standard~Species_name+Extraction, scales="free_y")
  dev.off()
  
  
  #
  # Use drug concentration correction from plate #1 to correct plate #2 because later plate didn't have it's own concentration correction data
  #
  data_long.0x = data_long.0 %>% 
    dplyr::filter(Batch=="Batch 3" & Plate=="Plate 1") %>% 
    dplyr::mutate(Plate="Plate 2", BatchPlate="Plate: 1/2") %>% data.frame()
  data_long.0x = rbind(data_long.0, data_long.0x) %>%
    dplyr::rename(IntensityNorm.0="IntensityNorm")
  data_long.0x_mean = data_long.0x %>%
    dplyr::select(-Species_name) %>%
    dplyr::group_by(BatchPlate, Drug, Standard, Extraction, Concentration) %>%
    dplyr::summarise(IntensityNorm.0=median(IntensityNorm.0, na.rm=T))
  
  #
  # Density plot of internal standards
  #
  ggplot(data_long %>% dplyr::filter(!grepl("None", Standard))) +
    geom_density(aes(x=StandardIntensity, fill=factor(Standard)), alpha=0.3) +
    geom_vline(aes(xintercept=Intensity, color=factor(Concentration)), data=data_long.0 %>% dplyr::group_by(Concentration) %>% dplyr::summarise(Intensity=mean(Intensity))) + 
    facet_wrap(~Batch + BatchPlate, scales="free_y", ncol=1)
  
  #
  # Caclulate sum of concetration change between different extraction methods
  #
  data_long.1 = data_long %>% 
    dplyr::filter(Species_name!="no bug" & Drug==Drug_measured) %>%
    #dplyr::inner_join(data_long.0x, by=c("Batch", "Plate", "Drug_measured", "Incubated", "BatchPlate", "Drug", "Standard", "Extraction", "Concentration", "Species_name")) %>%
    dplyr::group_by(Batch, Plate, BatchPlate, Drug, Drug_measured, Incubated, Extraction, Species_name, Standard) %>%
    dplyr::do((function(z){
      if(nrow(z) < 3) return(data.frame(DetectedConcentration=NA))
      zz<<-z
      z.0 = data_long.0x %>% dplyr::filter(Batch==z$Batch[1] & Plate==z$Plate[1] & Drug==z$Drug[1] & Standard==z$Standard[1] & Extraction==z$Extraction[1])
      z.lm = lm(Concentration~IntensityNorm.0+0, data=z.0)
      z$DetectedConcentration = predict(z.lm, z %>% dplyr::mutate(IntensityNorm.0=IntensityNorm))
      z$InitialConcentration = z$Concentration
      z
    })(.)) %>%
    dplyr::ungroup()
  
  
  pdf("reports/exp10bioaccumulation_effect.pdf", width=10, height=8)
  data_long.1.ggplot = data_long.1 %>%
    dplyr::filter(dplyr::between(InitialConcentration, 0, 70) & Standard=="Fluoxetine_IS") %>%
    dplyr::mutate(InitialConcentration=factor(InitialConcentration)) %>%
    dplyr::select(Batch, Plate, Species_name, Drug, Extraction, Standard, Well, InitialConcentration, DetectedConcentration, IntensityNorm)
  data_long.1.ggplot_test = data_long.1.ggplot %>% 
    dplyr::group_by(Batch, Species_name, Drug, Standard) %>%
    dplyr::mutate(IntensityNormMax=max(DetectedConcentration), IntensityNormMax50=max(DetectedConcentration[InitialConcentration=="50"])) %>%
    dplyr::group_by(Batch, Species_name, Drug, InitialConcentration, Standard, IntensityNormMax, IntensityNormMax50) %>%
    dplyr::do((function(z) {
      zz<<-z
      pval = tryCatch({ 
        t.test(z$DetectedConcentration[z$Extraction=="total"], z$DetectedConcentration[z$Extraction=="supernatant"])$p.value 
      }, finally=function(t) NA_real_)
      data.frame(pval=pval)
    })(.)) %>% 
    dplyr::filter(pval<0.05) %>%
    dplyr::ungroup()
  
  f.concentrations = c("30", "40", "50", "70")
  data_long.1.ggplot.f = data_long.1.ggplot %>% 
    dplyr::filter(((Batch=="Batch 3" & grepl("sacch", Species_name)) | (Batch=="Batch 4" & grepl("IAI1", Species_name))) & InitialConcentration %in% f.concentrations & Drug=="Duloxetine") %>%
    dplyr::mutate(InitialConcentration=factor(InitialConcentration, f.concentrations))
  data_long.1.ggplot_test.f = data_long.1.ggplot_test %>% 
    dplyr::inner_join(data_long.1.ggplot.f %>% dplyr::distinct(Batch, Species_name, Drug, Extraction, Standard), by=c("Batch", "Species_name", "Drug", "Standard")) %>%
    dplyr::mutate(InitialConcentration=factor(InitialConcentration, f.concentrations))
  ggplot(data_long.1.ggplot.f) +
    geom_boxplot(aes(y=DetectedConcentration, x=InitialConcentration, fill=Extraction), width=.5) +
    #geom_errorbarh(aes(y=IntensityNormMax50*1.05, xmin=as.numeric(InitialConcentration)-0.1, xmax=as.numeric(InitialConcentration)+0.1), data=data_long.1.ggplot_test.f, height=0) + 
    labs(x="Initial duloxetine concentration", y="Normalized duloxetine peak area") +
    facet_wrap(Batch~Species_name, scales="free") +
    scale_fill_manual(values=c(supernatant="#d3d2d2", total="#656263")) +
    theme_classic(base_size=16) + 
    theme(strip.background=element_blank(), aspect.ratio=1)
  for(d in unique(data_long.1.ggplot$Drug)) {
    print(ggplot(data_long.1.ggplot %>% dplyr::filter(Drug==d)) +
      geom_boxplot(aes(y=DetectedConcentration, x=InitialConcentration, fill=Extraction), width=.5) +
      #geom_errorbarh(aes(y=IntensityNormMax*1.05, xmin=as.numeric(InitialConcentration)-0.1, xmax=as.numeric(InitialConcentration)+0.1), data=data_long.1.ggplot_test %>% dplyr::filter(Standard==Standard), height=0) + 
      labs(x=paste0("Initial ", d, " concentration"), y=paste0("Normalized ", d, " peak area")) +
      facet_grid(Standard~Batch+Species_name, scales="free") +
      scale_fill_manual(values=c(supernatant="#d3d2d2", total="#656263")) +
      theme_classic(base_size=16) + 
      theme(strip.background=element_blank(), aspect.ratio=1))
  }
  dev.off()
  
  
  
  pdf("reports/exp10bioaccumulation_standards_bias.pdf", width=10, height=10)
  ggplot(data_long.1, aes(x=InitialConcentration, y=DetectedConcentration)) + 
    # geom_boxplot(aes(x=factor(InitialConcentration), color=Extraction)) + 
    # geom_hline(yintercept=1, color="#666666", alpha=0.5, size=2) +
    geom_line(color="#666666", alpha=0.5, size=2, data=data.frame(InitialConcentration=sort(unique(data_long$Concentration)), DetectedConcentration=sort(unique(data_long$Concentration)))) +
    geom_line(aes(color=Extraction, group=paste(Extraction, Replicate, Standard)), alpha=0.3) +
    #geom_line(aes(color=Extraction, group=paste(Extraction, Standard)), size=2, data=data_long.1_mean) +
    geom_point(aes(fill=Plate), color="#FFFFFF00", shape=21, alpha=0.3) +
    facet_wrap(Species_name~Batch+Standard, scales="free_y")
  dev.off()
}
