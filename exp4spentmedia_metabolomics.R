dir.create("reports", showWarnings=F)
library(reshape2)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
source("functions.R")

exp4spentmedia.metabolomics = function()
{
  #################### Select Hits for Drug ########################
  drug = read.csv(file="data/exp4spentmedia/Salivarius_duloxetine.csv", header = T, sep =",", quote="", na.strings = "")
  
  ##### Data preprocessing #####
  
  ### Put the raw data in tidy format
  drug_comb = drug %>%
    dplyr::rename(X13_5="X13.5", X13_5.1="X13.5.1", X13_5.2="X13.5.2") %>%
    reshape2::melt(id.vars="mz", value.name="intensity") %>%
    dplyr::mutate(
      sample=gsub("X([^.]*).*", "\\1", variable),
      species.short=ifelse(sample %in% c("0", "10"), "S. salivarius", "E. rectale"),
      replicate=gsub("X([^.]*)\\.?(.*)", "\\2", variable),
      replicate=ifelse(replicate=="", "0", replicate))
  
  ### Calculate time point means across all replicates
  drug_comb.avg = drug_comb %>%
    dplyr::group_by(mz, sample, species.short) %>%
    dplyr::summarise(intensity=ifelse(any(intensity==0) & !all(intensity==0), mean(intensity[which(intensity>0)]), mean(intensity)))
  
  
  ##### Initial filtering #####
  
  ### Select only metabolites that get newly produced by the first species:
  # Select only ions that equal zero at time point 0 and not zero at time point 13_5
  drug_comb.filtZ = drug_comb.avg %>% dplyr::group_by(mz) %>% dplyr::filter(any(sample == "0" & intensity < 0.5))
  drug_comb.filtZ2 = drug_comb.filtZ %>% dplyr::group_by(mz) %>% dplyr::filter(any(sample == "13_5" & intensity > 0.5))
  
  ### Filter noisy data; metabolites with potential drop outs in the time series:
  # Filter ions which are zero at time point 30 and have the same (+/-10 percent) intensity in 
  # time points 13_5 and 36
  drug_comb.filtTP = drug_comb.filtZ2 %>%
    dplyr::group_by(mz) %>%
    dplyr::filter(!any((intensity[sample == "13_5"] + 0.1*intensity[sample == "13_5"]) > intensity[sample == "36"] &
                         (intensity[sample == "13_5"] - 0.1*intensity[sample == "13_5"]) < intensity[sample == "36"] &
                         (sample == "30" & intensity < 0.5) 
    ))
  
  # Filter ions which are zero at time point 36 and have the same (+/-10 percent) intensity in 
  # time points 30 and 47
  drug_comb.filtTP2 = drug_comb.filtTP %>%
    dplyr::group_by(mz) %>%
    dplyr::filter(!any((intensity[sample == "30"] + 0.1*intensity[sample == "30"]) > intensity[sample == "47"] &
                         (intensity[sample == "30"] - 0.1*intensity[sample == "30"]) < intensity[sample == "47"] &
                         (sample == "36" & intensity < 0.5)
    ))
  
  ##### Subsequent filtering to obtain the biologically relevant hits. ##### 
  ##### Hits are metabolites that get newly produced by the first and consumed by the second species. ##### 
  
  ### Filter metabolites that do not get consumed by the second species:
  # Filter ions that do not change over time
  drug_steady.ions = drug_comb.filtTP2 %>%
    dplyr::group_by(mz) %>%
    dplyr::filter(!any((abs(intensity[sample == "30"] - intensity[sample == "13_5"]) < 1) &
                         (abs(intensity[sample == "36"] - intensity[sample == "30"]) < 1) &
                         (abs(intensity[sample == "47"] - intensity[sample == "36"]) < 1) &
                         (intensity[sample == "36"] - intensity[sample == "13_5"] > -1) &
                         (intensity[sample == "47"] - intensity[sample == "30"] > -1) &
                         (intensity[sample == "47"] - intensity[sample == "13_5"] > -1 )))
  
  ### Filter metabolites that get only produced but not consumed by the second species:
  # Filter ions that increase after time point 13_5
  drug_increasing.ions = drug_steady.ions %>%
    dplyr::group_by(mz) %>%
    dplyr::mutate(is_hit=!any((intensity[sample == "30"] - intensity[sample == "13_5"]) > 1))
  
  # Filter ions that increase after time point 30 with no change at time point 13.5 
  drug_increasing.ions2 = drug_increasing.ions %>%
    dplyr::group_by(mz) %>%
    dplyr::mutate(is_hit=is_hit & !any(abs(intensity[sample == "30"] - intensity[sample == "13_5"]) < 1 &
                         (intensity[sample == "36"] - intensity[sample == "30"]) > 1))
  
  # Filter ions that increase after time point 36 with no change at time point 13.5 and 30
  drug_increasing.ions3 = drug_increasing.ions2 %>%
    dplyr::group_by(mz) %>%
    dplyr::mutate(is_hit=is_hit & !any(abs(intensity[sample == "30"] - intensity[sample == "13_5"]) < 1 &
                         abs(intensity[sample == "36"] - intensity[sample == "30"]) < 1 &
                         (intensity[sample == "36"] - intensity[sample == "13_5"]) > -1 &
                         (intensity[sample == "47"] - intensity[sample == "36"]) > 1))
  
  ### Filter any metabolite that gets produced by the second species even if it was consumed before:
  # Filter all ions that reach zero at time point 30 or time point 36 and increase afterwards again
  # to obtain the final hits
  drug_finalHits = drug_increasing.ions3 %>%
    dplyr::group_by(mz) %>%
    dplyr::mutate(is_hit=is_hit & !any((intensity[sample == "30"] == 0 &
                          (intensity[sample == "36"] > 0 | intensity[sample == "47"] > 0)) |
                         (intensity[sample == "36"] == 0 & intensity[sample == "47"] > 0) ))
  
  
  pdf("reports/exp3spentmedia_metabolomics.pdf")
  ggplot(drug_finalHits %>% dplyr::filter(is_hit), aes(x=sample, y=intensity, group=mz, color=species.short)) + 
    geom_line(alpha=0.4, size=0.5) +
    scale_color_manual(values=c("S. salivarius"="#4B7FAE", "E. rectale"="#E24F52")) +
    facet_wrap(~is_hit) +
    myTheme
  dev.off()
  
  file_ToSave = drug_finalHits %>% 
    dplyr::select(mz, is_hit, sample, intensity) %>%
    tidyr::spread(key=sample, value=intensity) %>% 
    dplyr::mutate(is_hit=ifelse(is_hit, "Yes", "No"))
  
  write.table(file_ToSave, file = "reports/exp4spentmedia_metabolomics.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
}