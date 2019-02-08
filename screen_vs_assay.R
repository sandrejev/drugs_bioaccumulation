library(ggplot2)
library(xlsx)
library(dplyr)
source("functions.R")

bug_map.new = readr::read_delim("data/bug_map.tsv", "\t") %>% dplyr::filter(!is.na(species.NT_ID)) %>% dplyr::select(species.short, species.long, species.code, species.color)
drug_map.new = readr::read_delim("data/drug_map.tsv", "\t") %>% dplyr::select(drug.code, drug.short, drug.long, drug.uplc_excluded, drug.color3)

load("data/exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
data.assaylong <- read.csv("~/Workspace/170610_save/computational stuff/UPLC analysis/data.depletionmodeassay_long.csv")
data.assaylong = data.assaylong %>% 
  dplyr::mutate(Ctrl=gsub("smpl","sample",Ctrl)) %>%
  dplyr::left_join(drug_map.new, by=c("Method"="drug.code")) %>% 
  dplyr::left_join(bug_map.new, by=c("Bugs"="species.code")) %>% 
  dplyr::filter(is.na(drug.uplc_excluded))

data.screen = data.clean %>% 
  dplyr::filter(Status!="GMM") %>%
  dplyr::left_join(drug_map.new, by=c("GroupName"="drug.short")) %>% 
  dplyr::left_join(bug_map.new, by=c("Species.x"="species.long")) %>% 
  dplyr::filter(is.na(drug.uplc_excluded)) %>%
  dplyr::group_by(drug.long, species.short, Status) %>%
  dplyr::mutate(MedianDiff=median(DiffCtrlSample,na.rm=T), SEMDiff=sd(DiffCtrlSample,na.rm=T)/sqrt(sum(!is.na(DiffCtrlSample)))) %>%
  data.frame()

data.SN = data.assaylong %>%
  dplyr::filter(Extraction=="Supernatant" & Ctrl!="zero" & Condition!="med") %>%
  dplyr::group_by(drug.long, species.short, Ctrl) %>%
  dplyr::mutate(MedianDiff=median(DiffCtrlSample,na.rm=T), SEMDiff=sd(DiffCtrlSample,na.rm=T)/sqrt(sum(!is.na(DiffCtrlSample)))) %>%
  data.frame()

# 'data.merge should have all possible sample replicate combinations, not just means for biological replicates'
data.merge = data.screen %>%
  dplyr::inner_join(data.SN, by=c("species.short", "drug.long", "Status"="Ctrl")) %>% 
  dplyr::filter(Status=="sample") %>%
  setNames(gsub("(\\.x)$", ".assay", names(.), perl=T)) %>%
  setNames(gsub("(\\.y)$", ".screen", names(.), perl=T)) %>%
  dplyr::mutate(MedianDiff.assay=MedianDiff.assay*100, MedianDiff.screen=MedianDiff.screen*100, SEMDiff.assay=SEMDiff.assay*100, SEMDiff.screen=SEMDiff.screen*100)

#
# Assay vs screen (all together)
#
corscreenassay = with(data.merge, round(cor(MedianDiff.screen, MedianDiff.assay, use="complete.obs", method="pearson"), digits=2))

pdf("reports/screen_vs_assay.pdf", height=6, width=8)
myColors = drug_map.new$drug.color3; names(myColors) = drug_map.new$drug.long
ggplot(data.merge, aes(x=MedianDiff.screen,y=MedianDiff.assay, color=drug.long)) +
  geom_point(size=0.5) +
  geom_errorbarh(aes(xmax = MedianDiff.screen+SEMDiff.screen, xmin = MedianDiff.screen-SEMDiff.screen), height=2.5, size=0.1) +
  geom_errorbar(aes(ymax = MedianDiff.assay+SEMDiff.assay, ymin = MedianDiff.assay-SEMDiff.assay), width=2.5, size=0.1) +
  coord_cartesian(ylim=c(-105,105), xlim = c(-105,105)) +
  scale_color_manual(values=myColors) +
  annotate("label", x=60, y=90, label=paste0("Pearson Correlation: ",corscreenassay)) +
  labs(x="Depletion screen", y="Metabolomics assay")
  geom_abline(slope=corscreenassay,intercept=0, linetype="dotted", color="black") +
  myTheme
dev.off()

readr::write_tsv(data.merge, "reports/screen_vs_assay.tsv", col_names=T, na="")  