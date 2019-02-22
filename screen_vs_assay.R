library(ggplot2)
library(xlsx)
library(dplyr)
source("functions.R")

bug_map.new = readr::read_delim("data/bug_map.tsv", "\t") %>% dplyr::filter(!is.na(species.NT_ID)) %>% dplyr::select(species.short, species.long, species.code, species.color)
drug_map.new = readr::read_delim("data/drug_map.tsv", "\t") %>% dplyr::select(drug.code, drug.short, drug.long, drug.uplc_excluded, drug.color3)

load("data/exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData")
data.assaylong = read.csv("data/exp2metabolomics/data.depletionmodeassay_long.csv") # ~/Workspace/170610_save/computational stuff/UPLC analysis/data.depletionmodeassay_long.csv
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
  dplyr::mutate(Drug=drug.long, MedianDiff.assay=MedianDiff.assay*100, MedianDiff.screen=MedianDiff.screen*100, SEMDiff.assay=SEMDiff.assay*100, SEMDiff.screen=SEMDiff.screen*100) %>%
  dplyr::group_by(species.short, species.code, Drug, Status) %>%
  dplyr::summarise(MedianDiff.screen=unique(MedianDiff.screen), SEMDiff.screen=unique(SEMDiff.screen), MedianDiff.assay=unique(MedianDiff.assay), SEMDiff.assay=unique(SEMDiff.assay))
  

#
# Assay vs screen (all together)
#
corscreenassay.test = cor.test(data.merge$MedianDiff.screen, data.merge$MedianDiff.assay, use="complete.obs", method="pearson")
corscreenassay.estimate = round(corscreenassay.test$estimate, digits=2)
corscreenassay = paste0("R=", corscreenassay.estimate, " (pval < ", gsub("([0-9]e-0?)", "e-", sprintf("%.0e", corscreenassay.test$p.value)), ")")

pdf("reports/screen_vs_assay.pdf", height=12, width=9)
myColors = drug_map.new$drug.color3; names(myColors) = drug_map.new$drug.long
ggplot(data.merge, aes(x=MedianDiff.screen,y=MedianDiff.assay, color=Drug)) +
  geom_point(aes(fill=paste0(stringr::str_pad(paste0(species.code, ":"), width=7, side="right"), "\t", species.short)), size=0.6) +
  geom_errorbarh(aes(xmax = MedianDiff.screen+SEMDiff.screen, xmin = MedianDiff.screen-SEMDiff.screen), height=2, size=0.4) +
  geom_errorbar(aes(ymax = MedianDiff.assay+SEMDiff.assay, ymin = MedianDiff.assay-SEMDiff.assay), width=2, size=0.4) +
  geom_text(aes(label=species.code), size=2, hjust=-0.2, vjust=-0.5) +
  geom_label(x=60, y=90, label=corscreenassay, size=6, color="black") +
  geom_abline(slope=corscreenassay.estimate , intercept=0, linetype="dotted", color="black") +
  coord_cartesian(ylim=c(-105,105), xlim = c(-105,105)) +
  scale_fill_discrete('Species', guide=guide_legend(override.aes=list(alpha=0, size=1))) +
  scale_color_manual(values=myColors) +
  guides(fill=guide_legend(ncol=3)) +
  labs(x="Median difference to control (%)(screen)", y="Median difference to  (%) (validation)") +
  myTheme +
  theme(axis.text=element_text(size=18), legend.position="bottom", legend.box = "vertical")
dev.off()

readr::write_tsv(data.merge, "reports/screen_vs_assay.tsv", col_names=T, na="")  