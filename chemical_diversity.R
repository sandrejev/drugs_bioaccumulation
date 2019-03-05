dir.create("reports", showWarnings=F)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(widyr)
library(labdsv)
library(vegan)
library(data.table)

drug_classes = readr::read_delim("data/chemical_diversity/2064_Drug_class_combined.dms", "\t", col_names=c("drug.id", "drug.class"))
drug_ids = readr::read_delim("data/chemical_diversity/drugs.tsv", "\t")
drug_map = readr::read_delim("data/drug_map.tsv", "\t") %>% dplyr::mutate(in_study=T) %>% dplyr::select(drug.long, drug.uplc_excluded, in_study) %>% unique()

#
# Distance matrix (from long)
#
AT = read.table("data/chemical_diversity/2064_MCS.dms", sep = "\t", quote = "", header = T)
vals = sort(unique(c(as.character(AT$Drug1), as.character(AT$Drug2))))
AT.mat = matrix(NA, nrow=length(vals), ncol=length(vals), dimnames=list(vals, vals))
diag(AT.mat) = 1    
AT.mat[as.matrix(AT[, 1:2])] <- AT[,3]
AT.mat[as.matrix(AT[, 2:1])] <- AT[,3]
AT.mat = 1-AT.mat

#
# Prepare PCO data for visualization
#
AT.pco = labdsv::pco(AT.mat, k=2)
AT.eigenvalues = vegan::eigenvals(AT.pco) 
AT.variance = AT.eigenvalues / sum(AT.eigenvalues)
AT.pco_df = as.data.frame(AT.pco$points) %>%
  dplyr::transmute(drug.id=rownames(.), PC1=V1, PC2=V2) %>%
  dplyr::inner_join(drug_classes, by="drug.id") %>%
  dplyr::inner_join(drug_ids, by="drug.id") %>%
  dplyr::left_join(drug_map, by="drug.long") %>%
  dplyr::mutate(
    drug.study_name=ifelse(!is.na(in_study), drug.long, ""),
    drug.study_name=ifelse(!is.na(drug.uplc_excluded), paste(drug.study_name, "(excluded)"), drug.study_name), 
    in_study=ifelse(!is.na(in_study), "Selected", "Not selected")
  ) 

#
# Plot
#
pdf("reports/chemical_diversity.pdf", width=10, height=9)
ggplot(aes(PC1, PC2, shape=in_study, size=in_study, alpha=in_study, fill=drug.class, color=in_study), data=AT.pco_df) +
  geom_point() +
  geom_text(aes(label=drug.study_name), hjust=-0.12, vjust=0.1, col='black', size=4, show.legend=F) +
  theme_bw(base_size=18) +
  scale_size_manual(values=c("Selected"=5, "Not selected"=2.3)) +
  scale_shape_manual(values=c("Selected"=21, "Not selected"=21)) +
  scale_alpha_manual(values=c("Selected"=1, "Not selected"=0.25)) +  
  scale_colour_manual(values=c("Selected"="#000000", "Not selected"="#FFFFFF")) +
  labs(alpha="Part of the study", fill="ATC", xlab=paste0("PCO1 (",round(round(100*AT.variance[1])),"% variance explained)"), ylab=paste0("PCO2 (",round(round(100*AT.variance[2])),"% variance explained)")) +
  guides(shape=F, color=F, size=F, alpha=guide_legend(override.aes=list(shape=16, size=5, color="#FF0000")), fill=guide_legend(override.aes=list(shape=21, size=5, color="#FFFFFF00")))
dev.off()
