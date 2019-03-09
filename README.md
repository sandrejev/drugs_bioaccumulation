Analyzis from Klünemann et al. "Bioaccumulation of therapeutic drugs by human gut bacteria" publication.

## Docker image
Dockerfile provides a recepie to reproduce our setup. A functional docker image with newest code is also avalable from Docker cloud:

~~~~
docker pull sandrejev/drugs_bioaccumulation
docker run -it sandrejev/drugs_bioaccumulation bash # to inspect contents of the container
docker run -it sandrejev/drugs_bioaccumulation # To repeat the analysis
~~~~

## Data files

* **drug_map.csv** General information about all the drugs used in this study
* **bug_map.csv** General information about all the bacteria used in this study
* **go_terms.csv** List of all GO terms for enrichment testing
* **screenG_tax_info_specI_clusters.tab** Core gut microbiome [1] composition average abundance and prevalence
* **screenG_tax_input_specI_clusters.tab** Core gut microbiome [1] cluster-to-strain (taxid) mapping
* **kegg.db** SQLite file containing information about all human gut bacteria and enzymes present in them taken from KEGG

#### Growth assay
* **exp0growth/2016-11-28_curves_annotation.tab** Processed growth curves data (maxOD, lag, mu) collected for all bacteria/drug combinations
* **exp0growth/curves.rel_annotation_2016-11-28.tab** Relative (to no-drug control) growth curves data (maxOD, lag, mu) collected for all bacteria/drug combinations

#### Bacteria-drug interaction screen. 
* **exp1depletion/170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData** Data containing depletion of drugs from the media when a species is grown in presence (sample) and without (control) drug

#### Screen validation and bioaccumulation detection
* **exp2metabolomics/data.depletionmodeassay_long.csv** Drug UPLC measurements in absence (control) and in presence (sample) of species. Data collected from supernatant and total inoculum
* **exp2metabolomics/hits.csv** Relative drug UPLC measurements. Data collected from supernatant and total inoculum

#### NMR-based validation of duloxetine bioaccumulation and mass-spectrometry based secreted metabolite analysis
* **exp3metabolomics/extraxset.fillden.noInj1.RData** Untargeted metabolomics of extracellular composition
* **exp3metabolomics/lysxset.fillden.noInj1.RData** Untargeted metabolomics of lysed cells
* **exp3metabolomics/clickxset.fillden.noInj1.RData** Click-chemistry based pull down of duloxetine binding proteins
* **exp3metabolomics/KEGG_metabolites_bugs.csv**
* **exp3metabolomics/KEGG_pathways.csv**

#### Mass spectrometry-based protein identification
* **exp3proteomics/Data_all.txt** Quantitative abundance of etected proteins
* **exp3proteomics/uniprot_C_saccharolyticum_enzymes.tab** Uniprot annotation of *C. saccharolyticum* proteins

#### Untargeted metabolomics analysis of conditioned media assay
* **exp4spentmedia/Salivarius_duloxetine.csv** Metabolomics data from untargeted metabolomics analysis of *E. rectale* growth on *S. salivarius* spent media

#### Sensitivity of the five species used in the community assembly experiment
* **exp5concentration/161118_dilution24h_rep[1-3].xlsx** Growth data with different concentration of duloxetine
* **exp5concentration/well2species.tsv** Sample to condition mapping

#### Conditioned Media Assay an Community assembly assay
* **exp6transfers/160916_160907_160823_transferassaydulox.xlsx** Metabolomics data collected from Conditioned Media Assay
* **exp6transfers/otu_table_wo_mono_normalized_by_gene_copy_mod.csv** Metagenomic data for Community assembly assay
* **exp6transfers/Sample_map_mod.csv** Sample to condition mapping for Community assembly assay

#### References
[1] Mende, D. R., Sunagawa, S., Zeller, G. & Bork, P. Accurate and universal delineation of prokaryotic species. Nat. Methods 10, 881–884 (2013).


## Packages used in the analysis
* RSQLite:2.1.1
* ggrepel:0.8.0.9000
* scales:1.0.0
* stringi:1.2.4
* imputeLCMD:2.0
* pcaMethods:1.70.0
* norm:1.0-9.5
* tmvtnorm:1.4-10
* gmm:1.6-2
* sandwich:2.5-0
* Matrix:1.2-12
* impute:1.52.0
* gplots:3.0.1
* Mfuzz:2.38.0
* DynDoc:1.56.0
* widgetTools:1.56.0
* e1071:1.7-0
* rcdk:3.4.7.1
* rcdklibs:2.0
* rJava:0.9-10
* xcms:3.0.2
* MSnbase:2.4.2
* ProtGenerics:1.10.0
* mzR:2.12.0
* Rcpp:1.0.0
* qvalue:2.10.1
* mutoss:0.1-12
* mvtnorm:1.0-8
* purrr:0.2.5
* magrittr:1.5
* reshape2:1.4.3
* corpcor:1.6.9
* limma:3.34.9
* sva:3.26.0
* BiocParallel:1.12.0
* genefilter:1.60.0
* stringr:1.3.1
* locfit:1.5-9.1
* vsn:3.46.0
* DESeq2:1.18.1
* SummarizedExperiment:1.8.1
* DelayedArray:0.4.1
* matrixStats:0.54.0
* GenomicRanges:1.30.3
* GenomeInfoDb:1.14.0
* IRanges:2.12.0
* S4Vectors:0.16.0
* ade4:1.7-11
* boot:1.3-20
* LSD:4.0-0
* plyr:1.8.4
* multtest:2.34.0
* Biobase:2.38.0
* BiocGenerics:0.24.0
* readr:1.1.1
* networkD3:0.4
* rCharts:0.4.2
* data.table:1.11.4
* vegan:2.5-2
* lattice:0.20-35
* permute:0.9-4
* labdsv:1.8-0
* cluster:2.0.6
* MASS:7.3-50
* mgcv:1.8-23
* nlme:3.1-131
* widyr:0.1.1
* RColorBrewer:1.1-2
* tidyr:0.8.1
* bindrcpp:0.2.2
* dplyr:0.7.6
* xlsx:0.6.1
* ggplot2:3.0.0.9000
