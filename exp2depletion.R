library(ggplot2)
library(xlsx)
library(dplyr)
source("functions.R")

scatter.total2supernatant = function() 
{
  bug_map = readr::read_delim("data/bug_map.tsv", "\t") %>% dplyr::filter(!is.na(species.NT_ID))
  drug_map = readr::read_delim("data/drug_map.tsv", "\t")
  
  data.assaylong = readr::read_delim("data/exp2metabolomics/data.depletionmodeassay_long.csv", ",") %>%
    dplyr::mutate(Ctrl=gsub("smpl","sample", Ctrl))  %>%
    dplyr::filter(!(Method=="di" & Bugs=="Fn")) %>% 
    dplyr::left_join(bug_map, by=c("Bugs"="species.code")) %>% 
    dplyr::rename(species.code="Bugs") %>%
    dplyr::inner_join(drug_map, by=c("Method"="drug.code")) %>%
    dplyr::group_by(Method, test, species.code, Ctrl, Extraction) %>%
    dplyr::mutate(MedianDifftest=median(DiffCtrlSample,na.rm=T)) %>%
    dplyr::group_by(Method, species.code, Ctrl, Extraction) %>%
    dplyr::mutate(MedianDiff=median(DiffCtrlSample,na.rm=T),  SEMDiff=sd(DiffCtrlSample,na.rm=T)/sqrt(sum(!is.na(DiffCtrlSample)))) %>%
    data.frame()
  
  data.plot = data.assaylong %>% dplyr::mutate(MedianDiff=MedianDiff*100, SEMDiff=SEMDiff*100)
  data.plot.Super = data.plot %>% dplyr::filter(Ctrl!="zero" & Condition!="med" & Extraction=="Supernatant")
  data.plot.Total = data.plot %>% dplyr::filter(Ctrl!="zero" & Condition!="med" & Extraction=="Total")
  data.plot.rep = data.plot.Super %>% 
    dplyr::inner_join(data.plot.Total, by=c("Method","test","species.short","Replicate","Ctrl")) %>%
    setNames(gsub("(\\.x)$", ".SN", names(.), perl=T)) %>%
    setNames(gsub("(\\.y)$", ".To", names(.), perl=T))
    
  myColors = drug_map$drug.color3
  names(myColors) = drug_map$drug.long
  
  pdf("reports/exp2depletion_scatter_total2supernatant.pdf", height=6, width=8)
  ggplot(data.plot.rep, aes(x=MedianDiff.SN,y=MedianDiff.To, color=drug.long.SN, shape=Ctrl)) +
    geom_point() +
    geom_errorbarh(aes(xmax = MedianDiff.SN+SEMDiff.SN, xmin = MedianDiff.SN-SEMDiff.SN), height = 0.05) +
    geom_errorbar(aes(ymax = MedianDiff.To+SEMDiff.To, ymin = MedianDiff.To-SEMDiff.To), width = 0.05) +
    geom_abline(slope=1, linetype="dashed", color="grey50") +
    scale_color_manual(values = myColors) +
    coord_cartesian(ylim=c(-105,105), xlim = c(-105,105)) +
    labs(x="Supernatant: Median difference to control, %", y="Total: Median difference to control, %") +
    myTheme
  dev.off()
}




























