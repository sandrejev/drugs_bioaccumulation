library(ggplot2)
library(readr)
library(dplyr)

exp7spentmedia.Erectale_growth = function()
{
  bug_map = readr::read_delim("data/bug_map.tsv", "\t")
  data = readr::read_delim("data/exp7spentmedia_Erectale_growth/data.tsv", "\t", col_types=readr::cols(Annotation=readr::col_character()))  %>%
    dplyr::group_by(Plate, Well) %>% 
    dplyr::mutate(OD_corrected=OD578-OD578[Hour==0]) %>%
    dplyr::left_join(bug_map, by="species.code") %>%
    dplyr::arrange(match(species.code, c("Ctrl", "Bt", "Ss")), match(Condition, c("fresh", "conditioned", "spent"))) %>%
    dplyr::mutate(species.short=ifelse(is.na(species.short), "Control", species.short), Group=paste0(species.short, "\n", Condition))
  data$Group = factor(data$Group, unique(data$Group))
  
  pdf(file="reports/exp7spentmedia_Erectale_growth.pdf", width=9, height=6)
  data.24 = data %>% 
    dplyr::filter(Hour==24 & Concentration %in% c(0, 50) & BugsGrow=="Er") %>% 
    dplyr::select(BugsGrow, Group, Concentration, Plate, Replicate, OD_corrected, species.short)
  readr::write_tsv(data.24, file="reports/exp7spentmedia_Erectale_growth.tsv", col_names=T, na="") 
  ggplot(data.24, aes(y=OD_corrected, x=Group, fill=as.factor(Concentration))) +
    geom_boxplot() +
    geom_point(position=position_jitterdodge()) +
    scale_fill_manual(values=c("grey","firebrick4")) +
    labs(x="", y="Baseline corrected OD578",fill="Duloxetine \nconcentration") +
    facet_wrap(~BugsGrow) +
    myTheme
  dev.off()
}
