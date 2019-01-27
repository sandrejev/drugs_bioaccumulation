library(scales)
library(ggplot2)
library(ggrepel)

plot.piechart = function()
{
  blank_theme = theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  
  bug_map = readr::read_delim("data/bug_map.tsv", "\t") %>% dplyr::filter(!is.na(species.NT_ID))
  drug_map = readr::read_delim("data/drug_map.tsv", "\t")
  
  myBugColors = unique(bug_map[,c("species.group", "species.group_color")]) %>% dplyr::filter(!is.na(species.group_color)) %>% data.frame()
  myBugColors = myBugColors$species.group_color %>% setNames(myBugColors$species.group)
  
  bug_map.pie = bug_map[1:25,] %>%
    dplyr::group_by(species.group) %>%
    dplyr::summarise(n=length(species.group)) %>%
    dplyr::group_by() %>%
    dplyr::mutate(center=sum(n) - cumsum(n) + n/2)
  
  pdf("reports/bugs_piechart.pdf", width=6, height=6)
  ggplot(bug_map.pie)+
    geom_bar(aes(x=factor(1), y=n, fill=species.group), width=1, stat="identity") +
    ggrepel::geom_text_repel(aes(x=factor(1), y=center, label=paste0(species.group, " (", n, ")"))) +
    scale_fill_manual(values=myBugColors) +
    coord_polar("y") +
    blank_theme+
    theme(axis.text.y=element_blank(), axis.text.x=element_blank())
  dev.off()
    
  
  atc.codes = drug_map %>% 
    dplyr::group_by(drug.atc) %>% 
    dplyr::summarise(n=length(drug.atc)) %>% 
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::mutate(drug.atc_color=colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(atc.codes)))
  drug_map$ATC.code <- factor(drug_map$drug.atc, atc.codes)
  myGreys <- c(rep(c("#D9D9D9","#BDBDBD","#737373"),4),"#BDBDBD")
  
  myDrugColors = drug_map %>% dplyr::group_by(drug.atc) %>% dplyr::summarise(drug.atc_color=drug.atc_color[1])
  myDrugColors = setNames(myDrugColors$drug.atc_color, myDrugColors$drug.atc)
  ggplot(drug_map[!duplicated(drug_map$drug.long),], aes(x=factor(1), fill=drug.atc))+
    geom_bar(width = 1)+
    scale_fill_manual(values=myDrugColors) +
    coord_polar("y") +
    blank_theme
}

