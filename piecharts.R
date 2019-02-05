library(scales)
library(ggplot2)
library(ggrepel)
library(RSQLite)
library(scales)
library(tidyr)
source("functions.R")

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
}

cummulative_plot.enzymatic_coverage = function()
{
  top95clusters = read.table("data/screenG_tax_input_specI_clusters.tab", na.strings="", header=T, sep="\t", stringsAsFactors=F)
  
  db = dbConnect(SQLite(), dbname="data/kegg.db")
  all.ec = dbGetQuery(db, paste0(
    "SELECT DISTINCT 
       ec4 AS ec4 
     FROM 
       org2ec INNER JOIN 
       tax2org ON (tax2org.org=org2ec.org)
     WHERE 
       taxid IN  (", paste(na.omit(top95clusters$taxid), collapse=","), ")"))
  

  covered.ec = c(); bug_map = readr::read_delim("data/bug_map.tsv", "\t") %>%
    dplyr::slice(1:25) %>%
    dplyr::mutate(in_kegg=ifelse(is.na(species.keggorg), "Absent from KEGG", "Present in KEGG")) %>%
    dplyr::arrange(species.group_order, dplyr::desc(in_kegg)) %>%
    dplyr::mutate(species.short=factor(species.short, unique(species.short))) %>%
    dplyr::rowwise() %>% #(species.short, species.group_order, species.group, species.group_color) %>%
    dplyr::do((function(z){
      if(!is.na(z$species.keggorg)) {
        covered.ec.i = dbGetQuery(db, paste0("SELECT DISTINCT ec4 FROM org2ec WHERE org2ec.org IN ('", paste(z$species.keggorg, collapse="','"), "')"))$ec4
        covered.ec <<- unique(c(covered.ec, covered.ec.i))
        z$coverage = mean(all.ec$ec4 %in% covered.ec.i)
        z$cummulative = mean(all.ec$ec4 %in% covered.ec)
      } else {
        z$coverage=NA
        z$cummulative=NA
      }
      
      as.data.frame(z)
    })(.)) %>%
    data.frame() %>%
    dplyr::mutate(species.short.i=as.numeric(species.short) - 0.5) %>% 
    tidyr::fill(cummulative, .direction="down")
  
  myBugColors = unique(bug_map[,c("species.group", "species.group_color")]) %>% dplyr::filter(!is.na(species.group_color)) %>% data.frame()
  myBugColors = myBugColors$species.group_color %>% setNames(myBugColors$species.group)
  
  pdf("reports/enzyme_coverage.pdf", height=6, width=6)
  ggplot(bug_map) +
    geom_point(aes(x=species.short, y=cummulative), size=0.1) + 
    geom_step(aes(x=species.short.i, y=cummulative)) +
    geom_rect(aes(xmin=species.short.i+0.2, xmax=species.short.i+0.8, ymin=cummulative-0.01, ymax=cummulative+0.01, color=in_kegg, fill=species.group)) + 
    scale_y_continuous(labels=function(x) paste0(round(x*100),"%"), breaks=seq(0.3, 1, 0.1)) +
    scale_color_manual(values=c("Present in KEGG"="#FFFFFF00", "Absent from KEGG"="#000000")) +
    scale_fill_manual(values=myBugColors) +
    coord_cartesian(ylim=c(0.3, 1)) +
    myTheme + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    guides(fill=FALSE, color=FALSE) +
    labs(x="", y="cumulative enzyme coverage, %")
  
  dev.off()
}
