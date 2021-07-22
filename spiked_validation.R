library(dplyr)
library(ggplot2)
library(readr)
library(ggrepel)

exp14spiked_validation = function() 
{
  data = readr::read_tsv("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/data/some_data.tsv") %>%
    dplyr::mutate(SampleType=ifelse(grepl("spiked", Sample), "Spiked metabolite measurement", "Raw metabolite measurement")) %>%
    dplyr::mutate(species.short=gsub(".*\\(([^)]+)\\).*", "\\1", Sample)) %>%
    dplyr::arrange(SampleType, Metabolite) %>%
    dplyr::mutate(facet=paste0(Metabolite, "\n", SampleType), facet=factor(facet, unique(facet)))
  levels(data$facet)
  
  colors = c("#E41A1C", "#4DAF4A")
  names(colors) = c(unique(data$SampleType))
  
  data.summit = data %>%
    dplyr::group_by(facet) %>%
    dplyr::arrange(dplyr::desc(Count)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  for(m in unique(data$Metabolite)) {
    pdf(stringr::str_glue("C:/Users/sandrejev/Desktop/drugs_bioaccumulation/reports/exp14spiked_validation_{met}.pdf", met=m), width=12, height=10)
    p = ggplot(data %>% dplyr::filter(Metabolite==m), aes(x=RetentionTime, y=Count, color=SampleType, fill=SampleType)) +
      geom_line(size=0.1) +
      geom_area(alpha=0.5) +
      ggrepel::geom_text_repel(aes(x=RetentionTime+2, label=paste0(round(RetentionTime, 3), "\n", round(Count, 2))), data=data.summit %>% dplyr::filter(Metabolite==m)) +
      facet_wrap(~facet, scales="free", nrow=2) +
      scale_x_continuous(breaks=seq(0, 24, 2)) +
      scale_fill_manual(values=colors) +
      scale_color_manual(values=colors) +
      theme_classic(base_size=18) +
      guides(color=F) +
      theme(panel.border=element_blank(), strip.background=element_blank(), legend.position="bottom", legend.box="horizontal")
    print(p)
    dev.off()
  }
}
