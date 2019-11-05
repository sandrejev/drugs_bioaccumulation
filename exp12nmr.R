library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggbeeswarm)

data_long = readr::read_tsv("data/exp12nmr/integrated.tsv") %>%
  dplyr::mutate(Bug=factor(Bug, c("No bug", "S. salivarius", "IAI1", "ED1a", "C. saccharolyticum"))) %>%
  dplyr::mutate(Drug=factor(Drug, c("Duloxetine", "No drug"))) 
data_long.sum = data_long %>%
  dplyr::group_by(File, SampleName, Bug, Drug) %>%
  dplyr::summarise(Integral=mean(Integral))
readr::write_delim(data_long.sum, path="reports/exp12nmr_duloxetine.tsv", delim="\t", na="")

pdf("reports/exp12nmr_duloxetine.pdf", width=8, height=8)
ggplot(data_long.sum) +
  geom_boxplot(aes(x=Bug, y=Integral, fill=Drug), outlier.size=0, position=position_dodge2(0.75, preserve="single")) + 
  geom_point(aes(x=Bug, y=Integral, fill=Drug), position=position_dodge(0.75, preserve="total")) + 
  scale_fill_manual(values=c("Duloxetine"="#ef3628", "No drug"="#cacac8")) +  
  scale_y_continuous(breaks=seq(0, 1, 0.2)) +
  labs(x="", y="Normalized duloxetine peak\narea") +
  guides(color=F) +
  theme_classic(base_size=16) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    legend.justification=c(1,1), 
    legend.position=c(0.95, 0.95), 
    aspect.ratio=1)
dev.off()

