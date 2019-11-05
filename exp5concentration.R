dir.create("reports", showWarnings=F)
dir.create("tmp", showWarnings=F)
dir.create("tmp/exp5concentration", showWarnings=F)
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(xlsx)
source("functions.R")

exp5concentration.preprocess = function()
{
  data.rep1  =  read.xlsx("data/exp5concentration/161118_dilution24h_rep1.xlsx",1,startRow = 40, endRow = 65)
  data.rep2  =  read.xlsx("data/exp5concentration/161118_dilution24h_rep2.xlsx",1,startRow = 40, endRow = 65)
  data.rep3  =  read.xlsx("data/exp5concentration/161118_dilution24h_rep3.xlsx",1,startRow = 40, endRow = 65)
  map  =  readr::read_delim("data/exp5concentration/well2species.tsv", "\t")
  data = rbind(data.rep1, data.rep2, data.rep3) %>%
    dplyr::mutate(Time=rep(0:24, 3), Replicate=rep(1:3, each=25)) %>%
    dplyr::select(-T..578) %>%
    reshape2::melt(id.vars=c("Time", "Replicate"), variable.name="Well", value.name="OD") %>%
    dplyr::inner_join(map, by="Well") %>%
    dplyr::select(Well, Species, Concentration, Replicate, Time, OD)
  readr::write_tsv(data, path="tmp/exp5concentration/data.tsv", na="")
}

exp5concentration.analyze = function()
{
  if(!file.exists("tmp/exp5concentration/data.tsv")) {
    exp5concentration.preprocess()
  }
  
  data = readr::read_delim("tmp/exp5concentration/data.tsv", "\t")
  pdf(file="reports/exp5concentration_curves.pdf", width=11, height=11)
  ggplot(data) + 
    geom_line(aes(x=Time, y=OD, color=Concentration, group=paste(Concentration, Replicate)))+
    facet_wrap(~Species)
  dev.off()
  
  #####
  data.approx = data %>%
    dplyr::group_by(Species, Concentration, Replicate) %>%
    dplyr::do((function(z){
      z.time = seq(0, max(z$Time), by=0.1)
      data.frame(ODfit=predict(loess(z, formula=OD~Time, span=0.5), z.time), Time=z.time)
    })(.)) %>%  
    data.frame()
  
  data.stats = data.approx %>% dplyr::group_by(Species, Replicate) %>%
    dplyr::summarise(
      maxOD.species=max(ODfit[Concentration==0]),
      maxODtime.species=Time[which.max(ODfit[Concentration==0])],
      minOD.species=min(ODfit[Concentration==1000])) %>%
    data.frame()
  
  data.sum = data.approx %>%
    dplyr::inner_join(data.stats, by=c("Species", "Replicate")) %>%
    dplyr::arrange(Species, Concentration, Replicate, Time) %>%
    dplyr::group_by(Species, Concentration, Replicate, maxOD.species, minOD.species) %>%
    dplyr::summarise(
      maxOD=ODfit[Time==maxODtime.species],
      maxODnorm=(maxOD-minOD.species[1])/(maxOD.species[1]-minOD.species[1])) %>%
    dplyr::group_by(Species, Concentration) %>%
    dplyr::summarise(
      n=length(Concentration),
      maxODnorm.se=sd(maxODnorm)/length(Replicate),
      maxODnorm=mean(maxODnorm),
      maxOD.se=sd(maxOD)/length(Replicate),
      maxOD=mean(maxOD)
    ) %>%
    dplyr::select(Species, Concentration, maxOD, maxOD.se, maxODnorm, maxODnorm.se, n)
  readr::write_tsv(data.sum, "reports/exp5concentration_growth.tsv", col_names=T, na="")  
  
  pdf(file="reports/exp5concentration_growth.pdf", width=11, height=11)
  ggplot(data.sum %>% dplyr::filter(grepl("thetaiotaomicron|rectale|gasseri|torques|salivarius", Species))) + 
    geom_line(aes(x=Concentration, y=maxODnorm, color=Species)) +
    geom_errorbar(aes(x=Concentration, ymin=maxODnorm-maxODnorm.se, ymax=maxODnorm+maxODnorm.se, color=Species), width=0.033) +
    scale_x_log10(breaks=c(0, 5,10,50,100,500,1000)) +
    scale_y_continuous(labels=scales::percent) +
    scale_color_manual(values=c("B. thetaiotaomicron"="#4B7FAE", "E. rectale"="#D60307", "L. gasseri"="#634FA2", "R. torques"="#FB775A", "S. salivarius"="#D4CFE0")) +  
    labs(xlab="Duloxetine concentration, μM", ylab="Relative maximum OD578, %") +
    myTheme
  dev.off()
}
