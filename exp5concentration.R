# data.rep1  =  read.xlsx("data/exp5concentration/161118_dilution24h_rep1.xlsx",1,startRow = 40, endRow = 65)
# data.rep2  =  read.xlsx("data/exp5concentration/161118_dilution24h_rep2.xlsx",1,startRow = 40, endRow = 65)
# data.rep3  =  read.xlsx("data/exp5concentration/161118_dilution24h_rep3.xlsx",1,startRow = 40, endRow = 65)
# map  =  readr::read_delim("data/exp5concentration/well2species.tsv", "\t")
# data = rbind(data.rep1, data.rep2, data.rep3) %>%
#   dplyr::mutate(Time=rep(0:24, 3), Replicate=rep(1:3, each=25)) %>%
#   dplyr::select(-T..578) %>%
#   reshape2::melt(id.vars=c("Time", "Replicate"), variable.name="Well", value.name="OD") %>%
#   dplyr::inner_join(map, by="Well") %>%
#   dplyr::select(Well, Species, Concentration, Replicate, Time, OD)
# readr::write_tsv(data, path="data/exp5concentration/data.tsv", na="")

data = readr::read_delim("data/exp5concentration/data.tsv", "\t")
ggplot(data) + 
  geom_line(aes(x=Time, y=OD, color=Concentration, group=paste(Concentration, Replicate)))+
  facet_wrap(~Species)

data_sum = data %>% 
  dplyr::group_by(Species, Replicate) %>%
  dplyr::mutate(maxOD=max(OD), scaledOD=OD/maxOD) %>%
  dplyr::group_by(Species, Concentration, Replicate) %>%
  dplyr::summarise(scaledOD=max(scaledOD)) %>%
  dplyr::group_by(Species, Concentration) %>%
  dplyr::summarise(n=length(scaledOD), scaledOD.se=mean(scaledOD)/length(scaledOD), scaledOD=mean(scaledOD))

ggplot(data_sum) + 
  geom_line(aes(x=Concentration, y=scaledOD, color=Species)) +
  scale_x_log10(breaks=c(0, 5,10,50,100,500,1000)) 


data_sum = data %>% 
  dplyr::arrange(Species, Concentration, Replicate, Time) %>%
  dplyr::group_by(Species, Replicate) %>%
  dplyr::mutate(
    maxOD=max(OD[Concentration==0]),
    minOD=min(OD[Concentration==1000]),
    halfmax=maxOD/2,
    halfmax_time=Time[which(OD[Concentration==0]>halfmax)[1]],
    I50=(halfmax-minOD)/2+minOD) %>%
  dplyr::group_by(Species, Concentration, Replicate) %>%
  dplyr::mutate(ODhmt=OD[which(Time==halfmax_time[1])]) 

data_sum_sum = data_sum %>%
  dplyr::mutate(Time=paste0("OD", Time)) %>%
  dplyr::filter(Time %in% c("OD6", "OD8", "OD12")) %>%
  reshape2::dcast(Replicate+Species+Concentration+maxOD+minOD+halfmax+halfmax_time+ODhmt+I50~Time, value.var="OD") %>%
  dplyr::group_by(Species, Concentration) %>%
  dplyr::summarise(n=length(ODhmt), ODhmt.se=mean(ODhmt)/length(ODhmt), ODhmt=mean(ODhmt))

ggplot(data_sum_sum %>% dplyr::filter(!(Species %in% c("empty","C. saccharolyticum"))), aes(x=Concentration, y=ODhmt, color=Species)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=ODhmt-ODhmt.se, ymax=ODhmt+ODhmt.se)) +
  scale_x_log10(breaks=c(5,10,50,100,500,1000)) +
  theme_minimal(base_size = 18) +
  labs(x="Duloxetine concentration in µM", y="OD578 at time of half the maximum OD")+
  theme(panel.background = element_rect(color = "grey90"), 
        panel.grid.minor = element_line(colour = "grey95"), 
        panel.grid.major = element_line(colour = "grey50"),
        axis.ticks.length = unit(0, "lines"))









data.loess = data %>%
  dplyr::group_by(Species, Concentration) %>%
  dplyr::do((function(z) {
    z.time = sort(unique(z$Time))
    z.loess = loess(z, formula=OD~Time, span=0.5)
    data.frame(Time=z.time, ODfit=predict(z.loess, ret$Time))
  })(.)) %>% 
  dplyr::arrange(Species, Concentration, Time) %>%
  dplyr::group_by(Species) %>%
  dplyr::mutate(
    maxOD=max(ODfit[Concentration==0]),
    minOD=min(ODfit[Concentration==1000]),
    halfmax=maxOD/2,
    halfmax_time=Time[which(ODfit[Concentration==0]>halfmax)[1]],
    I50=(halfmax-minOD)/2+minOD) %>%
  dplyr::group_by(Species, Concentration) %>%
  dplyr::mutate(ODhmt=ODfit[which(Time==halfmax_time[1])]) 

data.loess_sum = data.loess %>%
  dplyr::mutate(Time=paste0("OD", Time)) %>%
  dplyr::filter(Time %in% c("OD6", "OD8", "OD12")) %>%
  reshape2::dcast(Species+Concentration+maxOD+minOD+halfmax+halfmax_time+ODhmt+I50~Time, value.var="ODfit")

ggplot(data.loess_sum %>% dplyr::filter(!(Species %in% c("empty","C. saccharolyticum"))), aes(x=Concentration, y=ODhmt)) +
  geom_point() +
  geom_smooth(span=0.5) +
  #geom_hline(aes(yintercept = I50), linetype="dashed") +
  scale_x_log10(breaks=c(5,10,50,100,500,1000)) +
  facet_wrap(~Species, scales = "free") +
  theme_minimal(base_size = 18) +
  labs(x="Duloxetine concentration in µM", y="OD578 at time of half the maximum OD")+
  theme(panel.background = element_rect(color = "grey90"), 
        panel.grid.minor = element_line(colour = "grey95"), 
        panel.grid.major = element_line(colour = "grey50")) +
  theme(axis.ticks.length = unit(0, "lines")) 










data.fitted  =  dlply(data, .(Species,Concentration), loess, formula=OD~Time, span=0.5)
data.fitted  =  ldply(data.fitted, "fitted")


data.growth  =  rbind(data.rep1,data.rep2,data.rep3)
#glimpse(data.growth)

map  =  read.csv("WellToSpecies.csv",stringsAsFactors=F)
data.growth  =  data.growth[,-c(1,2)]
data.growth  =  gather(data.growth,key=Well, value=OD, A1:H12)
data.growth  =  merge(data.growth, map, all.x=T)

mods  =  dlply(data.growth, .(Species,Concentration), loess, formula=OD~time, span=0.5)
fitted  =  ldply(mods,"fitted")

fitted  =  fitted[,c(1,2,seq(3,75,3))] # This is not correct, replicates are not ordered every 3rd value but instead come one after another 3:27, 28:52, 53:77 (table(data.fitted[1, 3:27]==data.fitted[1, 28:52]))
colnames(fitted)  =  c("Species","Concentration", paste0(rep("t",25),0:24))
fitted.long  =  gather(fitted, key=time, value=ODfit, t0:t24)
fitted.long$time  =  as.double(str_split_fixed(fitted.long$time,"t",2)[,2])


cbind(fitted[1:2], do.call(cbind, fitted[-(1:2)]))
