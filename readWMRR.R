library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Hmisc)

wmrr <- read_tsv("../ossData/HanyWMRR.msmt",
                 col_types = "cccciiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii")
summary(wmrr)

wmrr$Service <- as.factor(wmrr$Service)
wmrr$Quantity <- as.factor(wmrr$Quantity)

wmrr <- wmrr %>% select(-(Bin17:Bin37),-CellName)
colnames(wmrr) <- c("CellUserLabel","Service","Quantity","Bin00","Bin01",         
                    "Bin02","Bin03","Bin04","Bin05","Bin06",         
                    "Bin07","Bin08","Bin09","Bin10","Bin11",        
                    "Bin12","Bin13","Bin14","Bin15","Bin16")


summary(wmrr)
str(wmrr)

wmrrGr <- wmrr %>%
  select (-CellUserLabel) %>%
  group_by(Service,Quantity) %>%
  summarize_all(funs(mean(.,na.rm=TRUE)))

wmrrTidy <- wmrrGr %>%
  gather(Bin,Value,3:19)

wmrrCheck <- wmrrTidy %>%
  select (-Bin) %>%
  group_by(Service,Quantity) %>%
  summarize_all(funs(sum(.,na.rm=TRUE)))

wmrrTidy <- wmrrTidy %>%
  filter(Service!="PS Conversational Speech EUL/HS")

rm(wmrrCheck)

wmrrSmall <- wmrrTidy %>%
  filter(Service=="Speech AMR NB 12.2")
  select(-Service)

ggplot(wmrrSmall,aes(x=Bin,y=Value)) +
  geom_col() +
  facet_grid(.~Quantity)

wmrrPerCell <- wmrr %>%
  gather(Bin,Value,4:20) %>%
  filter(Service=="Speech AMR NB 12.2") %>%
  select(-Service) %>%
  rename(Cell=CellUserLabel) %>%
  group_by(Cell,Quantity) %>%
  summarise(avg=mean(Value,na.rm=TRUE),sd=sd(Value,na.rm=TRUE),
            min=min(Value),max=max(Value),
            Q90=quantile(Value,probs=0.9))

