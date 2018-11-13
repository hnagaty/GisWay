library(readxl)
library(readr)
library(dplyr)
library(tidyr)


gAreas <- read_excel("D:/Optimisation/~InProgress/201806_GisFramework/gAreas.xlsx", 
                    col_types = c("skip", "text", "skip","text", "skip", "numeric", "numeric"),
                    col_names = c("Cell","NbrCell","Area","Perimeter"),
                    skip=1)

gsmHo <- read_delim("D:/Optimisation/~InProgress/201806_GisFramework/ossData/GsmHo.txt",
                    "\t", escape_double = FALSE,
                    skip=1,
                    col_names = c("Weekday","AdjType","BSC","CellName","NbrCellName","PhysicalSite","NbrPhysicalSite","Daytime","DayType","CellTraffic","HoCnt","HoSucc"),
                    col_types = cols_only(
                      Weekday = col_integer(),
                      AdjType = col_factor(levels = c("I","E")),
                      CellName = col_character(), 
                      NbrCellName = col_character(),
                      Daytime = col_factor(levels = c("Day","Night")),
                      HoCnt = col_double(),
                      HoSucc = col_double(),
                      CellTraffic= col_double(),
                      DayType = col_factor(levels = c("Workday","Friday"))),
                    na = c("NA","#EMPTY"),
                    trim_ws = TRUE)


gsmHo <- gsmHo %>%
  select(-Weekday) %>%
  group_by(CellName,NbrCellName,AdjType,Daytime,DayType) %>%
  summarise(HoCnt=sum(HoCnt),HoSucc=sum(HoSucc),CellTraffic=mean(CellTraffic)) %>%
  ungroup()

gsmHoS <- gsmHo %>%
  select(CellName,NbrCellName,HoCnt) %>%
  filter(substr(CellName,1,1) != "C",substr(NbrCellName,1,1) != "C") %>%
  group_by(CellName,NbrCellName) %>%
  summarise(HoCnt=sum(HoCnt)) %>%
  group_by(CellName) %>%
  mutate(HoPcnt = prop.table(HoCnt)) %>%
  arrange(CellName,desc(HoCnt))

write_csv(gsmHo,"GsmHoGrouped.csv")

gAreasS <- gAreas %>%
  group_by(Cell) %>%
  mutate(AreaPcnt=prop.table(Area)) %>%
  select(-Perimeter) %>%
  arrange(Cell,desc(AreaPcnt))

write_csv(gAreasS,"Areas.csv")

combined <- full_join(gAreasS,gsmHoS,by=c("Cell"="CellName","NbrCell"="NbrCellName")) %>%
  semi_join(gAreasS,by="Cell") %>%
  select(-Area,-HoCnt) %>%
  mutate_all(funs(replace(., is.na(.), 0))) %>%
  mutate(Site=substr(Cell,1,nchar(Cell)-1))

siteExample <- combined %>%
  filter(Site=="4182") %>%
  rename(Cellname=Cell) %>%
  mutate(Cell=substring(Cellname,nchar(Site)+1)) %>%
  arrange(Cellname,NbrCell) %>%
  select(-Cell,-Site)

siteHo <- siteExample %>%
  select(-AreaPcnt) %>%
  spread (NbrCell,HoPcnt) %>%
  mutate(Sector=substring(Cellname,nchar(Cellname)))
siteArea <- siteExample %>%
  select(-HoPcnt) %>%
  spread (NbrCell,AreaPcnt) %>%
  mutate(Sector=substring(Cellname,nchar(Cellname)))
