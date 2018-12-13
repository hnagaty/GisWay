# Reads MRR files
# Summarizes them
# Do some plots on idividual cells or groups


# Load libraries ----------------------------------------------------------
source("00_globalVars.R")
library(Hmisc)



# Read in the WCDMA MRR files -----------------------------------------------
wmrr <- readWmrrFiles(paths$wmrrDir,getMrrFiles(paths$wmrrDir))

qtls <- seq(0.15,0.9,0.15) # the stops
qtlsnames <- paste0("Q",qtls*100)


#Combine all instances into one row
wmrr <- wmrr %>%
  group_by(CellName,CellUserLabel,Service,Quantity) %>%
  summarise_all(sum) %>%
  ungroup()

# remove DL power from wmrr, as it's missing in many cells
# I may choose to add it later
wmrr <- wmrr %>%
  filter(!startsWith(Quantity,"DLT"))

# filter rows with all zeros
wmrr <- wmrr %>% 
  mutate(binsSum = rowSums(select(.,starts_with("Bin")), na.rm = TRUE)) %>%
  filter (binsSum !=0)

wmrrTest <- wmrr %>%
  separate(CellName, into = c("RNC", "C"), sep = "/") %>%
  select(-C) %>%
  rename(CellName = CellUserLabel) %>%
  filter(CellName == "G43481") %>%
  gather(BinC, Value, starts_with("Bin")) %>%
  filter (!is.na(Value)) %>%
  mutate(Bin = as.integer(sub("Bin","", BinC))) %>%
  select(-BinC) %>%
  arrange(Quantity,Bin)

# this one for wmrr features
wmrrFeatures <- wmrr %>% # I re-used some code from gsm mrr; the code for computing quantiles
  select(-CellName) %>%
  rename(CellName = CellUserLabel) %>%
  #filter (CellName=="GD3153A") %>%
  mutate_all(funs(replace(., is.na(.), 0))) %>%
  unite(cnt, starts_with("Bin"), sep=",") %>%
  mutate(val = paste0(seq (0,37, 1),collapse=",")) %>% # this step is not fool proof, if no. of bins is not 37
  rowwise() %>%
  mutate(q=paste0(wtd.quantile(as.numeric(strsplit(val,",",fixed=TRUE)[[1]]),
                               weights=as.numeric(strsplit(cnt,",",fixed=TRUE)[[1]]),
                               probs=qtls),
                  collapse=",")) %>%
  ungroup %>%
  select(-val,-cnt) %>%
  separate(q,sep=",",into=qtlsnames,convert=TRUE) %>%
  gather(key=qtl,value=value,-(CellName:Quantity)) %>%
  unite(ftr,Service, Quantity,qtl,sep="-") %>%
  spread(key=ftr,value = value)


# clustering
# identifier & notes for the model
data <- wmrrFeatures[,2:ncol(wmrrFeatures)]
labels <- wmrrFeatures[1]
cClass <- "WMRR"
cNotes <- "1st trial on WMRR"
cVersion <- "v0.01"
cDate <- as.Date("2018-12-13")
cList <- list(class=cClass, version=cVersion, date=cDate, notes=cNotes)
tic("Clustering Fit")
hClustFit <- hClustMrr(data, labels, cList, k=8)
toc()

# Explore other cuts
k <- NULL
h <- 400
newClusters <- tryHClust(model = hClustFit$model, labeled = hClustFit$clusters,
                         h = h, k =k)


#===========================================================================================
