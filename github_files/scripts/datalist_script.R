# Clear the workspace
rm(list=ls())

# Load required libraries
library(rjags)
library(mcmcplots)
library(tidyr)
library(dplyr)
library(jagsUI)

# Load JAGS module for computing the deviance information criterion (DIC)
load.module("dic")

# Define a standardizing function (tidyr using group and scale has issues)
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Load data files
load("./data/climate/climdata.Rdata")
load("./data/rw/rw_final.Rdata")

# Convert data frames to appropriate formats (this is the climate data and rw 
# data for the species with over 100 trees)
climdata<- data.frame(climdata)
rw_final <- data.frame(rw_final)

## species abrv names
spc_code <- sort(unique(rw_final$Species))
            
spc<-c("abal", "abla", "auch", "piab", "pien", "pisy",
       "potr", "psme", "quru")

## Living ----------------------------------------------------------------------

for (k in 1:length(spc)){
k=1
rw <- rw_final %>%
  filter(Species==spc_code[k], status=="Living")  %>%
  arrange(SiteNum) %>%
  mutate(id = as.numeric(factor(CoreID, levels = unique(CoreID)))) %>%
  mutate(SiteID = as.numeric(factor(SiteNum, levels = unique(SiteNum))))


dat.widths <- rw %>% 
  arrange(Year) %>%
  mutate(yrID = as.numeric(factor(Year, levels = unique(Year)))) %>% 
  arrange(id)

# if(min(dat.widths$Year)=="1905"){
#   print(paste0(spc[k], " is good"))
# }else{
#   print(paste0(spc[k], " is NOT good"))
# }
# }

# Get unique combinations of site and individual tree IDs
trees<- dat.widths %>% 
  select(SiteID, id) %>% 
  distinct() 

# Get unique site numbers
snum<-unique(dat.widths$SiteNum)
n=length(unique(dat.widths$SiteNum))

# Filter climate data for the selected sites
climate <- climdata %>% 
  filter(site_num %in% snum)

# Create lagged ring widths, AR1 (only prior year's ring width)
dat.widths$AR1 = dat.widths$Width

for(i in 2:dim(dat.widths)[1]){
  # If prior value in width dataset is for same core, get prior value.
  if(dat.widths$CoreID[i]==dat.widths$CoreID[i-1]) {
    dat.widths$AR1[i] = dat.widths$Width[i-1]
  }else{
    # If prior value belongs to a different core, set AR1 = NA.
    dat.widths$AR1[i] = NA
  }
}

# Scale climate variables by site
clim <-
  climate %>%
  group_by(site_num) %>%
  mutate(SCsummer_pre=scale_this(summer_ppt), SCwinter_pre=scale_this(winter_ppt),
         SCsummer_temp=scale_this(summer_tmp), SCwinter_temp=scale_this(winter_tmp)) %>% 
  ungroup()

# Scale ring widths by status and join with climate data for lm and "smart inits"
rw_clim <-
  dat.widths %>% 
  group_by(CoreID) %>%
  mutate(logAR1=log(AR1+1)) %>%
  mutate(SClogAR1=scale_this(logAR1), SCage=scale_this(Age)) %>% 
  left_join(clim, by=c("SiteNum"="site_num", "Year"="year")) %>% 
  ungroup() %>% 
  data.frame()


ppt.summer <- climate %>% 
  select(summer_ppt, site_num, year) %>% 
  mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>% 
  select(-site_num) %>%
  # filter(year>=1960) %>%
  pivot_wider(names_from = "SiteID", values_from = "summer_ppt") %>% 
  select(-year) %>% 
  data.frame()


ppt.winter <- climate %>% 
  select(winter_ppt, site_num, year) %>% 
  mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>%
  select(-site_num) %>%
  # filter(year>=1960) %>% 
  pivot_wider(names_from = "SiteID", values_from = "winter_ppt") %>% 
  select(-year) %>% 
  data.frame()

temp.summer <- climate %>% 
  select(summer_tmp, site_num, year) %>% 
  mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>% 
  select(-site_num) %>%
  # filter(year>=1960) %>% 
  pivot_wider(names_from = "SiteID", values_from = "summer_tmp") %>% 
  select(-year) %>% 
  data.frame()

temp.winter <- climate %>% 
  select(winter_tmp, site_num, year) %>% 
  mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>% 
  select(-site_num) %>%
  # filter(year>=1960) %>% 
  pivot_wider(names_from = "SiteID", values_from = "winter_tmp") %>% 
  select(-year) %>% 
  data.frame()

ppt.season.full<-array(data=NA, dim = c(dim(ppt.summer)[1],n,2))
temp.season.full<-array(data=NA, dim = c(dim(temp.summer)[1],n,2))

Nyears.climate = dim(ppt.summer)[1]
for(y in 1:Nyears.climate){
  for(s in 1:n){
    ppt.season.full[y,s,1] = (ppt.winter[y,s]-mean(ppt.winter[,s]))/sd((ppt.winter[,s]))
    ppt.season.full[y,s,2] = (ppt.summer[y,s]-mean(ppt.summer[,s]))/sd((ppt.summer[,s]))
    temp.season.full[y,s,1] = (temp.winter[y,s]-mean(temp.winter[,s]))/sd((temp.winter[,s]))
    temp.season.full[y,s,2] = (temp.summer[y,s]-mean(temp.summer[,s]))/sd((temp.summer[,s]))
  }
}
start.year = 1901
years.prior = 4
which.year = min(dat.widths$Year)
which.index = which.year - start.year + 1
index.range = (which.index - years.prior):dim(ppt.summer)[1]


ppt.season<-array(data=NA, dim = c(length(index.range),n,2))
temp.season<-array(data=NA, dim = c(length(index.range),n,2))

ppt.season <- ppt.season.full[index.range,1:n,]
temp.season <- temp.season.full[index.range,1:n,]


# Preparing data for JAGS; climate data are not standardize, standardization
# occurs in JAGS "data block":
dat.list = list(N = dim(df.datalist)[1], Nmain = 6, Nint = 15,
                v1 = c(1,1,1,1,1,2), s1 = c(1,1,1,2,2,1), v2 = c(1,2,2,2,2,2),
                s2 = c(2,1,2,1,2,2), Nlags = 5,
                year.start = 5, 
                Nyears.climate = dim(ppt.season)[1], 
                Ncores = max(rw_clim$id), 
                Nyears = max(rw_clim$yrID),
                AR1 = rw_clim$SClogAR1,
                Age = rw_clim$SCage,
                LogWidth = log(rw_clim$Width + 1),
                CoreID = rw_clim$id, Year = rw_clim$yrID,
                ppt.season = ppt.season,
                temp.season = temp.season,
                Nsites = n, siteID = trees$SiteID,
                site = rw_clim$SiteID)

save(dat.list, file = paste0("./data/data.list/", "L_", spc[k], "_datalist.Rdata"))
save(rw_clim, file = paste0("./data/rw/", "L_", spc[k], "_rw_clim.Rdata"))

}

## Dead ------------------------------------------------------------------------

for (k in 1:length(spc)){
  # k=8
  rw <- rw_final %>%
    filter(Species==spc_code[k], status=="DEAD")  %>%
    arrange(SiteNum) %>%
    mutate(id = as.numeric(factor(CoreID, levels = unique(CoreID)))) %>%
    mutate(SiteID = as.numeric(factor(SiteNum, levels = unique(SiteNum))))
  
  
  dat.widths <- rw %>% 
    arrange(Year) %>%
    mutate(yrID = as.numeric(factor(Year, levels = unique(Year)))) %>% 
    arrange(id)
  
  # if(min(dat.widths$Year)=="1905"){
  #   print(paste0(spc[k], " is good"))
  # }else{
  #   print(paste0(spc[k], " is NOT good"))
  # }
  # }
  
  # Get unique combinations of site and individual tree IDs
  trees<- dat.widths %>% 
    select(SiteID, id) %>% 
    distinct() 
  
  # Get unique site numbers
  snum<-unique(dat.widths$SiteNum)
  n=length(unique(dat.widths$SiteNum))
  
  # Filter climate data for the selected sites
  climate <- climdata %>% 
    filter(site_num %in% snum)
  
  # Create lagged ring widths, AR1 (only prior year's ring width)
  dat.widths$AR1 = dat.widths$Width
  
  for(i in 2:dim(dat.widths)[1]){
    # If prior value in width dataset is for same core, get prior value.
    if(dat.widths$CoreID[i]==dat.widths$CoreID[i-1]) {
      dat.widths$AR1[i] = dat.widths$Width[i-1]
    }else{
      # If prior value belongs to a different core, set AR1 = NA.
      dat.widths$AR1[i] = NA
    }
  }
  
  # Scale climate variables by site
  clim <-
    climate %>%
    group_by(site_num) %>%
    mutate(SCsummer_pre=scale_this(summer_ppt), SCwinter_pre=scale_this(winter_ppt),
           SCsummer_temp=scale_this(summer_tmp), SCwinter_temp=scale_this(winter_tmp)) %>% 
    ungroup()
  
  # Scale ring widths by status and join with climate data for lm and "smart inits"
  rw_clim <-
    dat.widths %>% 
    group_by(CoreID) %>%
    mutate(logAR1=log(AR1+1)) %>%
    mutate(SClogAR1=scale_this(logAR1), SCage=scale_this(Age)) %>% 
    left_join(clim, by=c("SiteNum"="site_num", "Year"="year")) %>% 
    ungroup() %>% 
    data.frame()
  
  
  ppt.summer <- climate %>% 
    select(summer_ppt, site_num, year) %>% 
    mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>% 
    select(-site_num) %>%
    # filter(year>=1960) %>%
    pivot_wider(names_from = "SiteID", values_from = "summer_ppt") %>% 
    select(-year) %>% 
    data.frame()
  
  
  ppt.winter <- climate %>% 
    select(winter_ppt, site_num, year) %>% 
    mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>%
    select(-site_num) %>%
    # filter(year>=1960) %>% 
    pivot_wider(names_from = "SiteID", values_from = "winter_ppt") %>% 
    select(-year) %>% 
    data.frame()
  
  temp.summer <- climate %>% 
    select(summer_tmp, site_num, year) %>% 
    mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>% 
    select(-site_num) %>%
    # filter(year>=1960) %>% 
    pivot_wider(names_from = "SiteID", values_from = "summer_tmp") %>% 
    select(-year) %>% 
    data.frame()
  
  temp.winter <- climate %>% 
    select(winter_tmp, site_num, year) %>% 
    mutate(SiteID = as.numeric(factor(site_num, levels = unique(site_num)))) %>% 
    select(-site_num) %>%
    # filter(year>=1960) %>% 
    pivot_wider(names_from = "SiteID", values_from = "winter_tmp") %>% 
    select(-year) %>% 
    data.frame()
  
  ppt.season.full<-array(data=NA, dim = c(dim(ppt.summer)[1],n,2))
  temp.season.full<-array(data=NA, dim = c(dim(temp.summer)[1],n,2))
  
  Nyears.climate = dim(ppt.summer)[1]
  for(y in 1:Nyears.climate){
    for(s in 1:n){
      ppt.season.full[y,s,1] = (ppt.winter[y,s]-mean(ppt.winter[,s]))/sd((ppt.winter[,s]))
      ppt.season.full[y,s,2] = (ppt.summer[y,s]-mean(ppt.summer[,s]))/sd((ppt.summer[,s]))
      temp.season.full[y,s,1] = (temp.winter[y,s]-mean(temp.winter[,s]))/sd((temp.winter[,s]))
      temp.season.full[y,s,2] = (temp.summer[y,s]-mean(temp.summer[,s]))/sd((temp.summer[,s]))
    }
  }
  start.year = 1901
  years.prior = 4
  which.year = min(dat.widths$Year)
  which.index = which.year - start.year + 1
  index.range = (which.index - years.prior):dim(ppt.summer)[1]
  
  
  ppt.season<-array(data=NA, dim = c(length(index.range),n,2))
  temp.season<-array(data=NA, dim = c(length(index.range),n,2))
  
  ppt.season <- ppt.season.full[index.range,1:n,]
  temp.season <- temp.season.full[index.range,1:n,]
  
  
  # Preparing data for JAGS; climate data are not standardize, standardization
  # occurs in JAGS "data block":
  dat.list = list(Nring = dim(rw_clim)[1], Nvars = 2, Nseasons = 2, Ninteractions = 6,
                  v1 = c(1,1,1,1,1,2), s1 = c(1,1,1,2,2,1), v2 = c(1,2,2,2,2,2),
                  s2 = c(2,1,2,1,2,2), Nlags = 5,
                  year.start = 5, 
                  Nyears.climate = dim(ppt.season)[1], 
                  Ncores = max(rw_clim$id), 
                  Nyears = max(rw_clim$yrID),
                  AR1 = rw_clim$SClogAR1,
                  Age = rw_clim$SCage,
                  LogWidth = log(rw_clim$Width + 1),
                  CoreID = rw_clim$id, Year = rw_clim$yrID,
                  ppt.season = ppt.season,
                  temp.season = temp.season,
                  Nsites = n, siteID = trees$SiteID,
                  site = rw_clim$SiteID)
  
  save(dat.list, file = paste0("./data/data.list/", "D_", spc[k], "_datalist.Rdata"))
  save(rw_clim, file = paste0("./data/rw/", "D_", spc[k], "_rw_clim.Rdata"))
  
}
