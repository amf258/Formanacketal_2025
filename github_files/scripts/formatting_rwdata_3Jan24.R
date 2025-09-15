library(dplR)
library(dplyr)
library(tidyr)
library(ggplot2)

# Housekeeping ------------------------------------------------------------

rm(list=ls())

# setwd("/Users/Alicia1/Documents/LiveDead/01_final_top_species")

rw<-read.csv("./data/rw/RingWidths.csv", header=T, sep=",")

## Find species with more than 100 trees for live and dead combined ----
species_sum <- rw %>% 
  group_by(Species) %>% 
  summarise(trees=length(unique(TreeNum))) %>% 
  filter(trees>100)

## species codes
species<-species_sum$Species

## species abrv names
spc<-c("abal", "abla", "auch", "piab", "pien", "pisy",
      "potr", "psme", "quru")

# looking at species with more than 100 trees -----------------------------

# read in ring widths for individuals from that site

rw <- rw %>%
  filter(Species %in% species) %>%
  select(SiteNum, Year, TreeNum, Width, YearNum, Age, Species)
new_colnames = c("SiteNum", "Year", "CoreID", "Width", "yrID", "Age", "Species")
colnames(rw)<- new_colnames

sites1 <- unique(rw$SiteNum)


## temperature data
tmp<-read.csv("./data/climate/tmp_2017_11_28.csv")

## finding long-term average of monthly temp at each site 
avg_tmp <- tmp %>% 
  filter(site_num %in% sites1) %>% 
  group_by(site_num, month) %>%
  summarise(avg_tmp=mean(tmp))

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# summer=2 and winter=1
season <- function(x) {
  ifelse(x > 0, 2, 1)}

## scale tmp data by site
site_season_key<- avg_tmp %>% 
  group_by(site_num) %>% 
  mutate(scale_temp = scale_this(avg_tmp)) %>% 
  mutate(season=season(scale_temp)) %>% 
  select(site_num, month, season)


# Add necessary columns ---------------------------------------------------

## Step 1: Load file that contains site names and their location (latitude/longitude)
site.dat <- read.csv("./data/rw/MortalitySites.csv", header = T, stringsAsFactors = F) 

## Step 2: Only keep the columns that we want to add to the rw file
site.dat <- 
  site.dat %>%
  filter(SiteNum %in% sites1) %>% 
  # select the columns that we want to keep
  select(SiteNum, site, longitude, latitude)

## Step 3: Add this data to the rw file
rw_sites <- 
  rw %>% 
  # we will use left join to add the relevant columns from site.dat to rw
  left_join(site.dat, by="SiteNum")

## Step 4: Add a column for status (live/dead)
tree.table<-
  read.csv("./data/rw/TreeTable2.csv", header = T, stringsAsFactors = F) %>% 
  select(TreeNum, site, status) 

rw_final<-
  rw_sites %>% 
  left_join(tree.table, by = c("CoreID"= "TreeNum" ,"site" = "site")) 


# Add climate data (new) -------------------------------------------------
tmp<-read.csv("./data/climate/tmp_2017_11_28.csv")
precip<-read.csv("./data/climate/pre_2017_11_28.csv")
# pet<-read.csv("/data/environ_data/pet_2017_11_28.csv")
# frs<-read.csv("/data/environ_data/frs_2017_11_28.csv")

temperature <-
  tmp %>%
  filter(site_num %in% sites1) %>%
  left_join(site_season_key, by= c('site_num', 'month'))

precipitation <- precip %>% 
  filter(site_num %in% sites1) %>%
  left_join(site_season_key, by= c('site_num', 'month'))


# Create a "water year" based on site_season_key to connect seasons that span different years
foo <- temperature %>% 
  as_tibble() %>% 
  select(site_num, year, season, month) %>% 
  group_by(site_num) %>% 
  mutate(season2=season-lag(season)) %>% 
  filter(season2!=0) %>% 
  mutate(t.state = ifelse(season2==1, "break1", "break2")) %>% 
  ungroup() %>% 
  select(site_num, t.state, t.month=month) %>% 
  distinct() %>% 
  pivot_wider(names_from = t.state,
              values_from = t.month) %>% 
  right_join(temperature, by="site_num") %>% 
  group_by(site_num) %>% 
  mutate(water_yr = ifelse(month<break2, year, year+1)) %>% 
  select(-break1, -break2)


summer_tmp <- foo %>% 
  filter(season=='2') %>% 
  group_by(year, site_num) %>%
  summarise(summer_tmp = mean(tmp))

winter_tmp <- foo %>% 
  filter(season=='1') %>% 
  group_by(year, site_num) %>%
  summarise(winter_tmp = mean(tmp))

foo2 <-precipitation %>% 
  as_tibble() %>% 
  select(site_num, year, season, month) %>% 
  group_by(site_num) %>% 
  mutate(season2=season-lag(season)) %>% 
  filter(season2!=0) %>% 
  mutate(t.state = ifelse(season2==1, "break1", "break2")) %>% 
  ungroup() %>% 
  select(site_num, t.state, t.month=month) %>% 
  distinct() %>% 
    pivot_wider(names_from = t.state,
                values_from = t.month) %>% 
  right_join(precipitation, by="site_num") %>% 
  group_by(site_num) %>% 
  mutate(water_yr = ifelse(month<break2, year, year+1))
    
summer_ppt <- foo2 %>% 
  filter(season=='2') %>% 
  group_by(year, site_num) %>%
  summarise(summer_ppt = sum(pre))

climdata<- foo2 %>% 
  filter(season=='1') %>% 
  group_by(year, site_num) %>%
  summarise(winter_ppt = sum(pre)) %>% 
  left_join(summer_ppt, by=c("site_num","year")) %>% 
  left_join(winter_tmp, by=c("site_num","year")) %>% 
  left_join(summer_tmp, by=c("site_num","year")) 

save(climdata, file = "./data/climate/climdata.Rdata")
save(rw_final, file = "./data/rw/rw_final.Rdata")
load(file = "./data/rw/rw_final.Rdata")
aspen <- rw_final %>% 
  filter(Species==17)
saveRDS(aspen, "./data/rw/aspen_rw.RDS")
write.csv(aspen, file="./data/rw/aspen_rw.csv", row.names = FALSE)
## Just stats to double check
status_sum <- rw_final %>% 
  group_by(Species) %>% 
  select(CoreID,status) %>% 
  distinct() %>%  
  summarise(live=sum(status=="LIVING"), dead=sum(status=="DEAD"))

plots_sum <- rw_final %>% 
  group_by(Species) %>% 
  select(Species,SiteNum) %>% 
  distinct() %>% 
  summarise(plots=length(unique(SiteNum)))

## making the final dataset with site data ----
## REDUNDANT
# tree.table<-
#   read.csv("./TreeTable2.csv", header = T, stringsAsFactors = F) %>% 
#   select(TreeNum, site, status) 
# 
# rw_final<-
#   rw_subset %>% 
#   left_join(tree.table, by = "TreeNum") 
# 
# ## unique site numbers for subset
# sites<- unique(rw_final$SiteNum)


# smoothed plots ----------------------------------------------------------
library(gratia)
library(dplR)
library(tidyr)
library(mgcv)
par(mfrow=c(2,1))


############################dead Abies alba

d.abies_alba <- rw_final %>% 
  filter(Species==1, status=="DEAD") 

yrID <- sort(unique(d.abies_alba$yrID))

rwl <- d.abies_alba %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(d.abies_alba$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

d.abal <- d.abies_alba %>% 
  left_join(rwi, by= c("CoreID", "yrID"))

############################live

l.abies_alba <- rw_final %>% 
  filter(Species==1, status=="LIVING") 

yrID <- sort(unique(l.abies_alba$yrID))

rwl <- l.abies_alba %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(l.abies_alba$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

l.abal <- l.abies_alba %>% 
  left_join(rwi, by= c("CoreID", "yrID"))

# GAM models 
dead <- gam(RWI ~ s(Year),
            data = d.abal, method = "REML")
live <- gam(RWI ~ s(Year),
            data = l.abal, method = "REML")

# create and object that contains the info to compare smooths
comp <- compare_smooths(dead, live)

# Save the plot as a PNG file
png("abal_rwi_smooth.png", width = 700, height = 600) 
# plot
draw(comp)

dev.off()

# ggsave("abal_rwi_smooth.png", plot = comp, width = 10, height = 6, units = "in")
############################dead Picea engelmannii

d.picea_eng <- rw_final %>% 
  filter(Species==11, status=="DEAD")

yrID <- sort(unique(d.picea_eng$yrID))

rwl <- d.picea_eng %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(d.picea_eng$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

d.pien <- d.picea_eng %>% 
  left_join(rwi, by= c("CoreID", "yrID"))



############################live

l.picea_eng <- rw_final %>% 
  filter(Species==11, status=="LIVING")

yrID <- sort(unique(l.picea_eng$yrID))

rwl <- l.picea_eng %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(l.picea_eng$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

l.pien <- l.picea_eng %>% 
  left_join(rwi, by= c("CoreID", "yrID"))


# GAM models 
dead <- gam(RWI ~ s(Year),
            data = d.pien, method = "REML")
live <- gam(RWI ~ s(Year),
            data = l.pien, method = "REML")

# create and object that contains the info to compare smooths
comp <- compare_smooths(dead, live)

# Save the plot as a PNG file
png("pien_rwi_smooth.png", width = 700, height = 600) 
# plot
draw(comp)

dev.off()

############################dead Quercus rubra

d.queru <- rw_final %>% 
  filter(Species==23, status=="DEAD")

yrID <- sort(unique(d.queru$yrID))

rwl <- d.queru %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(d.queru$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

d.quru <- d.queru %>% 
  left_join(rwi, by= c("CoreID", "yrID"))

############################live 

l.queru <- rw_final %>% 
  filter(Species==23, status=="LIVING")

yrID <- sort(unique(l.queru$yrID))

rwl <- l.queru %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(l.queru$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

l.quru <- l.queru %>% 
  left_join(rwi, by= c("CoreID", "yrID"))



# GAM models 
dead <- gam(RWI ~ s(Year),
            data = d.quru, method = "REML")
live <- gam(RWI ~ s(Year),
            data = l.quru, method = "REML")

# create and object that contains the info to compare smooths
comp <- compare_smooths(dead, live)

# Save the plot as a PNG file
png("quru_rwi_smooth.png", width = 700, height = 600) 
# plot
draw(comp)

dev.off()
############################dead Pinus sylvestris

d.pinusy <- rw_final %>% 
  filter(Species==16, status=="DEAD")

yrID <- sort(unique(d.pinusy$yrID))

rwl <- d.pinusy %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(d.pinusy$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

d.pisy <- d.pinusy %>% 
  left_join(rwi, by= c("CoreID", "yrID"))

############################live

l.pinusy <- rw_final %>% 
  filter(Species==16, status=="LIVING")

yrID <- sort(unique(l.pinusy$yrID))

rwl <- l.pinusy %>% 
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-
  rwl[order(rwl$yrID),] %>% 
  select(-yrID) %>% 
  dplR::detrend(., method = "ModNegExp", verbose = F) %>% 
  cbind(yrID) %>%
  pivot_longer(1:length(unique(l.pinusy$CoreID)),
               names_to = "CoreID", values_to = "RWI")

rwi$CoreID <- as.integer(rwi$CoreID)

l.pisy <- l.pinusy %>% 
  left_join(rwi, by= c("CoreID", "yrID"))

# GAM models 
dead <- gam(RWI ~ s(Year),
            data = d.pisy, method = "REML")
live <- gam(RWI ~ s(Year),
            data = l.pisy, method = "REML")

# create and object that contains the info to compare smooths
comp <- compare_smooths(dead, live)

# Save the plot as a PNG file
png("pisy_rwi_smooth.png", width = 700, height = 600) 
# plot
draw(comp)

dev.off()

