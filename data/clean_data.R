library(tidyverse)
library(here)
library(sf)
library(readxl)

# read in raw geopackage
dat <- read_sf(here("data", "MB_03-31-2026.gpkg")) %>%
  st_transform(crs = 3113)

# create new sediment type variable
dat <- dat %>%
  mutate(sed.type = case_when(
    sed.type %in% c("M", "MR") ~ "Mud",
    sed.type %in% c("SM", "RSM") ~ "Sandy.Mud",
    sed.type %in% c("MS", "MSR") ~ "Muddy.Sand",
    sed.type %in% c("S", "G", "SR") ~ "Sand",
    sed.type %in% c("R", "RM", "RMS", "RS", "SH") ~ "Rock",
    TRUE ~ "Other"
  ))  %>%
  # create month, year, growing season variables
  mutate(month = str_sub(date, 6, 7)) %>%
  mutate(year = str_sub(date, 1, 4)) %>%
  mutate(growing = if_else(month %in% 9:12, "Growing", "Not-Growing"))

var_names <- colnames(dat)

# define variables to keep
vars <-  c("sed.type", "depth.range", "temp.90d.mn", "turb.90d.mn",
           "turb.fld.mn", "ds.flood", "sal.fld.mn", "fy", "month", "year",
           "dom.spp", "ehmp.subzone", "ehmp.zone", "zone.type",
           "site.id", "nearest.river", "dist.km", "growing", "date")

# remove data with missing eReefs variables
dat <- dat %>%
  dplyr::select(c(total.cover, all_of(vars))) %>%
  na.omit()

# add days since flood, years since flood, an after flood effect,
# a decayign flood effect, and truncate the turbidity and salinity flood 
# variables to be zero before the flood
dat <- dat %>%
  mutate(ds.flood = if_else(ds.flood < 0, 0, ds.flood)) %>%
  mutate(yr.flood = ds.flood / 365) %>%
  mutate(after.flood = if_else(ds.flood > 0, 1, 0)) %>%
  mutate(yr.flood.recip2 = if_else(yr.flood == 0, 0, 1/yr.flood^2)) %>%
  mutate(turb.fld.mn = if_else(yr.flood <= 0, 0, turb.fld.mn)) %>%
  mutate(sal.fld.mn = if_else(yr.flood <= 0, 0, sal.fld.mn))

vars <- colnames(dat)
dat_copy <- dat
dat_vars <- vars[! vars %in% c("dom.spp", "total.cover")]
dat <- dat %>%
  group_by(across(all_of(dat_vars))) %>%
  summarize(total.cover = sum(total.cover)) %>%
  ungroup()


# there are several ties here, pick one randomly
dat_spp <- dat_copy %>%
  group_by(across(all_of(dat_vars))) %>%
  slice_max(total.cover, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site.id, dom.spp) %>%
  st_drop_geometry() 

dat <- left_join(dat, dat_spp, by = "site.id")

dat_cov <- read_xlsx(here("data", "MBsub_03-31-2026_covariate_data.xlsx"), sheet = "MBsub_compiled (2)")
dat_cov <- dat_cov %>%
  select(-fid, -date, -dom.spp) %>%
  group_by(site.id) %>%
  distinct() %>%
  ungroup()

dat <- left_join(dat, dat_cov, by = "site.id")

st_write(dat, here("data", "MBsub_03-31-2026.gpkg"))


dat <- read_sf(here("data", "WB_03-31-2026.gpkg")) %>%
  st_transform(crs = 3113)


dat <- dat %>%
  mutate(across(c(dist, sal.90d.30.freq:turb.fld2.max), as.numeric))

dat <- dat %>%
  mutate(turb.fld.mn = (turb.fld1.mn + turb.fld2.mn)/2) %>%
  mutate(sal.fld.mn = (sal.fld1.mn + sal.fld2.mn)/2)

dat <- dat %>%
  mutate(fy = str_sub(date, 1, 4)) %>%
  filter(fy > 2002) %>%
  mutate(fy = str_sub(date, 1, 7)) %>%
  mutate(month = str_sub(date, 6, 7)) %>%
  mutate(year = str_sub(date, 1, 4))

dat <- dat %>%
  rename(tidal = tidal_3cat, substrate = Substrate.5cat) %>%
  mutate(fy = case_when(
    fy %in% "2022-05" ~ "2021-2022",
    fy %in% c("2022-10", "2022-11") ~ "2022-2023",
    .default = "2023-2024"
  ))

dat$yr.flood1 <- dat$ds.flood1 / 365.25
dat$ds.flood2 <- ifelse(dat$ds.flood2 < 0, 0, dat$ds.flood2)
dat$yr.flood2 <- dat$ds.flood2 / 365.25

dat <- dat %>%
  mutate(yr.flood1 = ds.flood1 / 365.25) %>%
  mutate(yr.flood2 = ds.flood2 / 365.25) %>%
  mutate(yr.flood1.recip2 = 1/yr.flood1^2) %>%
  mutate(yr.flood2.recip2 = 1/yr.flood2^2) %>%
  rename(sed.type = substrate, depth.range = tidal)

zones <- read_csv(here("data", "WB_zone.name.csv"))

dat <- dat %>% left_join(zones, by = "zone.name")

vars <-  c("depth.range", "sed.type", "seagrass.b", "seagrass.p", "temp.90d.mn", "turb.90d.mn", 
           "zone.type", "yr.flood1", "yr.flood2", "turb.fld.mn", "sal.fld.mn", "zone.class", "fy", "month", "year",
           "yr.flood1.recip2", "yr.flood2.recip2", "dom.spp", "total.cover", "uniq.id")

dat <- dat %>%
  dplyr::select(c(seagrass.b, all_of(vars))) %>%
  na.omit() 

dat2 <- read_xlsx(here("data", "WBsub_03-31-2026_covariate_data.xlsx"), sheet = "WBsub_compiled") %>%
  select(-fid, -date, -dom.spp)
dat <- left_join(dat, dat2, by = "uniq.id")
st_write(dat, here("data", "WBsub_03-31-2026.gpkg"))
