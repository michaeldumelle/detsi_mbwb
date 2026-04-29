library(tidyverse)
library(here)
library(sf)
library(ranger)
library(spmodel)

##################
### Moreton Bay
##################
dat <- read_sf(here("data", "MBsub_03-31-2026.gpkg"))

#########
### Spatial Linear Model
#########
mb_tc_splm_form <- total.cover ~ sed.type + depth.range + temp.90d.mn + turb.90d.mn +
  turb.fld.mn:yr.flood.recip2


# no spatial covariance
mb_tc_lm_mod <- splm(
  formula = mb_tc_splm_form,
  data = dat,
  spcov_type = "none"
)
saveRDS(mb_tc_lm_mod, file = here("models", "mb_tc_lm_mod.rds"))
# leave one out
mb_tc_lm_loocv <- loocv(mb_tc_lm_mod, local = FALSE)
saveRDS(mb_tc_lm_loocv, file = here("models", "mb_tc_lm_loocv.rds"))

# spatial covariance
set.seed(0)
mb_tc_splm_mod <- splm(
  formula = mb_tc_splm_form,
  data = dat,
  spcov_type = "exponential",
  local = list(size = 500),
  random = ~ ehmp.subzone
)
saveRDS(mb_tc_splm_mod, file = here("models", "mb_tc_splm_mod.rds"))
# leave one out
mb_tc_splm_loocv <- loocv(mb_tc_splm_mod, local = FALSE)
saveRDS(mb_tc_splm_loocv, file = here("models", "mb_tc_splm_loocv.rds"))

# block kriging
get_bk_grid <- function(rvar, gvar1, gvar2, data, model) {
  
  if (missing(gvar2)) {
    data$all <- "all"
    gvar2 <- "all"
  }
  
  pred_grid <- data %>%
    st_drop_geometry() %>%
    group_by(!!sym(gvar1), !!sym(gvar2)) %>% # alternatively, .data[[gvar1]]
    summarize(n = n(), raw_mean = mean(!!sym(rvar))) %>%
    arrange(!!sym(gvar1), !!sym(gvar2)) %>%
    ungroup()
  
  pred_list <- map2(pred_grid[[gvar1]], pred_grid[[gvar2]], function(x, y) {
    newdat <- data %>%
      filter(!!sym(gvar1) == x, !!sym(gvar2) == y)
    predict(model, newdat, block = TRUE, local = TRUE, interval = "prediction")
  })
  
  pred_grid <- pred_grid %>%
    mutate(fit = map_dbl(pred_list, \(x) x[, "fit"])) %>%
    mutate(fit = if_else(fit < 0, 0, fit)) %>%
    mutate(lwr = map_dbl(pred_list, \(x) x[, "lwr"])) %>%
    mutate(lwr = if_else(lwr < 0, 0, lwr)) %>%
    mutate(upr = map_dbl(pred_list, \(x) x[, "upr"])) %>%
    mutate(upr = if_else(upr > 100, 100, upr))
}

bk_mb_ehmp <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "ehmp.subzone",
  gvar2 = "fy",
  data = dat,
  model = mb_tc_splm_mod
)
saveRDS(bk_mb_ehmp, file = here("models", "bk_mb_ehmp.rds"))

bk_mb_zone <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "zone.type",
  gvar2 = "fy",
  data = dat,
  model = mb_tc_splm_mod
)
saveRDS(bk_mb_zone, file = here("models", "bk_mb_zone.rds"))

bk_mb_tidal <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "depth.range",
  gvar2 = "fy",
  data = dat,
  model = mb_tc_splm_mod
)
saveRDS(bk_mb_tidal, file = here("models", "bk_mb_tidal.rds"))

bk_mb_all <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "fy",
  data = dat,
  model = mb_tc_splm_mod
)
saveRDS(bk_mb_all, file = here("models", "bk_mb_all.rds"))

#########
### Spatial Generalized Linear Model
#########

# beta 
tc <- dat$total.cover
y <- tc/100
y <- pmax(0.01, pmin(0.99, y))
dat$total.cover <- y
set.seed(0)
mb_tc_spglm_beta_mod <- spglm(
  formula = mb_tc_splm_form,
  data = dat,
  family = "beta",
  spcov_type = "exponential",
  local = TRUE,
  random = ~ ehmp.subzone
)
saveRDS(mb_tc_spglm_beta_mod, file = here("models", "mb_tc_spglm_beta_mod.rds"))
# leave one out
mb_tc_spglm_beta_loocv <- loocv(mb_tc_spglm_beta_mod, local = FALSE)
saveRDS(mb_tc_spglm_beta_loocv, file = here("models", "mb_tc_spglm_beta_loocv.rds"))

# Gamma
y <- pmax(0.1, tc)
dat$total.cover <- y
set.seed(0)
mb_tc_spglm_gamma_mod <- spglm(
  formula = mb_tc_splm_form,
  data = dat,
  family = "Gamma",
  spcov_type = "exponential",
  local = TRUE,
  random = ~ ehmp.subzone
)
saveRDS(mb_tc_spglm_gamma_mod, file = here("models", "mb_tc_spglm_gamma_mod.rds"))
mb_tc_spglm_gamma_loocv <- loocv(mb_tc_spglm_gamma_mod, local = FALSE)
# leave one out
saveRDS(mb_tc_spglm_gamma_loocv, file = here("models", "mb_tc_spglm_gamma_loocv.rds"))
# reset total cover variable for use with formula
dat$total.cover <- tc

#########
### Random Forest
#########
mb_tc_rf_form <- total.cover ~ sed.type + depth.range + temp.90d.mn + turb.90d.mn +
  turb.fld.mn + yr.flood.recip2 + ehmp.subzone
set.seed(0)
mb_tc_rf_mod <- ranger(
  formula = mb_tc_rf_form,
  data = st_drop_geometry(dat),
  importance = "impurity",
  keep.inbag  = TRUE
)
saveRDS(mb_tc_rf_mod, file = here("models", "mb_tc_rf_mod.rds"))
# out of bag (akin to leave one out)
oob_errors <- dat$total.cover - mb_tc_rf_mod$predictions
mb_tc_rf_oob <- tibble(
  bias = mean(oob_errors),
  MSPE = mean(oob_errors^2),
  RMSPE = sqrt(MSPE),
  cor2 = cor(dat$total.cover, mb_tc_rf_mod$predictions)^2
)
saveRDS(mb_tc_rf_oob, file = here("models", "mb_tc_rf_oob.rds"))

##################
### Wide Bay
##################
dat <- read_sf(here("data", "WBsub_03-31-2026.gpkg"))

#########
### Spatial Linear Model
#########
wb_tc_splm_form <- total.cover ~ depth.range + sed.type + temp.90d.mn +
  turb.90d.mn + turb.fld.mn:yr.flood1.recip2


# no spatial covariance
wb_tc_lm_mod <- splm(
  formula = wb_tc_splm_form,
  data = dat,
  spcov_type = "none"
)
saveRDS(wb_tc_lm_mod, file = here("models", "wb_tc_lm_mod.rds"))
# leave one out
wb_tc_lm_loocv <- loocv(wb_tc_lm_mod, local = FALSE)
saveRDS(wb_tc_lm_loocv, file = here("models", "wb_tc_lm_loocv.rds"))

# spatial covariance
set.seed(0)
wb_tc_splm_mod <- splm(
  formula = wb_tc_splm_form,
  data = dat,
  spcov_type = "exponential",
  local = list(size = 500),
  random = ~ zone.class
)
saveRDS(wb_tc_splm_mod, file = here("models", "wb_tc_splm_mod.rds"))
# leave one out
wb_tc_splm_loocv <- loocv(wb_tc_splm_mod, local = FALSE)
saveRDS(wb_tc_splm_loocv, file = here("models", "wb_tc_splm_loocv.rds"))

bk_wb_class <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "zone.class",
  gvar2 = "fy",
  data = dat,
  model = wb_tc_splm_mod
)
saveRDS(bk_wb_class, file = here("models", "bk_wb_class.rds"))

bk_wb_zone <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "zone.type",
  gvar2 = "fy",
  data = dat,
  model = wb_tc_splm_mod
)
saveRDS(bk_wb_zone, file = here("models", "bk_wb_zone.rds"))

bk_wb_tidal <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "depth.range",
  gvar2 = "fy",
  data = dat,
  model = wb_tc_splm_mod
)
saveRDS(bk_wb_tidal, file = here("models", "bk_wb_tidal.rds"))

bk_wb_all <- get_bk_grid(
  rvar = "total.cover",
  gvar1 = "fy",
  data = dat,
  model = wb_tc_splm_mod
)
saveRDS(bk_wb_all, file = here("models", "bk_wb_all.rds"))

#########
### Spatial Generalized Linear Model
#########

# beta 
tc <- dat$total.cover
y <- tc/100
y <- pmax(0.01, pmin(0.99, y))
dat$total.cover <- y
set.seed(0)
wb_tc_spglm_beta_mod <- spglm(
  formula = wb_tc_splm_form,
  data = dat,
  family = "beta",
  spcov_type = "exponential",
  local = TRUE,
  random = ~ zone.class
)
saveRDS(wb_tc_spglm_beta_mod, file = here("models", "wb_tc_spglm_beta_mod.rds"))
# leave one out
wb_tc_spglm_beta_loocv <- loocv(wb_tc_spglm_beta_mod, local = FALSE)
saveRDS(wb_tc_spglm_beta_loocv, file = here("models", "wb_tc_spglm_beta_loocv.rds"))


# Gamma
y <- pmax(0.1, tc)
dat$total.cover <- y
set.seed(0)
wb_tc_spglm_gamma_mod <- spglm(
  formula = wb_tc_splm_form,
  data = dat,
  family = "Gamma",
  spcov_type = "exponential",
  local = TRUE,
  random = ~ zone.class
)
saveRDS(wb_tc_spglm_gamma_mod, file = here("models", "wb_tc_spglm_gamma_mod.rds"))
wb_tc_spglm_gamma_loocv <- loocv(wb_tc_spglm_gamma_mod, local = FALSE)
# leave one out
saveRDS(wb_tc_spglm_gamma_loocv, file = here("models", "wb_tc_spglm_gamma_loocv.rds"))

# binom
y <- ifelse(tc > 0, 1, 0)
dat$total.cover <- y
set.seed(0)
wb_tc_spglm_binom_mod <- spglm(
  formula = wb_tc_splm_form,
  data = dat,
  family = "binomial",
  spcov_type = "exponential",
  local = TRUE,
  random = ~ zone.class
)
saveRDS(wb_tc_spglm_binom_mod, file = here("models", "wb_tc_spglm_binom_mod.rds"))
wb_tc_spglm_binom_loocv <- loocv(wb_tc_spglm_binom_mod, local = FALSE)
# leave one out
saveRDS(wb_tc_spglm_binom_loocv, file = here("models", "wb_tc_spglm_binom_loocv.rds"))
# reset total cover variable for use with formula
dat$total.cover <- tc

#########
### Random Forest
#########
wb_tc_rf_form <- total.cover ~ sed.type + depth.range + temp.90d.mn + turb.90d.mn +
  turb.fld.mn + yr.flood1 + zone.class
set.seed(0)
wb_tc_rf_mod <- ranger(
  formula = wb_tc_rf_form,
  data = st_drop_geometry(dat),
  importance = "impurity",
  keep.inbag  = TRUE
)
saveRDS(wb_tc_rf_mod, file = here("models", "wb_tc_rf_mod.rds"))
# out of bag (akin to leave one out)
oob_errors <- dat$total.cover - wb_tc_rf_mod$predictions
wb_tc_rf_oob <- tibble(
  bias = mean(oob_errors),
  MSPE = mean(oob_errors^2),
  RMSPE = sqrt(MSPE),
  cor2 = cor(dat$total.cover, wb_tc_rf_mod$predictions)^2
)
saveRDS(wb_tc_rf_oob, file = here("models", "wb_tc_rf_oob.rds"))

