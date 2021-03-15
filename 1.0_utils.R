#=======================#
# Functions             #
# 1. summary statistics #
# 2. model perf         #
# 3. kappa statistics   #
#=======================#

rm(list = ls())

library(data.table)
library(haven)
library(dplyr)
library(survival)
library(purrr)
library(sandwich)
library(ggplot2)

# summary statistics: mean and sd for cont, count and pct for cat
summ_var <- function(dt){
  cont_mean_sd <- lapply(cont_list, function(var){
    dt_var <- dt %>% select("ffs", !!sym(var))
    
    dt_na <- dt_var %>% filter(is.na(!!sym(var)))
    na_n <- nrow(dt_na)
    
    dt_comp <- dt_var %>% filter(!is.na(!!sym(var)))
    dt_ffs1 <- dt_comp[dt_comp$"ffs" == 1, 2, drop = TRUE] 
    dt_ffs0 <- dt_comp[dt_comp$"ffs" == 0, 2, drop = TRUE] 
    
    mean_ffs1 <- mean(dt_ffs1)
    sd_ffs1 <- sd(dt_ffs1)
    mean_ffs0 <- mean(dt_ffs0)
    sd_ffs0 <- sd(dt_ffs0)
    
    t_test <- t.test(dt_ffs1, dt_ffs0)
    t.test(dt_ffs1, dt_ffs0)$p.value
    p_val <- t_test$p.value
    
    c("var" = var, "total_n" = length(dt_ffs0) + length(dt_ffs1) + na_n, 
      "na_n" = na_n, 
      "ffs0_n" = length(dt_ffs0), 
      "ffs0_mean_ct" =  mean_ffs0, "ffs0_sd_prop" = sd_ffs0,
      "ffs1_n" = length(dt_ffs1),
      "ffs1_mean_ct" =  mean_ffs1, "ffs1_sd_prop" = sd_ffs1,
      "p_val" = p_val)
  })
  
  bin_ct_prop <- lapply(bin_list, function(var){
    dt_var <- dt %>% select("ffs", !!sym(var))
    
    dt_na <- dt_var %>% filter(is.na(!!sym(var)))
    na_n <- nrow(dt_na)
    
    dt_comp <- dt_var %>% filter(!is.na(!!sym(var)))
    colnames(dt_comp) <- c("ffs", "var")
    dt_ffs1 <- dt_comp[dt_comp$"ffs" == 1, 2, drop = TRUE] 
    dt_ffs0 <- dt_comp[dt_comp$"ffs" == 0, 2, drop = TRUE] 
    
    ct_ffs1 <- sum(dt_ffs1 == 1)
    prop_ffs1 <- sum(dt_ffs1 == 1)/length(dt_ffs1)
    ct_ffs0 <- sum(dt_ffs0 == 1)
    prop_ffs0 <- sum(dt_ffs0 == 1)/length(dt_ffs0)
    
    chi_test <- chisq.test(table(dt_comp$ffs,dt_comp$var))
    p_val <- chi_test$p.value
    c("var" = var, "total_n" = length(dt_ffs0) + length(dt_ffs1) + na_n, 
      "na_n" = na_n, 
      "ffs0_n" = length(dt_ffs0), 
      "ffs0_mean_ct" =  ct_ffs0, "ffs0_sd_prop" = prop_ffs0,
      "ffs1_n" = length(dt_ffs1),
      "ffs1_mean_ct" =  ct_ffs1, "ffs1_sd_prop" = prop_ffs1,
      "p-value" = p_val)
  })
  rbind(do.call(rbind, cont_mean_sd), do.call(rbind, bin_ct_prop))
}

# assess model performance 
eval_perf <- function(dt){
  n <- nrow(dt)
  cms_y_whi_y <- sum(dt$crc_yy == "X")
  cms_n_whi_n <- sum(dt$crc_nn == "X")
  cms_n_whi_y <- sum(dt$crc_ny == "X")
  cms_y_whi_n <- sum(dt$crc_yn == "X")
  
  # sensitivity: TP/(TP + FN)
  sens <- cms_y_whi_y/(cms_y_whi_y + cms_n_whi_y)
  # specificity: TN/(TN + FP)
  spec<- cms_n_whi_n/(cms_n_whi_n + cms_y_whi_n)
  # ppv: TP/(TP + FP)
  ppv <- cms_y_whi_y/(cms_y_whi_y + cms_y_whi_n)
  # npv: TN/(TN + FN)
  npv <- cms_n_whi_n/(cms_n_whi_n + cms_n_whi_y)
  c("cms_y_whi_y" = cms_y_whi_y, "cms_y_whi_y_pct" = cms_y_whi_y/nrow(dt),
    "cms_n_whi_n" = cms_n_whi_n, "cms_n_whi_n_pct" = cms_n_whi_n/nrow(dt),
    "cms_y_whi_n" = cms_y_whi_n, "cms_y_whi_n_pct" = cms_y_whi_n/nrow(dt),
    "cms_n_whi_y" = cms_n_whi_y, "cms_n_whi_y_pct" = cms_n_whi_y/nrow(dt),
    "sens" = sens, "spec" = spec, "ppv" = ppv, "npv" = npv)
}

# calculate kappa
calc_kappa <- function(dt){
  n <- nrow(dt)
  
  cms_y_whi_y <- sum(dt$crc_yy == "X")
  cms_n_whi_n <- sum(dt$crc_nn == "X")
  cms_n_whi_y <- sum(dt$crc_ny == "X")
  cms_y_whi_n <- sum(dt$crc_yn == "X")
  
  p0 = (cms_y_whi_y + cms_n_whi_n)/n
  pyes = ((cms_y_whi_y + cms_y_whi_n)/n) * ((cms_y_whi_y + cms_n_whi_y)/n)
  pno = ((cms_n_whi_y + cms_n_whi_n)/n) * ((cms_y_whi_n + cms_n_whi_n)/n)
  pe = pyes + pno
  kappa = (p0-pe)/(1-pe)
  kappa.se <- sqrt(p0 * (1 - p0)/((1-pe)^2))/sqrt(n)
  
  c("kappa" = kappa,
    "kappa.lower" = kappa - 1.96 * kappa.se, 
    "kappa.upper" = kappa + 1.96 * kappa.se)
}




