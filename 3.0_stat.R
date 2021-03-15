#=======================#
# Results               #
# 1. baseline all var   #
# 2. algo perf all ffs  #
# 3. kappa by baseline  #
# 4. sens vs spec       #
# 5. cox                #
#=======================#

source("2.0_pop.R")

# Table 1: Baseline Clinical Characteristics of Women in WHI
# baseline characteristics for ffs vs non-ffs
cont_list <- c("AGE", "meno_year")
bin_list <- c("large_bmi", "SMOKEVR", "ALC12DR", "TOTH", "DIABTRT", "HYPT", "HICHOLRP")

table_1 <- summ_var(pop1)
table_1

# Table 2: Summary of Algorithm Performance among all eligible patients
table_2 <- c(eval_perf(pop2), calc_kappa(pop2))
table_2

# Table 3A: Baseline Characteristics and predictive performance (N = 107,743)
# kappa stratified by baseline
var_list <- c("old_age", "long_meno", "TOTH", "large_bmi", "SMOKEVR", "ALC12DR", 
              "DIABTRT", "HYPT", "HICHOLRP")

set.seed(47)
boot_kappa <- lapply(var_list, function(var){
  dt_na <- pop3 %>% filter(is.na(!!sym(var)))
  
  grp1 <- pop3 %>% filter((!is.na(!!sym(var))) & (!!sym(var)) == 1)
  grp1_kappa <- calc_kappa(grp1)
  grp1_n = nrow(grp1)
  
  grp0 <- pop3 %>% filter((!is.na(!!sym(var))) & (!!sym(var)) == 0)
  grp0_kappa <- calc_kappa(grp0)
  grp0_n = nrow(grp0)
  
  k_diff <- grp1_kappa[1] - grp0_kappa[1]
  boot_kappa_diff <- sapply(1:1000, function(i){
    boot1 = grp1[sample(grp1_n, replace = T),]
    boot0 = grp0[sample(grp0_n, replace = T),]
    
    boot1_kappa <- calc_kappa(boot1)
    boot0_kappa <- calc_kappa(boot0)
    
    boot_kappa_diff <-  boot1_kappa[1] - boot0_kappa[1]
  })
  
  c("var" = var, "total_n" = nrow(dt_na) + grp1_n + grp0_n, "na_n" = nrow(dt_na),
    "grp1_n" = grp1_n, "grp1_kappa" = grp1_kappa, 
    "grp0_n" = grp0_n,"grp0_kappa" = grp0_kappa,
    "diff_kappa" = k_diff, 
    "diff_kappa_lower" = quantile(boot_kappa_diff, 0.025), 
    "diff_kappa_upper" = quantile(boot_kappa_diff, 0.975))
})
table_3 <- do.call(rbind, boot_kappa)


# table(pop3$whi_y)
# table(pop3[pop3$whi_y == 1,]$cms_y)
# sens <- glm(1-cms_y ~ AGE, data = pop3[pop3$whi_y == 1,], family = "poisson")
# coeftest(mod,vcov=sandwich)[5]
# 
# sens_coef <- coef(sens)[2]
# sens_se <- sqrt(diag(vcovHC(sens, type = "HC0")))[2]
# sens_low <- sens_coef - 1.96 * sens_se
# sens_up <- sens_coef + 1.96 * sens_se
# sens_p_val <- 2 * pnorm(abs(sens_coef/sens_se), lower.tail = FALSE)


# Table 3B: Predictive performance by age group relative risk ratio for old vs young
spec_age_bin <- glm(cms_y ~ old_age, data = pop3[pop3$whi_y == 0,], family = "poisson")
exp(coef(spec_age_bin))[2]

# Table 4: Relative risk ratio for univariate regression models
# relative risk ratio regression
sig_var <- c("AGE", "meno_year", "HYPT", "BMI")
table_4 <- sapply(c(cont_list, bin_list), function(i){
  
  sens_fml <- as.formula(paste0("cms_n ~ ", i))
  sens <- glm(sens_fml, data = pop3[pop3$whi_y == 1,], family = "poisson")
  sens_coef <- coef(sens)[2]
  sens_se <- sqrt(diag(vcovHC(sens, type = "HC0")))[2]
  sens_low <- sens_coef - 1.96 * sens_se
  sens_up <- sens_coef + 1.96 * sens_se
  sens_p_val <- 2 * pnorm(abs(sens_coef/sens_se), lower.tail = FALSE)
  
  spec_fml <- as.formula(paste0("cms_y ~ ", i)) 
  spec <- glm(spec_fml, data = pop3[pop3$whi_y == 0,], family = "poisson")
  spec_coef <- coef(spec)[2]
  spec_se <- sqrt(diag(vcovHC(spec, type = "HC0")))[2]
  spec_low <- spec_coef - 1.96 * spec_se
  spec_up <- spec_coef + 1.96 * spec_se
  spec_p_val <- 2 * pnorm(abs(spec_coef/spec_se), lower.tail = FALSE)
  
  c("sens_rr" = exp(sens_coef), "sens_low" = exp(sens_low), "sens_up" = exp(sens_up), 
    "sens_p_val" = sens_p_val,
    "spec_rr" = exp(spec_coef), "spec_low" = exp(sens_low), "spec_up" = exp(sens_up), 
    "spec_p_val" = spec_p_val)
})

# confirmed
spec_young <- eval_perf(pop3[pop3$old_age == 0,])["spec"]
spec_old <- eval_perf(pop3[pop3$old_age == 1,])["spec"]
(1-spec_old)/(1-spec_young)

## by age group
spec_age_all_grp <- glm(cms_y ~ as.factor(age_grp), data = pop3[pop3$whi_y == 0,], family = "poisson")
spec_age_all_grp2 <- as.data.frame(coeftest(spec_age_all_grp,vcov=sandwich)[-1,1:2])
colnames(spec_age_all_grp2) <- c("est", "robust_se")
spec_age_all_grp2$rr_low <- exp(spec_age_all_grp2$est - 1.96 * spec_age_all_grp2$robust_se)
spec_age_all_grp2$rr_high <- exp(spec_age_all_grp2$est + 1.96 * spec_age_all_grp2$robust_se)

spec_age_all_grp2$grp <- c("54-58", "59-63",
                           "64-68", "69-73",
                           "74-78", "79-81")

# Figure 2: Relative risk ratio of misclassification for age groups among patients without WHI-adjudicated outcomes (N = 106,559)

# one reference group
ggplot(spec_age_all_grp2, aes(x = grp, y = exp(est))) + 
  geom_point(color = "deepskyblue4") +
  geom_errorbar(aes(ymin = rr_low, ymax = rr_high), color = "deepskyblue2", width = .2) +
  geom_hline(yintercept = 1, color = "coral") +
  scale_y_continuous(name = "Relative risk ratio") + xlab("Age Group")

## consecutive groups
rr_age_grp <- t(sapply(1:6, function(i){
  spec_age_grp <- glm(cms_y ~ as.factor(age_grp), data = pop3[pop3$whi_y == 0 & (pop3$age_grp == i | pop3$age_grp == i+1),], family = "poisson")
  c(i, exp(coef(spec_age_grp)[2]), 
    coeftest(spec_age_grp,vcov=sandwich)[4],
    exp(coef(spec_age_grp)[2] - 1.96 * coeftest(spec_age_grp,vcov=sandwich)[4]),
    exp(coef(spec_age_grp)[2] + 1.96 * coeftest(spec_age_grp,vcov=sandwich)[4]),
    coeftest(spec_age_grp,vcov=sandwich)[8])
}))
colnames(rr_age_grp) = c("ref_grp", "rr", "robust_se", "rr_low", "rr_high", "p_val")

ggplot(as.data.frame(rr_age_grp), aes(x = as.factor(ref_grp), y = rr)) + 
  geom_point(color = "deepskyblue4") +
  geom_errorbar(aes(ymin = rr_low, ymax = rr_high), color = "deepskyblue2", width = .2) +
  geom_hline(yintercept = 1, color = "coral") +
  scale_x_discrete(breaks = c(1,2,3,4,5,6), 
                   labels = c("Age 54-58 \nvs 49-53", "Age 59-63 \nvs 54-58",
                              "Age 64-68 \nvs 59-63", "Age 69-73 \nvs 64-68",
                              "Age 74-78 \nvs 69-73", "Age 79-81 \nvs 74-78")) + 
  theme(axis.title.x = element_blank()) +
  scale_y_continuous(name = "Relative risk ratio")


# Table 6: Agreement between WHI Adjudicated vs. Medicare Defined Diagnoses of the Hormone  Trials among Patients who already had Medicare enrollment when started in WHI (N = 8,148)
pop6_hr_trt1 = pop6[pop6$hr_trt == 1,]
pop6_hr_trt0 = pop6[pop6$hr_trt == 0,]

eval_perf(pop6_hr_trt1)
eval_perf(pop6_hr_trt0)


# Table 7: Risk Decrease and Hazard Ratio of Colorectal Cancer using WHI vs. Medicare diagnosis among Patients who already had Medicare enrollment when started in WHI (N = 8,148)

# hr among all hormone patients for comparison
mod <- coxph(Surv(whi_time, whi_diag) ~ hr_trt, data = pop4)
exp(mod$coefficients) # HR with whi for all hormone patients

# true hazard ratio
mod <- coxph(Surv(whi_time, whi_diag) ~ hr_trt, data = pop6)
(whi_hr_true <- exp(mod$coefficients)) # HR with whi

mod <- coxph(Surv(cms_time, cms_diag)~hr_trt, data = pop6)
(cms_hr_true <- exp(mod$coefficients)) # HR with cms

cms_hr_true/whi_hr_true


# true risk difference
n_hr_trt0 = nrow(pop6_hr_trt0)
n_hr_trt1 = nrow(pop6_hr_trt1)

(cms_rd_true <- sum(pop6_hr_trt1$cms_diag == 1)/n_hr_trt1 - 
    sum(pop6_hr_trt0$cms_diag == 1)/n_hr_trt0)
(whi_rd_true <- sum(pop6_hr_trt1$whi_diag == 1)/n_hr_trt1 - 
    sum(pop6_hr_trt0$whi_diag == 1)/n_hr_trt0)
cms_rd_true - whi_rd_true

# bootstrap
set.seed(47)
boot_hr <- lapply(1:1000, function(i){
  # risk decrease
  boot1 = pop6_hr_trt1[sample(n_hr_trt1, replace = T),]
  boot0 = pop6_hr_trt0[sample(n_hr_trt0, replace = T),]
  whi_rd = sum(boot1$whi_diag == 1)/n_hr_trt1 -
    sum(boot0$whi_diag == 1)/n_hr_trt0
  cms_rd = sum(boot1$cms_diag == 1)/n_hr_trt1 -
    sum(boot0$cms_diag == 1)/n_hr_trt0
  boot = pop6[sample(n_hr, replace = T),]
  mod <- coxph(Surv(whi_time, whi_diag) ~ hr_trt, data = boot)
  whi_hr <- exp(mod$coefficients)[1]
  mod <- coxph(Surv(cms_time, cms_diag) ~ hr_trt, data = boot)
  cms_hr <- exp(mod$coefficients)[1] # HR with cms
  
  c(whi_rd, cms_rd,  whi_hr, cms_hr, 
    cms_rd - whi_rd, cms_hr/whi_hr)
})
table_7 <- do.call(rbind, boot_hr)
colnames(table_7) <- c("whi_rd", "cms_rd",  "whi_hr", "cms_hr", "rd_diff", "hr_ratio")
sapply(1:ncol(table_7), function(i) quantile(table_7[,i], c(0.025, 0.975)))
mean(table_7[,"rd_diff"] < cms_rd_true - whi_rd_true)
mean(table_7[,"hr_ratio"] < cms_hr_true/whi_hr_true)




