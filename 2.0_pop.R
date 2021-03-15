#=======================#
# Population selectopn  #
# 1. all 160K           #
# 2. FFS vs non-FFS     #
# 3. Hormone trial      #
# 4. >= 65 years old    #
#=======================#

source("1.0_utils.R")

# all WHI participants
pop0 <-fread("I:/2020-09-30/Demographics/dem_ctos_inv.dat")[, c("ID", "CTFLAG", "HRTFLAG", "DMFLAG", "OSFLAG", "AGE", "HRTARM")]
nrow(pop0) #161,808


# all baseline tables
diab <-fread("I:/2020-09-30/Demographics/f2_ctos_inv.dat")[, c("ID", "DIABTRT")]
hype <-fread("I:/2020-09-30/Med History/f30_ctos_inv.dat")[, c("ID", "HYPT", "HICHOLRP")]
meno <- fread("I:/2020-09-30/Med History/f31_ctos_inv.dat")[, c("ID", "MENO")]
alc_smok <- fread("I:/2020-09-30/Psychosocial/f34_ctos_inv.dat")[, c("ID", "ALC12DR", "SMOKEVR")]
horm <- fread("I:/2020-09-30/Medications/f43_ctos_inv.dat")[, c("ID", "TOTH")]
measure <- fread("I:/2020-09-30/Measurements/f80_ctos_inv.dat")
bmi <- measure %>% group_by(ID) %>% filter(F80DAYS == min(F80DAYS), row_number()==1) %>% select(c("ID", "BMI"))


# population with at least one baseline
pop1 <- list(pop0, meno, diab, hype, alc_smok, horm, bmi) %>% reduce(left_join, by = "ID")
pop1$meno_year = pop1$AGE - pop1$MENO  # age at enrollment - age of meno
sum(pop1$meno_year < 0, na.rm = TRUE) # 63 partcipants at the enrollment with no menopause yet
pop1$meno_year = ifelse(pop1$AGE - pop1$MENO < 0, NA, pop1$AGE - pop1$MENO)

pop1$ht_e = ifelse(pop1$HRTARM == 1 | pop1$HRTARM == 2, 1, 0)
pop1$ht_ep = ifelse(pop1$HRTARM == 3 |pop1$HRTARM == 4, 1, 0)

pop1$large_bmi = ifelse(pop1$BMI >= 25, 1, 0)


# subset to those in Roberta's algorithm for colorectal analysis
agr<-read_sas('U:/crctable_vde.sas7bdat')
colnames(agr)[1] <- "ID"
agr_incl <- agr[agr$crc_na !='X',]
colnames(agr_incl)
nrow(agr_incl) #107,743

pop2 <- agr_incl %>% select(-crc_na) %>% left_join(pop0, by = "ID")
nrow(pop2) #107,743

sum(pop2$HRTFLAG) #19204 hormone
sum(pop2$DMFLAG) #33,519 dietary
sum(pop2$HRTFLAG & pop2$DMFLAG) #5,630 in both
sum(pop2$OSFLAG) #60650 os
sum(pop2$HRTFLAG) + sum(pop2$DMFLAG) - sum(pop2$HRTFLAG & pop2$DMFLAG) + sum(pop2$OSFLAG) 

pop1$ffs <- ifelse(pop1$ID %in% pop2$ID, 1, 0)

# population dichotomized by baseline
pop3 <- list(agr_incl, pop1) %>% reduce(left_join, by = "ID")
median(pop3$meno_year, na.rm = TRUE)
pop3$long_meno <- ifelse(pop3$meno_year <= median(pop3$meno_year, na.rm = TRUE), 0, 1)

median(pop3$AGE, na.rm = TRUE)
median(pop3$AGE, na.rm = TRUE)
pop3$old_age <- ifelse(pop3$AGE <= median(pop3$AGE, na.rm = TRUE), 0, 1)
table(pop3$old_age, useNA = "always")

pop3$cms_y = ifelse(pop3$crc_yy == "X" | pop3$crc_yn == "X", 1, 0)
pop3$whi_y = ifelse(pop3$crc_yy == "X" | pop3$crc_ny == "X", 1, 0)

# age group
pop3$age_grp <- ifelse(pop3$AGE < min(pop3$AGE) + 5, 1, 
                       ifelse(pop3$AGE < min(pop3$AGE) + 5 * 2, 2,
                              ifelse(pop3$AGE < min(pop3$AGE) + 5 * 3, 3, 
                                     ifelse(pop3$AGE < min(pop3$AGE) + 5 * 4, 4,
                                            ifelse(pop3$AGE < min(pop3$AGE) + 5 * 5, 5, 
                                                   ifelse(pop3$AGE < min(pop3$AGE) + 5 * 6, 6,7))))))
table(pop3[pop3$whi_y == 0,]$age_grp)
table(pop3[pop3$whi_y == 0 & pop3$age_grp == 1,]$AGE)
table(pop3[pop3$whi_y == 0 & pop3$age_grp == 2,]$AGE)
table(pop3[pop3$whi_y == 0 & pop3$age_grp == 3,]$AGE)
table(pop3[pop3$whi_y == 0 & pop3$age_grp == 4,]$AGE)
table(pop3[pop3$whi_y == 0 & pop3$age_grp == 5,]$AGE)
table(pop3[pop3$whi_y == 0 & pop3$age_grp == 6,]$AGE)
table(pop3[pop3$whi_y == 0 & pop3$age_grp == 7,]$AGE)
nrow(pop3[pop3$whi_y == 0,])
table(pop3[pop3$whi_y == 1,]$large_bmi)
table(pop3[pop3$whi_y == 1,]$cms_n)

pop3$cms_n = 1 - pop3$cms_y

# relative risk ratio regression
summary(pop3[pop3$large_bmi == 1,]$BMI)


# hr among all hormone patients
pop4 <- pop1[pop1$HRTFLAG == 1,]
outcome <- fread("I:/2020-09-30/Outcomes/outc_ct_os_inv.dat")[, c("ID", "COLORECTAL", "COLORECTALDY")]

enr_date <- read.table("U:/ppt_rand_enroll_dates.dat", header = TRUE, fill = TRUE)
enr_date$enrldate_str <- as.character(enr_date$ENRLDATE)
enr_date$enrldate_fmt <- sapply(enr_date$enrldate_str, function(i) ifelse(nchar(i) == 7, paste0("0",i), i))
enr_date$whi_st <- as.Date(enr_date$enrldate_fmt, format = "%m%d%Y")

pop4 <- pop1[pop1$HRTFLAG == 1,] %>% left_join(outcome, by = "ID") %>% left_join(enr_date, by = "ID")
colnames(pop4)
pop4$trial_end <- as.Date(ifelse(pop4$HRTARM == 1 | pop4$HRTARM == 2, "2004-03-31", "2002-07-31"))
pop4$whi_diag_dt <- pop4$whi_st + pop4$COLORECTALDY

pop4$whi_diag <- ifelse(pop4$COLORECTAL == 1 & pop4$whi_diag_dt <= pop4$trial_end, 1, 0)
pop4$whi_time <- ifelse(pop4$COLORECTAL == 1 & pop4$whi_diag_dt <= pop4$trial_end, 
                        pop4$trial_end - pop4$whi_diag_dt, 
                        pop4$trial_end - pop4$whi_st)

pop4$hr_trt <- ifelse(pop4$HRTARM == 1 | pop4$HRTARM == 3, 1, 0)

## FFS vs non-FFS for hormone trials
# pop5 <- pop4[pop4$AGE >= 65,]
pop5 <- pop4[pop4$AGE >= 65,]
table(pop5$ffs)
table_5 <- summ_var(pop5)
table_5

#  secondary
date <- read.csv("U:/capstone_perturbed_crcdate.csv")
colnames(date)[1] <- "ID"
date$whi_st <- as.Date(date$WHISTARTDT, format = "%m/%d/%Y")
date$cms_st <- as.Date(date$CMSSTARTDT, format = "%m/%d/%Y")
date$cms_end <- as.Date(date$CMSSTOPDT, format = "%m/%d/%Y")
date$cms_diag_dt <- as.Date(date$P_CMSCRCDT, format = "%m/%d/%Y")

sum(!is.na(date$cms_diag_dt) & date$cms_diag_dt >= date$cms_end) # all cms diagnosis date before cms end
date$cms_diag <- (date$CMSCRC == "Yes")
same_date = date[date$cms_st == date$whi_st,]

## HT patients in roberta's algorithm with same cms and whi start date
pop6 <- pop4  %>% filter(ffs == 1) %>% inner_join(same_date, by = "ID")