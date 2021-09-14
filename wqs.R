##########################################################################################
#                                                                                        #
# R code for generating WQS index in Kang et al., 2021, Science of the Total Environment #
#                  Contact: Habyeong Kang (habyeongkang@gmail.com)                      #
#                                                                                        #
##########################################################################################


# load the file
setwd("C:\\file path")
library(haven)
nh <- read_sas("nhanes.sas7bdat")

dim(nh)
names(nh)

# select relevant variables
library(dplyr)
nh1 <- nh %>% select(SEQN, cycle, RIAGENDR, RIDAGEYR, RIDRETH1, smk, incm_g, bmi, act, lnucr, lnacr, egfr,
                    MEP3, MBP3, MIB3, MZP3, MC13, ECP3, MHH3, MOH3, COP3, CNP3, BP33, BPH3, MPB3, PPB3,
                    MEP1, MBP1, MIB1, MZP1, MC11, ECP1, MHH1, MOH1, COP1, CNP1, BP31, BPH1, MPB1, PPB1
                    )
dim(nh1)
names(nh1)

# make dummy variables for the categorical variables 
nh1$cycle2=ifelse(nh1$cycle==708,1,0)
nh1$cycle3=ifelse(nh1$cycle==910,1,0)
nh1$cycle4=ifelse(nh1$cycle==1112,1,0)
nh1$cycle5=ifelse(nh1$cycle==1314,1,0)
nh1$cycle6=ifelse(nh1$cycle==1516,1,0)
nh1$race2=ifelse(nh1$RIDRETH1==2,1,0)
nh1$race3=ifelse(nh1$RIDRETH1==3,1,0)
nh1$race4=ifelse(nh1$RIDRETH1==4,1,0)
nh1$race5=ifelse(nh1$RIDRETH1==5,1,0)
nh1$smk2=ifelse(nh1$smk==2,1,0)
nh1$smk3=ifelse(nh1$smk==3,1,0)
nh1$incm_g2=ifelse(nh1$incm_g==2,1,0)
nh1$incm_g3=ifelse(nh1$incm_g==3,1,0)
nh1$act2=ifelse(nh1$act==2,1,0)
nh1$act3=ifelse(nh1$act==3,1,0)
attach(nh1)

# input name of the chemical variables into "mix_name"
mix_name<-names(nh1)[13:26]

# WQS index for positive association
library(gWQS)
gwqs1 <- gwqs(egfr ~ wqs + factor(cycle)+RIAGENDR+RIDAGEYR+factor(RIDRETH1)+factor(smk)+factor(incm_g)+bmi+factor(act)+lnucr,
              mix_name = mix_name, data = nh, 
              q = 10, validation = 0.6, b = 100, b1_pos = T, 
              b1_constr = T, family = "gaussian", seed = 2021)
summary(gwqs1)

# WQS index for negative association
gwqs2 <- gwqs(egfr ~ wqs + factor(cycle)+RIAGENDR+RIDAGEYR+factor(RIDRETH1)+factor(smk)+factor(incm_g)+bmi+factor(act)+lnucr,
              mix_name = mix_name, data = nh, 
              q = 10, validation = 0.6, b = 100, b1_pos = F, 
              b1_constr = T, family = "gaussian", seed = 2021)
summary(gwqs2)

# bar plot
gwqs_barplot(gwqs2)
# scatter plot y vs wqs
gwqs_scatterplot(gwqs2)
# scatter plot residuals vs fitted values
gwqs_fitted_vs_resid(gwqs2)

