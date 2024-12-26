#######################################################################
##        Code for: "Limited adaptive responses in safety traits     ##
## support greater hydraulic risk under drier conditions"            ##
#######################################################################

library(ggplot2)
library(plyr)
library(ggspatial)
library(metafor)
library (car)
library(extrafont)
library(boot)
library(progress)
library(stringr)
library(lme4)
library(MuMIn)
library(emmeans)
library(smatr)
library(dplyr)
library(broom)
library(tidyr)
library(scales)
library(ggridges)

######################################################

library(dplyr)
library(ggplot2)
library(readr)

# Define geometric mean function ------------------------------------------

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Read master version of the data

hydplas.data.master <- read_tsv('HydPlasticity_final_2024-12-25.tsv')

# Check units

hydplas.data.master %>% 
  group_by(MeasName) %>% 
  distinct(MeasUnit) %>% 
  arrange(MeasName) %>% View()

# Set N to integer

hydplas.data.master$MeasR_N <- as.integer(hydplas.data.master$MeasR_N)
hydplas.data.master$MeasT_N <- as.integer(hydplas.data.master$MeasT_N)

# Variables

all_vars <- unique(hydplas.data.master$MeasName)[!is.na(unique(hydplas.data.master$MeasName))]

wp_vars<- c('Ypd','Ymd','Ptlp','P0','P50','Pc','Pe')
k_vars <- grep('k',all_vars,value=TRUE,ignore.case = TRUE)
size_vars <- grep('^S',all_vars,value=TRUE)
plc_rwc_vars <- c('RWC','PLC')
# Add RWD
plc_rwc_rwd_vars <- c(plc_rwc_vars,'RWD')
hv_struc_vars <- c('HV','R_L')
hsm_vars <- 'HSM'

# Variables that will be imputed
imput_vars <- c('MeasR_SE','MeasR_N','MeasT_SE','MeasT_N')
se_vars<- c('MeasR_SE','MeasT_SE')
n_vars<- setdiff(imput_vars,se_vars)
seperc_vars <- c('SE_perc_ref','SE_perc_treat')

# Calculate mean

hydplas_imput_simple <- hydplas.data.master %>% 
  # group by variable MeasName
  group_by(MeasName) %>% 
  # calculate means for variables of interest
  mutate(across(all_of(c(imput_vars,seperc_vars)),
                mean,na.rm=TRUE,.names="{.col}_mn")) %>%
  mutate(across(all_of(c(seperc_vars)),
                gm_mean,na.rm=TRUE,.names="{.col}_gm")) %>%
  # for variables with heterogeneous units, use relative SE in % to estimate SE  
  mutate(
    MeasR_SE_mn=ifelse(!MeasName%in%c(wp_vars,plc_rwc_rwd_vars,hsm_vars),
                       SE_perc_ref_mn/100*abs(MeasR_mean),MeasR_SE_mn),
    MeasT_SE_mn=ifelse(!MeasName%in%c(wp_vars,plc_rwc_rwd_vars,hsm_vars),
                       SE_perc_treat_mn/100*abs(MeasT_mean),MeasT_SE_mn)
  ) %>% 
  # round N imputations
  mutate(across(contains('_N'),round)) %>%
  # replace and indicate whether data are imputed (TRUE)
  mutate(
    MeasR_SE_imp=ifelse(is.na(MeasR_SE),MeasR_SE_mn,MeasR_SE),
    MeasR_N_imp=ifelse(is.na(MeasR_N),MeasR_N_mn,MeasR_N),
    MeasT_SE_imp=ifelse(is.na(MeasT_SE),MeasT_SE_mn,MeasT_SE),
    MeasT_N_imp=ifelse(is.na(MeasT_N),MeasT_N_mn,MeasT_N),
    MeasR_SE_ind=ifelse(is.na(MeasR_SE),TRUE,FALSE),
    MeasR_N_ind=ifelse(is.na(MeasR_N),TRUE,FALSE),
    MeasT_SE_ind=ifelse(is.na(MeasT_SE),TRUE,FALSE),
    MeasT_N_ind=ifelse(is.na(MeasT_N),TRUE,FALSE)
  )


# Write  

write_tsv(hydplas_imput_simple ,
          paste0('HydPlasticity_',Sys.Date(),'_simple_impute.tsv'))
saveRDS(hydplas_imput_simple,
        paste0('HydPlasticity_',Sys.Date(),'_simple_impute.rds'))
