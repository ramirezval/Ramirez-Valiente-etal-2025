#setwd("F:/Dropbox/Plasticity in hydraulics review/data/5_final_dataset") #Original directory
#setwd("F:/Dropbox/Paper Database plasticity") #My laptop
setwd("C:/Users/j.ramirez/Dropbox/Paper Database plasticity") #CREAF PC and laptop
setwd("C:/Users/jarv/Dropbox/Paper Database plasticity") #My miniPC
setwd("C:/Users/José Alberto/Dropbox/Paper Database plasticity") #PC Portátil INIA
setwd("C:/Users/jose.ramirez/Dropbox/Paper Database plasticity") #PC INIA
setwd("C:/Users/j.ramirez/OneDrive - CREAF/Plasticity database paper/Datasets") #OneDrive CREAF PC and laptop
setwd("C:/Users/jarv//OneDrive - CREAF/Plasticity database paper/Datasets") #OneDrive My miniPC
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


################################
### Loading the datasets    #####
################################

plas <- read_tsv('HydPlasticity_alltraits_imp_final_2024-12-25.tsv')
#plas2<-read.table("Alldata.csv", header=T, sep=",", fill = TRUE) #all data including imputed missing information
plas$SpeciesFix2 <- gsub(" ", "_", plas$SpeciesFix)

#########################
######## MAP ############
#########################

# Set the base theme
theme_set(theme_bw())

# Create a world map layer
mapWorld <- borders("world", colour = "black", fill = "white")

# Create the base plot
mp <- ggplot() + mapWorld +
  geom_point(aes(x = plas$LocLong, y = plas$LocLat, colour = plas$VarFactorAgg),
             shape = 21, size = 4, stroke = 1, fill = NA, show.legend = TRUE) +
  geom_point(colour = "white", size = 1.5)+
  scale_color_manual(
    values = c("#9C9C9C", "#cccc00", "#00c800", "#FF0000", "#008DF2", "black"),
    labels = c(expression(CO[2]),  # Subscript for CO₂
               "Light", 
               "Nutrient", 
               "Temperature", 
               "Water availability")  # Clean labels
  )+
  guides(colour = guide_legend(
    override.aes = list(fill = "white",  # Ensure legend circles have a white interior
                        colour = c("#9C9C9C", "#cccc00", "#00c800", "#FF0000", "#008DF2"),  # Outline colors for levels
                        shape = 21, stroke = 1, size = 5)  # Increase the size of the circles
  ))

# Restrict the map to specific latitudes and longitudes
mp <- mp + 
  coord_cartesian(ylim = c(-54, 79), xlim = c(-164, 179)) +  # Set latitude and longitude limits
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_rect(fill = "aliceblue"),  # Set map background to aliceblue
    legend.title = element_blank(),  # Remove the legend title
    legend.position = c(0.15, 0.18),  # Place the legend inside the map (adjust coordinates)
    legend.background = element_rect(fill = "white", colour = "black"),  # Add a white background with a border
    legend.text = element_text(size = 18, family = "Arial"),  # Set font to Arial and increase font size
    legend.key = element_rect(fill = "white", colour = NA),  # Ensure the legend key itself is white
    legend.key.height = unit(0.8, "cm"),  # Increase the spacing between legend levels
    axis.title.x = element_text(size = 29, margin = margin(t = 20)),  # Add more margin above x-axis title
    axis.title.y = element_text(size = 29, margin = margin(r = 18)),  # Add more margin to the right of y-axis title
    axis.text.x = element_text(size = 26, colour = "black", margin = margin(t = 10)),  # Add space between x-axis text and ticks
    axis.text.y = element_text(size = 26, colour = "black", margin = margin(r = 10))   # Add space between y-axis text and ticks
  ) +
  labs(x = "Longitude (º)", y = "Latitude (º)")  # Add axis titles with units

print(mp)
ggsave(filename = "Fig.S1.png", plot = mp, dpi = 300, width = 14, height = 7)




#############################################################
## Number of observations by trait and environmental factor #
#############################################################
table<-table(plas$MeasName,plas$VarFactorAgg)
data <- data.frame(table) 
data2<-subset(data, !(data$Var1 %in% c('Sa','Sl','Sv','Sw','R_S','R_L'))) # removing spatial studies
data2

height_P50 <- sum(data2$Freq[data2$Var1 == "P50"])
height_Ymd <- sum(data2$Freq[data2$Var1 == "Ymd"])
height_Ypd <- sum(data2$Freq[data2$Var1 == "Ypd"])
height_Ptlp <- sum(data2$Freq[data2$Var1 == "Ptlp"])
height_Ks <- sum(data2$Freq[data2$Var1 == "Ks"])
height_KL <- sum(data2$Freq[data2$Var1 == "KL"])
height_HV <- sum(data2$Freq[data2$Var1 == "HV"])

# Create the bar plot
ggp <- ggplot(data2, aes(x=reorder(Var1,desc(Freq)), y = Freq, fill = Var2, label = "n"))+ ylim(0,500)+ # Create stacked bar chart  # Create stacked bar chart
  # Set the y-axis limit
  ylim(0, 500) +
  
  # Add the bars to the plot
  geom_bar(stat = "identity") +
  
  # Add axis labels and title (currently no title or x-axis label)
  labs(title = "", x = "", y = "Number of observations") +
  
  # Customize the colors of the bars and legend labels
  scale_fill_manual(
    values = c("T_VPD" = "#FF0000",   # Red for Temperature/Vapor Pressure Deficit
               "C" = "#9C9C9C",       # Gray for CO2
               "N" = "#00c800",       # Green for Nutrients
               "L" = "#cccc00",       # Yellow for Light
               "W" = "#008DF2"),      # Blue for Water Availability
    labels = c(expression(CO[2]),    # CO₂ with subscript in the legend
               "Light", 
               "Nutrient", 
               "Temperature", 
               "Water availability")  # Custom legend labels
  ) +
  
  # Remove major and minor gridlines
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  # Customize the labels for the x-axis
  scale_x_discrete(labels = c('Ymd' = expression(Psi[MD]),  # Ψ_MD (Midday water potential)
                              'Ypd' = expression(Psi[PD]),  # Ψ_PD (Predawn water potential)
                              'Ptlp' = expression(Pi[TLP]), # π_TLP (Turgor Loss Point)
                              'P50' = expression(P[50]),    # P50 (Hydraulic threshold)
                              'Ks' = expression(K[s]),      # Ks (Specific conductivity)
                              'KL' = expression(K[L]),      # KL (Leaf-specific conductivity)
                              'HV' = expression(Hv),        # Hv (Vulnerability curve)
                              'kL' = expression(k[L]),      # kL
                              'P0' = expression(P[0]),      # P0
                              'Kt' = expression(K[t]),      # Kt
                              'ks' = expression(k[s]),      # ks
                              'Pe' = expression(P[e]),      # Pe
                              'Pc' = expression(P[c]),      # Pc
                              'HSM' = expression(SM))) +    # Safety Margin
  
  # Customize the panel background
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  
  # Adjust the plot margins
  theme(plot.margin = unit(c(0, 1.5, 0.5, 0.5), "cm")) +
  
  # Customize the legend
  theme(legend.title = element_blank(), 
        legend.position = c(0.85, 0.85),  # Position the legend inside the plot
        legend.background = element_blank(),  # Remove the background of the legend
        legend.key = element_rect(fill = "white", colour = NA),  # White background for legend keys
        legend.key.height = unit(0.8, "cm"),  # Increase the spacing between legend items
        legend.text = element_text(size = 16)) +  # Increase the font size of the legend text
  
  # Remove facet strip backgrounds (if facets are used)
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  
  # Customize axis labels and ticks
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(size = 17, color = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)),  # Customize x-axis text
        axis.text.y = element_text(size = 20, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 20)),  # Customize y-axis text
        axis.title.y = element_text(size = 24, color = "black", face = "plain", vjust = 2, margin = margin(t = 0, r = -5, b = 2, l = 0)),  # Adjust y-axis title
        axis.ticks.x = element_line(color = "black", size = 1),  # Customize x-axis ticks
        axis.ticks.y = element_line(color = "black", size = 1),  # Customize y-axis ticks
        axis.ticks.length = unit(0.25, "cm")) +  # Length of ticks
  
  # Customize plot border
  theme(plot.background = element_rect(fill = "white"), 
        panel.border = element_rect(colour = "black", fill = NA, size = 2))+
  
  # Add asterisks (*) above specific bars (e.g., significant results)
  geom_text(aes(x = "P50", y = height_P50 + 10, label = "*"), size = 13, color = "black", family = "Tahoma") +
  geom_text(aes(x = "Ymd", y = height_Ymd + 10, label = "*"), size = 13, color = "black", family = "Tahoma") +
  geom_text(aes(x = "Ypd", y = height_Ypd + 10, label = "*"), size = 13, color = "black", family = "Tahoma") +
  geom_text(aes(x = "Ptlp", y = height_Ptlp + 10, label = "*"), size = 13, color = "black", family = "Tahoma") +
  geom_text(aes(x = "Ks", y = height_Ks + 10, label = "*"), size = 13, color = "black", family = "Tahoma") +
  geom_text(aes(x = "KL", y = height_KL + 10, label = "*"), size = 13, color = "black", family = "Tahoma") +
  geom_text(aes(x = "HV", y = height_HV + 10, label = "*"), size = 13, color = "black", family = "Tahoma")

print(ggp)
ggsave(filename = "Fig.S2.png", plot = ggp, dpi = 300, width = 11, height = 8)

######################################################
### Important  first steps for the meta-analysis #####
######################################################

#file with measurements of selected traits in the same units
plastraits<- subset(plas, subset=(MeasName =="P50" |MeasName =="Ptlp" |MeasName =="HV"|
                                     MeasName =="KL"| MeasName =="Ks" | MeasName =="Ymd" | MeasName =="Ypd"))  #seleccionado solo el factor water
plastraits$SpeciesFix2 <- gsub(" ", "_", plastraits$SpeciesFix)
#creating a new dataset with only studies that use k max and not k native
#FlagK_def is a variable with the following categories: ktemporal_likely_native, max, NA, native, plc_p50_present_likely_max 
#We remove from the file the native categories
plastraits$FlagK_defNA <- addNA(plastraits$FlagK_def) # We use addNA function to consider NA as a factor in a new variable called FlagK_defNA
levels(plastraits$FlagK_defNA) #to check NA is now considered a factor level of FlagK_defNA
plastraits <- subset(plastraits, plastraits$FlagK_defNA != 'ktemporal_likely_native'& plastraits$FlagK_defNA !='native') #removing data of native ks
levels(plastraits$FlagK_defNA) #to check levels after removing k native. Subset function deletes the values, but does not remove the rows
droplevels(plastraits$FlagK_defNA) #To remove the empty rows we use droplevels function
levels(plastraits$FlagK_defNA) #checking the three factor levels in FlagK_defNA
plas=plastraits

###############################
### Calculating log ratio #####
###############################

plas$MeasR_SD_cor_imp<-plas$MeasR_SE_cor_imp*sqrt(plas$MeasR_N_cor_imp)
plas$MeasT_SD_cor_imp<-plas$MeasT_SE_cor_imp*sqrt(plas$MeasT_N_cor_imp)
plas <- subset(plas,plas$MeasR_SD_cor_imp>0) #deleting the negative value for MeasR_SD_imp
plas <- subset(plas,plas$MeasT_SD_cor_imp>0)#deleting the negative value for MeasT_SD_imp
plas$MeasR_var_cor_imp<-plas$MeasR_SD_cor_imp*plas$MeasR_SD_cor_imp
plas$MeasT_var_cor_imp<-plas$MeasT_SD_cor_imp*plas$MeasT_SD_cor_imp

#Estimating MD
plas<-escalc(measure="MD", n2i=plas$MeasR_N_cor_imp, n1i=plas$MeasT_N_cor_imp, m2i=plas$MeasR_mean_cor, m1i=plas$MeasT_mean_cor, sd2i=plas$MeasR_SD_cor_imp, sd1i=plas$MeasT_SD_cor_imp, 
             var.names=c("yi","vi"), add.measure=FALSE,
             append=TRUE, data=plas)
#Estimating log ratio
plas<-escalc(measure="ROM", n2i=plas$MeasR_N_cor_imp, n1i=plas$MeasT_N_cor_imp, m2i=plas$MeasR_mean_cor, m1i=plas$MeasT_mean_cor, sd2i=plas$MeasR_SD_cor_imp, sd1i=plas$MeasT_SD_cor_imp, 
             var.names=c("yi","vi"), add.measure=FALSE,
             append=TRUE, data=plas)

#######################################
# Calculating variances of log ratios #
#######################################

#Variance by MeasName for all studies
varall<-aggregate(yi ~ MeasName, plas, function(x) c(Var = var(x), Count = length(x)))
print(varall)

#Variance by MeasName considering only plasticity studies
plas.plas <- subset(plas, (plas$VarType %in% c('E')))
varplas<-aggregate(abs(yi) ~ MeasName, plas, function(x) c(Var = var(x), Count = length(x)))
print(varplas)

########################################
# Testing for differences in variances #
########################################
plas.traits <- subset(plas.plas, subset=(MeasName =="Ypd" |MeasName =="Ymd" |MeasName =="HV" |MeasName =="Ks" |MeasName =="KL"|MeasName =="Ptlp" |MeasName =="P50"))  #Selection 7 traits
fligner.test(yi ~ MeasName, data = plas.traits) #fligner test


################################################################################
#Density plots for log ratios for individual traits y differences in variance ##
################################################################################

plas2<- subset(plas, subset=(VarTypePost =="E"))  #seleccionando solo plasticity studies
plastraits<- subset(plas2, subset=(MeasName =="P50" |MeasName =="Ptlp" |MeasName =="HV"|
                                     MeasName =="KL"| MeasName =="Ks" | MeasName =="Ymd" | MeasName =="Ypd"))  #seleccionado solo el factor water
dens<-ggplot(plastraits, aes(x=yi, y=MeasName)) +
  geom_density_ridges(fill="lightblue",alpha=0.5)+scale_x_continuous(breaks = c(-1,0,1,2), limits= c(-1, 2))+
  geom_vline(xintercept=0, size=1, color="black",linetype="dashed")
dens
dens<-dens+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin = unit(c(1,1.5,0.5,0.5), "cm"))+
  labs(title=, subtitle=, caption=, tag=,y ="", x = "log-ratio")+
  theme(legend.title = element_blank(),legend.position = "none",strip.background = element_blank(),strip.text.x = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(size = 29, color = "black", margin = margin(t = 15, r = 80, b = 20, l = 55)), 
        axis.text.y = element_text(size = 36, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 23)),
        axis.title.x = element_text(size = 36, color = "black", face="plain", vjust=2, margin = margin(t = 10, r = 0, b = 2, l = 0)),
        axis.ticks.x=element_line(color = "black"),
        axis.ticks.y=element_line(color = "black"),
        axis.ticks.length = unit(0.25, "cm"))+  scale_y_discrete(labels = c('Ymd' = expression(Psi[MD]),
                                                                            'Ypd'   = expression(Psi[PD]),
                                                                            'Ptlp'= expression(Pi[TLP]),
                                                                            'P50'= expression(P[50]),
                                                                            'Ks'= expression(K[s]),
                                                                            'KL'= expression(K[L]),
                                                                            'HV'= expression(Hv)))
print(dens)
ggsave(filename = "Fig.1a.png", plot = dens, dpi = 300, width = 8, height = 11)

################################################################################
# Meta-analysis for absolute log ratios for individual traits across factors ###
################################################################################
#Meta-analysis with metafor and bootstrapping for confidence intervals

# Function to perform an analysis with bootstrap for a specific trait
run_analysis_with_bootstrap <- function(data, trait, tot_rep = 10000) {
  
  # Filter the dataset for the specific trait
  subset_data <- subset(data, MeasName == trait)
  
  # Fit a random-effects meta-analysis model using rma.mv
  model <- rma.mv(abs(yi), vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2), 
                  control = list(optimizer = "optim", optmethod = "Nelder-Mead"), 
                  data = subset_data)
  
  # Define the bootstrap function
  boot_func <- function(plasboot, indices, prog) {
    
    # Update progress bar with each iteration
    prog$tick()
    
    sel <- plasboot[indices, ]
    res <- try(suppressWarnings(rma.mv(abs(yi), vi, random = list(~ 1 | SpeciesFix, ~ 1 | PaperID), data = sel)), silent = TRUE)
    
    # Handle cases where the model fails to fit
    if (inherits(res, "try-error")) {
      NA
    } else {
      c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
    }
  }
  
  # Initialize a progress bar
  pb <- progress_bar$new(total = tot_rep + 1)
  
  # Run the bootstrap
  boot_results <- boot(subset_data, boot_func, tot_rep, prog = pb)
  
  # Calculate confidence intervals
  ci_results <- boot.ci(boot_results, index = 1:2)
  
  # Extract specific parameters
  estimate <- coef(model)[1]
  ci_bca <- ci_results$bca[4:5]  # Extract BCa confidence interval limits
  sample_size <- model$k  # Extract the sample size
  
  # Return the model, bootstrap results, confidence intervals, and extracted values
  list(model = model, 
       bootstrap = boot_results, 
       confidence_intervals = ci_results, 
       extracted_values = list(estimate = estimate, 
                               bca_interval = ci_bca,
                               sample_size = sample_size))
}

# Function to run the analysis for each level of MeasName and subsets of VarContextAgg
run_analysis_with_subsets <- function(data, tot_rep = 10000) {
  
  # Define the three subsets
  subsets <- list(
    "No Spatial" = subset(data, VarContextAgg != "spatial"),
    "Only Spatial" = subset(data, VarContextAgg == "spatial"),
    "All Studies" = data
  )
  
  # Initialize a list to store results
  all_results <- list()
  
  # Loop through each subset
  for (subset_name in names(subsets)) {
    cat("Analyzing subset:", subset_name, "\n")  # Print progress
    
    # Get the current subset of data
    subset_data <- subsets[[subset_name]]
    
    # Get unique levels of MeasName in the subset
    meas_levels <- unique(subset_data$MeasName)
    
    # Loop through each level of MeasName
    for (trait in meas_levels) {
      cat("  Analyzing level:", trait, "\n")
      
      # Run the analysis for the current level in the current subset
      res <- run_analysis_with_bootstrap(subset_data, trait, tot_rep)
      
      # Extract specific values from the result
      extracted_values <- res$extracted_values
      
      # Add trait name (MeasName) and subset name to the extracted values
      extracted_values$MeasName <- trait
      extracted_values$Subset <- subset_name
      
      # Create a data frame with the extracted values
      results_df <- data.frame(
        Subset = extracted_values$Subset,
        MeasName = extracted_values$MeasName,
        Estimate = extracted_values$estimate,
        Lower_CI = extracted_values$bca_interval[1],
        Upper_CI = extracted_values$bca_interval[2],
        Sample_Size = extracted_values$sample_size
      )
      
      # Add the result to the list of all results
      all_results[[paste(subset_name, trait, sep = "_")]] <- results_df
    }
  }
  
  # Combine all results into a single data frame
  final_results_df <- do.call(rbind, all_results)
  
  return(final_results_df)
}

# Call the function to create the dataset with all levels of MeasName and subsets of VarContextAgg
results_all_subsets_df <- run_analysis_with_subsets(plas, tot_rep = 10000)  # 'plas' is the dataset for analysis

# Display the final dataset
print(results_all_subsets_df)

abslogratio<-results_all_subsets_df

# Save the results_all_subsets_df dataset as a .tsv file
write.table(results_all_subsets_df, 
            file = "abslogratio.tsv",  # Specify the output file name
            sep = "\t",                           # Use tab as the separator
            row.names = FALSE,                    # Do not include row names
            quote = FALSE)




#Create Plot
abslogratio$MeasName <- factor(abslogratio$MeasName, levels = c("HV", "KL", "Ks", "P50", "Ptlp", "Ymd", "Ypd"))

###### 3 vertical points

# Add the "Significant" column to the dataframe
# It is defined as significant if the confidence interval does not include the value 0
abslogratio$Significant <- ifelse(abslogratio$Lower_CI > 0 | abslogratio$Upper_CI < 0, TRUE, FALSE)

# Create the plot with ggplot and inverted axes
p <- ggplot(abslogratio, aes(y = MeasName, x = Estimate, fill = factor(Subset))) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(shape = 21, size = 7, color = "black", stroke = 1, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("gray", "white", "black")) +
  labs(y = "", x = "|log-ratio|") +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text.y = element_text(size = 32, color = "black", margin = margin(t = 10)),  # Space between letters and ticks on the y-axis
        axis.text.x = element_text(size = 30, color = "black", margin = margin(r = 5)),  # Space between numbers and ticks on the x-axis
        axis.title.x = element_text(size = 38, color = "black", margin = margin(r = 10)),  # Space between the x-axis title and numbers
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.15, "cm")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(-0.4, 1, by = 0.2)) +
  # Add text with sample size, in bold if significant
  geom_text(aes(label = Sample_Size, 
                x = Upper_CI + 0.05, fontface = ifelse(Significant, "bold", "plain")), 
            position = position_dodge(width = 0.8), size = 7, color = "black", hjust = 0) +
  scale_y_discrete(labels = c("Hv", expression(K[L]), expression(K[S]), expression(P[50]), 
                              expression(Π[TLP]), expression(Ψ[MD]), expression(Ψ[PD])))
p

# Save the plot with ggsave at 300 dpi
ggsave(filename = "Fig.1b.png", plot = p, dpi = 300, width = 12, height = 8)

###### 2 puntos vertical




########################################################################
######### FIG. 2 #########
########################################################################



#Análisis of yi by Trait, VarContextAgg y VarFactorAgg

plasboot<-plas.C<- subset(plas, subset=(VarFactorAgg =="C"))  #Se puede automatizar para todo o si queremos se pueden hacer subsets para cada VarFactorAgg

run_analysis_with_bootstrap <- function(data, trait, factor_level, tot_rep = 1000, min_sample_size = 3) {
  
  # Filtrar el dataset para el rasgo específico y el nivel del factor
  subset_data <- subset(data, MeasName == trait & VarFactorAgg == factor_level)
  
  # Si no hay datos (sample size = 0), devolver NA para Estimate, Lower_CI y Upper_CI, y 0 para Sample_Size
  if (nrow(subset_data) == 0) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = 0
    ))
  }
  
  # Si el tamaño de muestra es insuficiente, devolver NA sin intentar ajustar el modelo
  if (nrow(subset_data) < min_sample_size) {
    cat("  Tamaño de muestra insuficiente para trait:", trait, "y Factor:", factor_level, "- Sample size:", nrow(subset_data), "\n")
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Ajustar el modelo con rma.mv
  model <- tryCatch(
    {
      rma.mv(abs(yi), vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2), 
             control = list(optimizer="optim", optmethod="Nelder-Mead"), 
             data = subset_data)
    },
    error = function(e) NULL
  )
  
  # Si el modelo falla, retornar un data.frame con valores NA
  if (is.null(model)) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Función de bootstrap
  boot_func <- function(plasboot, indices, prog) {
    prog$tick()
    sel <- plasboot[indices, ]
    res <- try(suppressWarnings(rma.mv(abs(yi), vi, random = list(~ 1 | SpeciesFix, ~ 1 | PaperID), data = sel)), silent = TRUE)
    if (inherits(res, "try-error")) NA else c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
  # Inicializar barra de progreso
  pb <- progress_bar$new(total = tot_rep + 1)
  
  # Ejecutar bootstrap
  boot_results <- tryCatch(
    boot(subset_data, boot_func, tot_rep, prog = pb),
    error = function(e) NULL
  )
  
  # Si el bootstrap falla, retornar un data.frame con valores NA
  if (is.null(boot_results)) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Calcular intervalos de confianza
  ci_results <- tryCatch(
    boot.ci(boot_results, type = "bca", index = 1:2),
    error = function(e) NULL
  )
  
  # Si el cálculo de intervalos de confianza falla, retornar un data.frame con valores NA
  if (is.null(ci_results)) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Extraer parámetros específicos
  estimate <- coef(model)[1]
  ci_bca <- ci_results$bca[4:5]
  sample_size <- model$k
  
  # Retornar resultados
  data.frame(
    Estimate = estimate,
    Lower_CI = ci_bca[1],
    Upper_CI = ci_bca[2],
    Sample_Size = sample_size
  )
}

run_analysis_with_subsets <- function(data, tot_rep = 1000) {
  
  # Definir los tres subconjuntos
  subsets <- list(
    "No Spatial" = subset(data, VarContextAgg != "spatial"),
    "Only Spatial" = subset(data, VarContextAgg == "spatial"),
    "All Studies" = data
  )
  
  all_results <- list()
  
  # Obtener todas las combinaciones únicas de MeasName y VarFactorAgg en el conjunto completo
  meas_levels <- unique(data$MeasName)
  factor_levels <- unique(data$VarFactorAgg)
  
  for (subset_name in names(subsets)) {
    cat("Analizando subconjunto:", subset_name, "\n")
    subset_data <- subsets[[subset_name]]
    
    # Iterar sobre todas las combinaciones posibles de MeasName y VarFactorAgg
    for (trait in meas_levels) {
      for (factor_level in factor_levels) {
        
        cat("  Analizando nivel:", trait, "y Factor:", factor_level, "\n")
        
        # Ejecutar el análisis solo si hay datos en el subconjunto para esa combinación
        res <- run_analysis_with_bootstrap(subset_data, trait, factor_level, tot_rep)
        
        # Crear el data.frame con los resultados y agregar Subset, MeasName y VarFactorAgg
        results_df <- data.frame(
          Subset = subset_name,
          MeasName = trait,
          VarFactorAgg = factor_level,
          Estimate = res$Estimate,
          Lower_CI = res$Lower_CI,
          Upper_CI = res$Upper_CI,
          Sample_Size = res$Sample_Size
        )
        
        all_results[[paste(subset_name, trait, factor_level, sep = "_")]] <- results_df
      }
    }
  }
  
  final_results_df <- do.call(rbind, all_results)
  return(final_results_df)
}

# Ejecutar la función
logratiofactors <- run_analysis_with_subsets(plasboot, tot_rep = 1000)

# Mostrar la base de datos final
print(logratiofactors)

# Save the results_all_subsets_df dataset as a .tsv file
write.table(logratiofactors, 
            file = "logratiofactors.tsv",  # Specify the output file name
            sep = "\t",                           # Use tab as the separator
            row.names = FALSE,                    # Do not include row names
            quote = FALSE)






###################################################

# Function to perform bootstrap analysis for a specific Trait and Factor level
run_analysis_with_bootstrap <- function(data, trait, factor_level, tot_rep = 1000, min_sample_size = 3) {
  
  # Filter the dataset for the specific Trait and Factor level
  subset_data <- subset(data, MeasName == trait & VarFactorAgg == factor_level)
  
  # If there are no data (sample size = 0), return NA for Estimate, Lower_CI, Upper_CI, and 0 for Sample_Size
  if (nrow(subset_data) == 0) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = 0
    ))
  }
  
  # If the sample size is insufficient, return NA without fitting the model
  if (nrow(subset_data) < min_sample_size) {
    cat("  Insufficient sample size for Trait:", trait, "and Factor:", factor_level, "- Sample size:", nrow(subset_data), "\n")
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Fit the random-effects meta-analysis model using rma.mv
  model <- tryCatch(
    {
      rma.mv(abs(yi), vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2), 
             control = list(optimizer = "optim", optmethod = "Nelder-Mead"), 
             data = subset_data)
    },
    error = function(e) NULL
  )
  
  # If model fitting fails, return a data frame with NA values
  if (is.null(model)) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Bootstrap function
  boot_func <- function(plasboot, indices, prog) {
    prog$tick()
    sel <- plasboot[indices, ]
    res <- try(suppressWarnings(rma.mv(abs(yi), vi, random = list(~ 1 | SpeciesFix, ~ 1 | PaperID), data = sel)), silent = TRUE)
    if (inherits(res, "try-error")) NA else c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
  # Initialize the progress bar
  pb <- progress_bar$new(total = tot_rep + 1)
  
  # Run the bootstrap
  boot_results <- tryCatch(
    boot(subset_data, boot_func, tot_rep, prog = pb),
    error = function(e) NULL
  )
  
  # If bootstrap fails, return a data frame with NA values
  if (is.null(boot_results)) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Calculate confidence intervals
  ci_results <- tryCatch(
    boot.ci(boot_results, type = "bca", index = 1:2),
    error = function(e) NULL
  )
  
  # If confidence interval calculation fails, return a data frame with NA values
  if (is.null(ci_results)) {
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Extract specific parameters
  estimate <- coef(model)[1]
  ci_bca <- ci_results$bca[4:5]
  sample_size <- model$k
  
  # Return the results
  data.frame(
    Estimate = estimate,
    Lower_CI = ci_bca[1],
    Upper_CI = ci_bca[2],
    Sample_Size = sample_size
  )
}

# Function to run analysis for each combination of Trait and Factor level
run_analysis_with_subsets <- function(data, tot_rep = 1000) {
  
  # Define three subsets
  subsets <- list(
    "No Spatial" = subset(data, VarContextAgg != "spatial"),
    "Only Spatial" = subset(data, VarContextAgg == "spatial"),
    "All Studies" = data
  )
  
  all_results <- list()
  
  # Get all unique combinations of MeasName and VarFactorAgg
  meas_levels <- unique(data$MeasName)
  factor_levels <- unique(data$VarFactorAgg)
  
  for (subset_name in names(subsets)) {
    cat("Analyzing subset:", subset_name, "\n")
    subset_data <- subsets[[subset_name]]
    
    # Iterate over all possible combinations of MeasName and VarFactorAgg
    for (trait in meas_levels) {
      for (factor_level in factor_levels) {
        
        cat("  Analyzing Trait:", trait, "and Factor:", factor_level, "\n")
        
        # Run the analysis only if there is data for that combination
        res <- run_analysis_with_bootstrap(subset_data, trait, factor_level, tot_rep)
        
        # Create a data frame with the results and add Subset, MeasName, and VarFactorAgg
        results_df <- data.frame(
          Subset = subset_name,
          MeasName = trait,
          VarFactorAgg = factor_level,
          Estimate = res$Estimate,
          Lower_CI = res$Lower_CI,
          Upper_CI = res$Upper_CI,
          Sample_Size = res$Sample_Size
        )
        
        all_results[[paste(subset_name, trait, factor_level, sep = "_")]] <- results_df
      }
    }
  }
  
  # Combine all results into a single data frame
  final_results_df <- do.call(rbind, all_results)
  return(final_results_df)
}

# Run the function
logratiofactors <- run_analysis_with_subsets(plas, tot_rep = 1000)

# Display the final dataset
print(logratiofactors)

# Save the results as a .tsv file
write.table(logratiofactors, 
            file = "logratiofactors.tsv",  # Specify the output file name
            sep = "\t",                   # Use tab as the separator
            row.names = FALSE,            # Do not include row names
            quote = FALSE)                # Do not add quotes


######################################
### Gráfico horizontal de figura 2 ###############
######################################



####################################
##### Figure 3a - queda por se completada
####################################

#Estimate global Hydraulic safety margin (HSM) as the difference in the bootstrap distribution of P50 and Ymd
write.csv(resbootmd.plas.w.P50[["t"]], file = "plas.w.P50.csv")
boot.plas.w.P50<-read.table("plas.w.P50.csv", header=T, sep=",", fill = TRUE)
write.csv(resbootmd.plas.w.Ymd[["t"]], file = "plas.w.Ymd.csv")
boot.plas.w.Ymd<-read.table("plas.w.Ymd.csv", header=T, sep=",", fill = TRUE)
diff=boot.plas.w.Ymd$V1-boot.plas.w.P50$V1
hist(diff)
median(diff, na.rm=TRUE)
quantile(diff,c(0.05), na.rm=TRUE)
View(diff)





#######################################
##### Análisis y Figura 3b
#######################################

library(dplyr)
library(stringr)

# Step 1: Create the 'Combined' variable in the original dataframe
plas$Combined <- str_c(plas$PaperID, "_", plas$SpeciesFix2, "_", plas$VarContextAgg)

# Step 2: Filter for VarFactorAgg equal to "W"
plasboot <- filter(plas, VarFactorAgg == "W")

# Step 3: Define a function to fit the model and extract estimates
get_estimates <- function(data, measure_col, var_col) {
  tryCatch({
    res <- rma(get(measure_col), get(var_col), mods = ~Combined - 1, 
               control = list(optimizer = "optim", optmethod = "Nelder-Mead"), data = data)
    coef(res)
  }, error = function(e) NA)  # Return NA in case of an error
}

# Step 4: Apply the function for different combinations and store the results
results <- list(
  P50_MeasT = get_estimates(filter(plasboot, MeasName == "P50"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  Ymd_MeasT = get_estimates(filter(plasboot, MeasName == "Ymd"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  P50_MeasR = get_estimates(filter(plasboot, MeasName == "P50"), "MeasR_mean_cor", "MeasR_var_cor_imp"),
  Ymd_MeasR = get_estimates(filter(plasboot, MeasName == "Ymd"), "MeasR_mean_cor", "MeasR_var_cor_imp")
)

# Step 5: Combine results for MeasT and MeasR into dataframes
combine_results <- function(P50, Ymd, environment) {
  df_P50 <- data.frame(Combined = names(P50), P50 = unname(P50))
  df_Ymd <- data.frame(Combined = names(Ymd), Ymd = unname(Ymd))
  merged <- full_join(df_P50, df_Ymd, by = "Combined") %>%
    mutate(Environment = environment)
  return(merged)
}

final_MeasT <- combine_results(results$P50_MeasT, results$Ymd_MeasT, "T")
final_MeasR <- combine_results(results$P50_MeasR, results$Ymd_MeasR, "R")

# Step 6: Add columns for 'VarContextAgg' and 'Species'
extract_context <- function(combined_name) {
  case_when(
    str_detect(combined_name, "spatial") ~ "spatial",
    str_detect(combined_name, "experimental") ~ "experimental",
    str_detect(combined_name, "trial") ~ "trial",
    str_detect(combined_name, "temporal") ~ "temporal",
    TRUE ~ NA_character_
  )
}

extract_species <- function(combined_name) {
  str_extract(combined_name, "[A-Za-z]+_[A-Za-z]+(?=_(spatial|experimental|trial|temporal))")
}

final_MeasT <- final_MeasT %>%
  mutate(VarContextAgg = sapply(Combined, extract_context),
         Species = sapply(Combined, extract_species))

final_MeasR <- final_MeasR %>%
  mutate(VarContextAgg = sapply(Combined, extract_context),
         Species = sapply(Combined, extract_species))

# Step 7: Assign unique IDs
assign_id <- function(df) {
  df %>%
    mutate(ID_base = str_extract(str_remove(Combined, "^Combined"), "^[^_]+_[^_]+_[^_]+"),
           ID = dense_rank(ID_base)) %>%
    select(-ID_base)
}

final_MeasT <- assign_id(final_MeasT)
final_MeasR <- assign_id(final_MeasR)

# Step 8: Combine both datasets into 'final'
final <- bind_rows(final_MeasT, final_MeasR)

# Step 9: Add 'HSM' as the difference between Ymd and P50
final <- final %>%
  mutate(HSM = if_else(!is.na(Ymd) & !is.na(P50), Ymd - P50, NA_real_))

# Step 10: Save the results
write.csv(final, "final_results.csv", row.names = FALSE)

# Display the final dataframe
print("Final combined dataframe with HSM:")
print(final)


# Modelo con lmer
test <- lmer(HSM ~ Environment + (1|ID) + (1|Sp) + (1|VarContextAgg), data = final)

# Análisis
Anova(test)
summary(test)


emmeans(test, pairwise~Environment)
ggplot(final, aes(x = Environment, y = HSM))+ geom_boxplot(aes(fill= Environment))









#HSM
plot=ggplot(final, aes(x = factor(Environment), y = HSM, color = factor(Environment))) + 
  geom_boxplot(width=2, lwd=1,) + 
  scale_color_manual(values=c("#00B4F7", "#FF0000"))+
  geom_point(position = "jitter",size=4,shape=20, stroke=1) +
  #geom_line(aes(group=paired),color="black")+
  stat_summary(fun=mean, geom="crossbar", size=1, color="black", fill="red", linetype = "dashed")+
  ylab(expression(paste('HSM ', (MPa))))+
  theme_bw()+theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(), panel.background = element_rect(fill = "white", colour = "black", size=2),
        axis.text.x = element_text(size = 54, color = "black", margin = margin(t = 5, r = 50, b = 20, l = 55)), 
        axis.text.y = element_text(size = 54, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 0)),
        axis.title.x = element_text(size = 58, color = "black", face="plain", vjust=2, margin = margin(t = 10, r = 0, b = 2, l = 0)),
        axis.title.y = element_text(size = 58, color = "black", face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x=element_line(color = "black", size=1),
        axis.ticks.y=element_line(color = "black",  size=1),
        axis.ticks.length = unit(0.25, "cm"), 
        legend.margin=margin(l = 0.6, unit='cm'))+theme(legend.position = "none")
plot



###############################################################################################
############### Fig. 4. Correlations of log ratios by paper and environmental factor ##################
###############################################################################################

library(metafor)
library(dplyr)
library(stringr)
library(lme4)
library(car)

# Crear la variable 'Combined' en el dataframe original
plas.w<- subset(plas, VarFactorAgg =="W")
plas.w$Combined <- str_c(plas.w$PaperID, "_", plas.w$SpeciesFix2, "_", plas.w$VarContextAgg,"_", plas.w$MeasName)

res<-rma(yi, vi, mods = ~Combined - 1, 
            control = list(optimizer="optim", optmethod="Nelder-Mead"),data=plas.w)
res
           
data <- data.frame(Moderator = names(coef(res)), Estimate = coef(res))               
# Reemplazar todos los espacios por guiones bajos en la columna 'Moderator'
data$Moderator <- gsub(" ", "_", data$Moderator)
# Eliminar el código "_et_al" de la columna Moderator
data$Moderator <- gsub("_et_al", "", data$Moderator)
data$Moderator <- gsub("CombinedCarins_Murphy", "CombinedCarins-Murphy", data$Moderator)
data$Moderator <- gsub("CombinedEpron_", "CombinedEpron_2015_", data$Moderator)
data$Moderator <- gsub("CombinedLamy_", "CombinedLamy_2014_", data$Moderator)

# Extraer el último valor (MeasName)
data$MeasName <- sapply(strsplit(data$Moderator, "_"), function(x) tail(x, 1))

# Extraer el antepenúltimo valor (VarContextAgg)
data$VarContextAgg <- sapply(strsplit(data$Moderator, "_"), function(x) x[length(x) - 1])

# Extraer los valores 4 y 5 desde el final de la columna Moderator y combinarlos para formar la columna Species
data$Species <- sapply(strsplit(data$Moderator, "_"), function(x) {
  if (length(x) >= 5) paste(x[length(x) - 3], x[length(x) - 2], sep = "_") else NA
})

# Extraer los valores 1 y 2 desde el inicio, eliminando la palabra "Combined"
data$PaperID <- sapply(strsplit(data$Moderator, "_"), function(x) {
  if (length(x) >= 3) paste(x[1:3], collapse = "_") else NA
})

# Eliminar "Combined" de la nueva columna
data$PaperID <- gsub("^Combined", "", data$PaperID)

# Verificar las primeras filas del dataframe con la nueva columna
head(data)

data <- data %>%
  pivot_wider(
    id_cols = c(VarFactorAgg, VarContextAgg, Species, PaperID), # Variables para agrupar
    names_from = MeasName,  # Los nombres de las columnas serán los niveles de MeasName
    values_from = Estimate,  # Los valores de las nuevas columnas provendrán de Estimate
    names_glue = "{MeasName}_yi"  # Agregar el sufijo "_yi" a los nombres de las columnas
  )

# Verificar el nuevo dataframe
head(data)


data.l<- subset(data, subset=(VarFactorAgg =="L"))
data.n<- subset(data, subset=(VarFactorAgg =="N"))
data.w<- subset(data, subset=(VarFactorAgg =="W"))
data.c<- subset(data, subset=(VarFactorAgg =="C"))
data.w.plas<- subset(data.w, VarContextAgg != "spatial")

fit1 <- lmer(P50_yi ~ Ypd_yi +(1|VarContextAgg)+(1|PaperID), data = data) #change variables to obtain the different regressions

Anova(fit1)
summary(fit1)

#Graphics
# Function to darken colors
darken_color <- function(color, amount = 0.3) {
  col <- col2rgb(color)
  col <- col * (1 - amount)
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}

# Filtering groups with less than 5 valid observations
data_filtered <- data %>%
  filter(!is.na(Ypd_yi) & !is.na(P50_yi)) %>%
  group_by(VarFactorAgg) %>%
  filter(n() >= 5) %>%
  ungroup()

# Step 1: Calculating correlations
correlations <- data_filtered %>%
  group_by(VarFactorAgg) %>%
  summarize(cor_test = list(tidy(cor.test(Ypd_yi, P50_yi)))) %>%
  unnest(cor_test) %>%
  ungroup()

# Step 2: Filtering significant groups (p-value < 0.05)
significant_groups <- correlations %>%
  filter(p.value < 0.05) %>%
  pull(VarFactorAgg)

# Step 3: Calculating global correlation and significance
global_cor_test <- tidy(cor.test(data_filtered$Ypd_yi, data_filtered$P50_yi))
global_p_value <- global_cor_test$p.value

# Definition of colors
colores <- c("#9C9C9C", "#CCCC00", "#00C800", "#FF0200", "#008DF2")
names(colores) <- c("C", "L", "N", "TVPD", "W")

colores_oscuros <- sapply(colores, darken_color, amount = 0.3)
names(colores_oscuros) <- names(colores)

# Spet 4: Creating the figure with `ggplot2`
plot <- ggplot(data, aes(x = Ypd_yi, y = P50_yi, fill = factor(VarFactorAgg))) +
  geom_point(shape = 21, size = 13, stroke = 1, color = "black") +
  scale_fill_manual(values = colores) +
  scale_color_manual(values = colores) +
  geom_abline(slope = 1, linetype = "dashed", size = 2) +
  xlab(expression(paste('log-ratio ',Psi[PD]))) +
  ylab(expression(paste('log-ratio ', P[50]))) +
  theme_minimal()
plot

# Adding significant regression lines
for (grupo in significant_groups) {
  plot <- plot + geom_smooth(data = filter(data_filtered, VarFactorAgg == grupo), method = "lm", size = 2, se = TRUE, color = colores_oscuros[grupo])
}
#  Adding significant global regression line. Change P following mixed model outcome
if (global_p_value < 0.05) {
  plot <- plot + geom_smooth(data = data_filtered, method = "lm", se = TRUE, linetype = "solid", size = 2, aes(group = 1), color = "black")
}

print(plot)

plot+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(), panel.background = element_rect(fill = "white", colour = "black", size=2),
        axis.text.x = element_text(size = 54, color = "black", margin = margin(t = 5, r = 50, b = 20, l = 55)), 
        axis.text.y = element_text(size = 54, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 0)),
        axis.title.x = element_text(size = 58, color = "black", face="plain", vjust=2, margin = margin(t = 10, r = 0, b = 2, l = 0)),
        axis.title.y = element_text(size = 58, color = "black", face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x=element_line(color = "black", size=1),
        axis.ticks.y=element_line(color = "black",  size=1),
        axis.ticks.length = unit(0.25, "cm"), 
        legend.margin=margin(l = 0.6, unit='cm')) + 
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-1,1.6), ylim = c(-1, 1.6)) + 
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5)) + scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5))



###############################################################
### Regressions between species means within environments #####
###############################################################

# Step 1: Create the 'Combined' variable in the original dataframe
plas$Combined <- str_c(plas$PaperID, "_", plas$SpeciesFix2)

# Step 3: Create the 'Combined' variable in the dataframe
plas$Combined <- str_c(plas$PaperID, "_", plas$SpeciesFix2)

# Step 4: Filter for VarFactorAgg equal to "W"
plasboot <- subset(plas, subset = (VarFactorAgg == "W"))

# Step 5: Define the function to fit the model and extract estimates
get_estimates <- function(data, measure_col, var_col) {
  tryCatch({
    res <- rma(
      yi = data[[measure_col]], 
      vi = data[[var_col]], 
      mods = ~Combined - 1, 
      control = list(optimizer = "optim", optmethod = "Nelder-Mead"), 
      data = data
    )
    coef(res)  # Return model coefficients
  }, error = function(e) NA)  # Return NA if an error occurs
}

# Step 6: Apply the function for different combinations and store the results
results <- list(
  Ypd_MeasT = get_estimates(filter(plasboot, MeasName == "Ypd"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  Ymd_MeasT = get_estimates(filter(plasboot, MeasName == "Ymd"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  P50_MeasT = get_estimates(filter(plasboot, MeasName == "P50"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  Ptlp_MeasT = get_estimates(filter(plasboot, MeasName == "Ptlp"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  Ks_MeasT = get_estimates(filter(plasboot, MeasName == "Ks"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  KL_MeasT = get_estimates(filter(plasboot, MeasName == "KL"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  HV_MeasT = get_estimates(filter(plasboot, MeasName == "HV"), "MeasT_mean_cor", "MeasT_var_cor_imp"),
  Ypd_MeasR = get_estimates(filter(plasboot, MeasName == "Ypd"), "MeasR_mean_cor", "MeasR_var_cor_imp"),
  Ymd_MeasR = get_estimates(filter(plasboot, MeasName == "Ymd"), "MeasR_mean_cor", "MeasR_var_cor_imp"),
  P50_MeasR = get_estimates(filter(plasboot, MeasName == "P50"), "MeasR_mean_cor", "MeasR_var_cor_imp"),
  Ptlp_MeasR = get_estimates(filter(plasboot, MeasName == "Ptlp"), "MeasR_mean_cor", "MeasR_var_cor_imp"),
  Ks_MeasR = get_estimates(filter(plasboot, MeasName == "Ks"), "MeasR_mean_cor", "MeasR_var_cor_imp"),
  KL_MeasR = get_estimates(filter(plasboot, MeasName == "KL"), "MeasR_mean_cor", "MeasR_var_cor_imp"),
  HV_MeasR = get_estimates(filter(plasboot, MeasName == "HV"), "MeasR_mean_cor", "MeasR_var_cor_imp")
)


# Step 5: Combine results for MeasT and MeasR into dataframes
combine_results <- function(Ypd, Ymd, P50, Ptlp, Ks, KL, HV, environment) {
  # Create individual dataframes for each MeasName
  df_Ypd <- data.frame(Combined = names(Ypd), Ypd = unname(Ypd))
  df_Ymd <- data.frame(Combined = names(Ymd), Ymd = unname(Ymd))
  df_P50 <- data.frame(Combined = names(P50), P50 = unname(P50))
  df_Ptlp <- data.frame(Combined = names(Ptlp), Ptlp = unname(Ptlp))
  df_Ks <- data.frame(Combined = names(Ks), Ks = unname(Ks))
  df_KL <- data.frame(Combined = names(KL), KL = unname(KL))
  df_HV <- data.frame(Combined = names(HV), HV = unname(HV))
  
  # Merge dataframes iteratively
  merged <- df_Ypd %>%
    full_join(df_Ymd, by = "Combined") %>%
    full_join(df_P50, by = "Combined") %>%
    full_join(df_Ptlp, by = "Combined") %>%
    full_join(df_Ks, by = "Combined") %>%
    full_join(df_KL, by = "Combined") %>%
    full_join(df_HV, by = "Combined") %>%
    mutate(Environment = environment)
  
  return(merged)
}

# Combine results for MeasT
final_MeasT <- combine_results(
  results$Ypd_MeasT, results$Ymd_MeasT, 
  results$P50_MeasT, results$Ptlp_MeasT, 
  results$Ks_MeasT, results$KL_MeasT, results$HV_MeasT, 
  "T"
)

# Combine results for MeasR
final_MeasR <- combine_results(
  results$Ypd_MeasR, results$Ymd_MeasR, 
  results$P50_MeasR, results$Ptlp_MeasR, 
  results$Ks_MeasR, results$KL_MeasR, results$HV_MeasR, 
  "R"
)

# View final dataframes
print("Final MeasT dataframe:")
print(final_MeasT)

print("Final MeasR dataframe:")
print(final_MeasR)


# Step 7: Combine both datasets into 'final'
final <- bind_rows(final_MeasT, final_MeasR)





#Regressions

reference <- subset(final, Environment == "R")
treatment <- subset(final, Environment == "T")
fit1 <- sma(log(abs(P50))~log(Ks+1), data=reference) #do it for each pair of traits
fit1 <- sma(log(abs(P50))~log(Ks+1), data=treatment) #do it for each pair of traits

summary(fit1)

#Plotting log P50/ logKs as an example

plot=ggplot(final, aes(x = log(abs(P50)), y = log(Ks+1), color=Environment))+
  geom_point(size=10,shape=21,stroke = 2)+
  stat_smooth(method = "lm", level=0.95)+
  xlab(expression(paste('log ',P[50])))+
  ylab(expression(paste('log ',K[s])))
plot
regression=plot+scale_color_manual(values=c("#0080ff", "red"))+
  geom_smooth(aes(color = Environment, fill = Environment), method = lm, size=3, alpha=0.1)+
  scale_fill_manual(values = c("#0080ff", "red"))+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(), panel.background = element_rect(fill = "white", colour = "black", size=2),
        axis.text.x = element_text(size = 58, color = "black", margin = margin(t = 5, r = 50, b = 20, l = 55)), 
        axis.text.y = element_text(size = 58, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 0)),
        axis.title.x = element_text(size = 68, color = "black", face="plain", vjust=2, margin = margin(t = 10, r = 0, b = 2, l = 0)),
        axis.title.y = element_text(size = 68, color = "black", face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x=element_line(color = "black", size=1),
        axis.ticks.y=element_line(color = "black",  size=1),
        axis.ticks.length = unit(0.25, "cm"), 
        legend.margin=margin(l = 0.6, unit='cm')) + 
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-0.5,3), ylim = c(0, 3)) + 
  scale_x_continuous(breaks=c(0,1,2,3)) + scale_y_continuous(breaks=c(0,1,2,3))

ggsave(filename = "regression.png", plot = regression, dpi = 300, width = 10, height = 10)
