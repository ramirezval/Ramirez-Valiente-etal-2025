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
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

######################################
##### Imputation procedure ###########
######################################

# Define geometric mean function ------------------------------------------

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Read master version of the data

hydplas.data.master <- read_tsv('HydPlasticity_final_2024-12-27.tsv')

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


#####################################
### Loading the imputated file  #####
#####################################

plas <- read_tsv('HydPlasticity_alltraits_imp_final_2024-12-27.tsv')
plas$SpeciesFix2 <- gsub(" ", "_", plas$SpeciesFix)

######################
######## MAP #########
######################
##### Fig. S1 #####
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
##### Fig. S2 #####
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



####################################################
### Important  first steps for selected traits #####
####################################################

#file with measurements of selected traits
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


####################################################
### Calculating effect sizes: MD and log-ratio #####
####################################################

plas$MeasR_SD_cor_imp<-plas$MeasR_SE_cor_imp*sqrt(plas$MeasR_N_cor_imp)
plas$MeasT_SD_cor_imp<-plas$MeasT_SE_cor_imp*sqrt(plas$MeasT_N_cor_imp)
plas <- subset(plas,plas$MeasR_SD_cor_imp>0) #deleting the negative value for MeasR_SD_imp
plas <- subset(plas,plas$MeasT_SD_cor_imp>0)#deleting the negative value for MeasT_SD_imp
plas$MeasR_var_cor_imp<-plas$MeasR_SD_cor_imp*plas$MeasR_SD_cor_imp
plas$MeasT_var_cor_imp<-plas$MeasT_SD_cor_imp*plas$MeasT_SD_cor_imp

#Estimate MD
plas<-escalc(measure="MD", n2i=plas$MeasR_N_cor_imp, n1i=plas$MeasT_N_cor_imp, m2i=plas$MeasR_mean_cor, m1i=plas$MeasT_mean_cor, sd2i=plas$MeasR_SD_cor_imp, sd1i=plas$MeasT_SD_cor_imp, 
             var.names=c("yi","vi"), add.measure=FALSE,
             append=TRUE, data=plas)
#Estimate log-ratio
plas<-escalc(measure="ROM", n2i=plas$MeasR_N_cor_imp, n1i=plas$MeasT_N_cor_imp, m2i=plas$MeasR_mean_cor, m1i=plas$MeasT_mean_cor, sd2i=plas$MeasR_SD_cor_imp, sd1i=plas$MeasT_SD_cor_imp, 
             var.names=c("yi","vi"), add.measure=FALSE,
             append=TRUE, data=plas)


#######################################
# Calculating variances of log-ratios #
#######################################

#Variance by MeasName for all studies
varall<-aggregate(yi ~ MeasName, plas, function(x) c(Var = var(x), Count = length(x)))
print(varall)

#Variance by MeasName considering only plasticity studies
plas.plas <- subset(plas, (plas$VarType %in% c('E')))
varplas<-aggregate(abs(yi) ~ MeasName, plas, function(x) c(Var = var(x), Count = length(x)))
print(varplas)

######################################################
# Testing for differences in variances of log-ratios #
######################################################
plas.traits <- subset(plas.plas, subset=(MeasName =="Ypd" |MeasName =="Ymd" |MeasName =="HV" |MeasName =="Ks" |MeasName =="KL"|MeasName =="Ptlp" |MeasName =="P50"))  #Selection 7 traits
fligner.test(yi ~ MeasName, data = plas.traits) #fligner test


################################################################################
#Density plots for log-ratios for individual traits y differences in variance ##
################################################################################

plas2<- subset(plas, subset=(VarTypePost =="E"))  #seleccionando solo plasticity studies
plastraits<- subset(plas2, subset=(MeasName =="P50" |MeasName =="Ptlp" |MeasName =="HV"|
                                     MeasName =="KL"| MeasName =="Ks" | MeasName =="Ymd" | MeasName =="Ypd"))  #seleccionado solo el factor water
##### Creating Fig. 1b #####
dens<-ggplot(plastraits, aes(x=yi, y=MeasName)) +
  geom_density_ridges(fill="lightblue",alpha=0.5)+scale_x_continuous(breaks = c(-1,0,1,2), limits= c(-1, 2))+
  geom_vline(xintercept=0, size=1, color="black",linetype="dashed")+ 
  coord_cartesian(ylim = c(1.5, 7.4)) 
dens
dens<-dens+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(plot.margin = unit(c(1,1.5,0.5,0.5), "cm"))+
  labs(title=, subtitle=, caption=, tag=,y ="", x = "log-ratio")+
  theme(legend.title = element_blank(),legend.position = "none",strip.background = element_blank(),strip.text.x = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(size = 29, color = "black", margin = margin(t = 15, r = 80, b = 20, l = 55)), 
        axis.text.y = element_text(size = 36, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 23), vjust=-0.6),
        axis.title.x = element_text(size = 36, color = "black", face="plain", vjust=2, margin = margin(t = 2, r = 0, b = 2, l = 0)),
        axis.ticks.x=element_line(color = "black"),
        axis.ticks.y=element_blank(),
        axis.ticks.length = unit(0.25, "cm"))+  scale_y_discrete(labels = c('Ymd' = expression(Psi[MD]),
                                                                            'Ypd'   = expression(Psi[PD]),
                                                                            'Ptlp'= expression(Pi[TLP]),
                                                                            'P50'= expression(P[50]),
                                                                            'Ks'= expression(K[S]),
                                                                            'KL'= expression(K[L]),
                                                                            'HV'= expression(Hv)))
print(dens)
ggsave(filename = "Fig.1a.png", plot = dens, dpi = 300, width = 8, height = 11)



################################################################################
# Meta-analysis for absolute log-ratios for individual traits across factors ###
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


##### Creating Fig. 1a ########

#Load the file with the absolute log-ratios created in the previous section
abslogratio <- read_tsv('abslogratio.tsv')

#Variable selection
abslogratio$MeasName <- factor(abslogratio$MeasName, levels = c("HV", "KL", "Ks", "P50", "Ptlp", "Ymd", "Ypd"))

# Add the "Significant" column to the dataframe
# It is defined as significant if the confidence interval does not include the value 0
abslogratio$Significant <- ifelse(abslogratio$Lower_CI > 0 | abslogratio$Upper_CI < 0, TRUE, FALSE)
#abslogratio$MeasName <- factor(abslogratio$MeasName,levels = rev(c("Hv", "KL", "Ks", "P50", "Ptlp","Ymd", "Ypd")))
abslogratio$Subset <- factor(abslogratio$Subset,levels = rev(c("No Spatial","All Studies","Only Spatial")))

# Create the plot with ggplot and inverted axes
p <- ggplot(abslogratio, aes(y = MeasName, x = Estimate, fill = factor(Subset))) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0, position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(shape = 21, size = 8, color = "black", stroke = 1, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("grey", "white", "black")) + coord_cartesian(ylim = c(1.1, 6.9))+
  labs(y = "", x = "|log-ratio|") +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text.y = element_text(size = 38, color = "black", margin = margin(t = 10, r=15)),  # Space between letters and ticks on the y-axis
        axis.text.x = element_text(size = 35, color = "black", margin = margin(r = 5, t=16)),  # Space between numbers and ticks on the x-axis
        axis.title.x = element_text(size = 42, color = "black", margin = margin(r = 10, t=15)),  # Space between the x-axis title and numbers
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.x=element_line(color = "black", size = 0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.length = unit(0.30, "cm")) +
  scale_x_continuous(limits = c(0, 0.8), breaks = seq(-0.4, 1, by = 0.2)) +
  # Add text with sample size if wanted
  #geom_text(aes(label = Sample_Size, x = Upper_CI + 0.05), position = position_dodge(width = 0.8), size = 7, color = "black", hjust = 0) +
  scale_y_discrete(labels = c("Hv", expression(K[L]), expression(K[S]), expression(P[50]), 
                              expression(Π[TLP]), expression(Ψ[MD]), expression(Ψ[PD])))
p

# Save the plot with ggsave at 300 dpi
ggsave(filename = "Fig.1b.png", plot = p, dpi = 300, width = 8, height = 12)


##### Creating Fig. S3a #####
#Load the file with the absolute log-ratios created in the previous section
abslogratio <- read_tsv('abslogratio.tsv')

# two subsets: "No Spatial" and "Spatial" using subset
no_spatial <- subset(abslogratio, Subset == "No Spatial", select = c(MeasName, Estimate, Lower_CI, Upper_CI))
spatial <- subset(abslogratio, Subset == "Only Spatial", select = c(MeasName, Estimate, Lower_CI, Upper_CI))

# Merge the two subsets by MeasName
merged_data <- merge(no_spatial, spatial, by = "MeasName", suffixes = c("_NoSpatial", "_Spatial"))

# Perform linear regression
model <- lm(Estimate_Spatial ~ Estimate_NoSpatial, data = merged_data)

r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, 4]

sma_result <- sma(Estimate_Spatial ~ Estimate_NoSpatial, data = merged_data)
summary(sma_result)

# Plot the data with regression line and confidence intervals
plot <- ggplot(merged_data, aes(x = Estimate_NoSpatial, y = Estimate_Spatial)) +
  geom_point(size = 7) +
  geom_errorbar(aes(ymin = Lower_CI_Spatial, ymax = Upper_CI_Spatial), width = 0.02, color = "black") +
  geom_errorbarh(aes(xmin = Lower_CI_NoSpatial, xmax = Upper_CI_NoSpatial), height = 0.02, color = "black") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", fill = "black", alpha=0.04)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1)
  
figs3a<-plot+ theme(panel.grid = element_blank(),plot.margin = unit(c(1,0.5,0.5,0.5), "cm"))+
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(), panel.background = element_rect(fill = "white", colour = "black", size=2),
        axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 5, r = 50, b = 20, l = 55)), 
        axis.text.y = element_text(size = 45, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 0)),
        axis.title.x = element_text(size = 47, color = "black", face="plain", vjust=2, margin = margin(t = 10, r = 0, b = 2, l = 0)),
        axis.title.y = element_text(size = 47, color = "black", face="plain", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.ticks.x=element_line(color = "black", size=1),
        axis.ticks.y=element_line(color = "black",  size=1),
        axis.ticks.length = unit(0.25, "cm"), 
        legend.margin=margin(l = 0.6, unit='cm'))+
  coord_cartesian(xlim = c(0,0.8), ylim = c(0,0.8)) + 
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8)) + scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8))+
  xlab("|log-ratio| Plasticity studies")+ylab("|log-ratio| Spatial studies")
figs3a
ggsave(filename = "Fig.S3a.png", plot = figs3a, dpi = 300, width = 10, height = 10)



##################################################################
# Meta-analysis for absolute log-ratios for experimental types ###
##################################################################

# Function to perform an analysis with bootstrap for a specific VarContextAgg level
run_analysis_with_bootstrap_varcontext <- function(data, varcontext, tot_rep = 100) {
  
  # Filter the dataset for the specific VarContextAgg level
  subset_data <- subset(data, VarContextAgg == varcontext)
  
  # Fit a random-effects meta-analysis model using rma.mv
  model <- rma.mv(abs(yi), vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2, ~ 1 | MeasName), 
                  control = list(optimizer = "optim", optmethod = "Nelder-Mead"), 
                  data = subset_data)
  
  # Define the bootstrap function
  boot_func <- function(plasboot, indices, prog) {
    
    # Update progress bar with each iteration
    prog$tick()
    
    sel <- plasboot[indices, ]
    res <- try(suppressWarnings(rma.mv(abs(yi), vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2, ~ 1 | MeasName), data = sel)), silent = TRUE)
    
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

# Function to run the analysis for each level of VarContextAgg
run_analysis_by_varcontext <- function(data, tot_rep = 100) {
  
  # Get unique levels of VarContextAgg
  varcontext_levels <- unique(data$VarContextAgg)
  
  # Initialize a list to store results
  all_results <- list()
  
  # Loop through each level of VarContextAgg
  for (varcontext in varcontext_levels) {
    cat("Analyzing VarContextAgg level:", varcontext, "\n")  # Print progress
    
    # Run the analysis for the current VarContextAgg level
    res <- run_analysis_with_bootstrap_varcontext(data, varcontext, tot_rep)
    
    # Extract specific values from the result
    extracted_values <- res$extracted_values
    
    # Add VarContextAgg level to the extracted values
    extracted_values$VarContextAgg <- varcontext
    
    # Create a data frame with the extracted values
    results_df <- data.frame(
      VarContextAgg = extracted_values$VarContextAgg,
      Estimate = extracted_values$estimate,
      Lower_CI = extracted_values$bca_interval[1],
      Upper_CI = extracted_values$bca_interval[2],
      Sample_Size = extracted_values$sample_size
    )
    
    # Add the result to the list of all results
    all_results[[varcontext]] <- results_df
  }
  
  # Combine all results into a single data frame
  final_results_df <- do.call(rbind, all_results)
  
  return(final_results_df)
}

# Call the function to create the dataset with all levels of VarContextAgg
results_by_varcontext_df <- run_analysis_by_varcontext(plas, tot_rep = 100)  # 'plas' is the dataset for analysis

# Display the final dataset
print(results_by_varcontext_df)

# Save the results_by_varcontext_df dataset as a .tsv file
write.table(abslogratiovarcontext, 
            file = "abslogratiovarcontext.tsv",  # Specify the output file name
            sep = "\t",                           # Use tab as the separator
            row.names = FALSE,                    # Do not include row names
            quote = FALSE)

##### Creating Fig. S3b #####
logratiofactors <- read_tsv('abslogratiovarcontext.tsv')



#######################################################################################################
######### Meta-analysis for log ratios for each combination of trait and environmental factor #########
#######################################################################################################

# Function to perform bootstrap analysis for all combinations of traits and factors
run_analysis_with_bootstrap <- function(data, trait, factor_level, tot_rep = 10000, min_sample_size = 3) {
  
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
  
  # Attempt to fit the random-effects meta-analysis model using rma.mv
  model <- tryCatch(
    {
      rma.mv(yi, vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2), 
             control = list(optimizer = "optim", optmethod = "Nelder-Mead"), 
             data = subset_data)
    },
    error = function(e) NULL
  )
  
  # If the first attempt fails, retry with a different optimizer
  if (is.null(model)) {
    cat("  Warning: Initial model fitting failed for Trait:", trait, "and Factor:", factor_level, "\n")
    model <- tryCatch(
      {
        rma.mv(yi, vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2), 
               control = list(optimizer = "nlminb"), 
               data = subset_data)
      },
      error = function(e) NULL
    )
  }
  
  # If model fitting still fails, return NA values
  if (is.null(model)) {
    cat("  Error: Model fitting failed for Trait:", trait, "and Factor:", factor_level, "even after retrying\n")
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Bootstrap function
  boot_func <- function(plasboot, indices, prog) {
    prog$tick()
    sel <- plasboot[indices, ]
    res <- try(suppressWarnings(rma.mv(yi, vi, random = list(~ 1 | SpeciesFix, ~ 1 | PaperID), data = sel)), silent = TRUE)
    if (inherits(res, "try-error")) NA else c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
  # Initialize the progress bar
  pb <- progress_bar$new(total = tot_rep + 1)
  
  # Run the bootstrap
  boot_results <- tryCatch(
    boot(subset_data, boot_func, tot_rep, prog = pb),
    error = function(e) NULL
  )
  
  # If bootstrap fails, return NA values
  if (is.null(boot_results)) {
    cat("  Error: Bootstrap failed for Trait:", trait, "and Factor:", factor_level, "\n")
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Calculate confidence intervals
  ci_results <- tryCatch(
    boot.ci(boot_results, type = "bca", index = 1:2),
    error = function(e) NULL
  )
  
  # If confidence interval calculation fails, return NA values
  if (is.null(ci_results)) {
    cat("  Error: Confidence interval calculation failed for Trait:", trait, "and Factor:", factor_level, "\n")
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
run_analysis_with_subsets <- function(data, tot_rep = 10000) {
  
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
logratiofactors <- run_analysis_with_subsets(plas, tot_rep = 10000)

# Display the final dataset
print(logratiofactors)

# Save the results as a .tsv file
write.table(logratiofactors, 
            file = "logratiofactors.tsv",  # Specify the output file name
            sep = "\t",                   # Use tab as the separator
            row.names = FALSE,            # Do not include row names
            quote = FALSE)                # Do not add quotes

##### Creating Fig. 2 ########
#Load the file with the log-ratios created in the previous section
logratiofactors <- read_tsv('logratiofactors.tsv')

#Plot log-ratios only for plasticity and all_studies subsets: datafigure2<- subset(logratiofactors, subset=(VarFactorAgg =="C"))
logratiofactors.2<-subset(logratiofactors, subset=(Subset =="No Spatial" |Subset =="All Studies")) #C,N,W,T_VPD,L

#Figures are created for each environmental factor separately
datafigure2<-subset(logratiofactors.2, subset=(VarFactorAgg =="C")) #C,N,W,T_VPD,L
datafigure2$MeasName <- factor(datafigure2$MeasName,levels = rev(c("Hv", "KL", "Ks", "P50", "Ptlp","Ymd", "Ypd")))
datafigure2$Subset <- factor(datafigure2$Subset,levels = rev(c("All Studies", "No Spatial")))
logratio <- ggplot(datafigure2, aes(x = MeasName, y = Estimate, fill = factor(Subset))) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0, position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(shape = 21, size = 9, color = "black", stroke = 1, position = position_dodge(width = 0.6))+
  scale_fill_manual(values = c("#9C9C9C","white"))+ #"T_VPD"="#FF0000", "C"="#9C9C9C", "N"="#00c800", "L"="#cccc00", "W"="#008DF2")
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "", y = "log-ratio") +
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(size = 32, color = "black", margin = margin(t = 10)),  # Espacio entre letras y ticks en el eje x
        axis.text.y = element_text(size = 30, color = "black", margin = margin(r = 5)),  # Espacio entre números y ticks en el eje y
        axis.title.y = element_text(size = 38, color = "black", margin = margin(r = 10)),  # Espacio entre el título del eje y y los números
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.15, "cm")) +
  scale_y_continuous(limits = c(-0.26, 1), breaks = seq(-0.4, 1, by = 0.2))+
  scale_x_discrete(labels = c(expression(Ψ[PD]),expression(Ψ[MD]),expression(Π[TLP]),expression(P[50]),
                              expression(K[S]),expression(K[L]),"Hv"))
logratio
ggsave(filename = "Fig2.C.png", plot = logratio, dpi = 300, width = 10, height = 5.5)

##### Creating Fig. S4 ########
#Load the file with the log-ratios created in the previous section
logratiofactors <- read_tsv('logratiofactors.tsv')

#Plot log-ratios for the three subsets
datafigure2<-subset(logratiofactors, subset=(VarFactorAgg =="T_VPD")) #C,N,W,T_VPD,L
datafigure2$MeasName <- factor(datafigure2$MeasName,levels = rev(c("Hv", "KL", "Ks", "P50", "Ptlp","Ymd", "Ypd")))
datafigure2$Subset <- factor(datafigure2$Subset,levels = rev(c("Only Spatial","All Studies","No Spatial")))
logratio3 <- ggplot(datafigure2, aes(x = MeasName, y = Estimate, fill = factor(Subset))) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0, position = position_dodge(width = 0.75), size = 0.5) +
  geom_point(shape = 21, size = 9, color = "black", stroke = 1, position = position_dodge(width = 0.75))+
  scale_fill_manual(values = c("#FF0000","white","#ffcccb"))+ #"T_VPD"="#FF0000", "C"="#9C9C9C", "N"="#00c800", "L"="#cccc00", "W"="#008DF2")
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "", y = "log-ratio") +
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(size = 32, color = "black", margin = margin(t = 10)),  # Espacio entre letras y ticks en el eje x
        axis.text.y = element_text(size = 30, color = "black", margin = margin(r = 5)),  # Espacio entre números y ticks en el eje y
        axis.title.y = element_text(size = 38, color = "black", margin = margin(r = 10)),  # Espacio entre el título del eje y y los números
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.15, "cm")) +
  scale_y_continuous(limits = c(-0.42, 1), breaks = seq(-0.4, 1, by = 0.2))+
  scale_x_discrete(labels = c(expression(Ψ[PD]),expression(Ψ[MD]),expression(Π[TLP]),expression(P[50]),
                              expression(K[S]),expression(K[L]),"Hv"))
logratio3

ggsave(filename = "Fig4S.T.png", plot = logratio3, dpi = 300, width = 12.5, height = 5.5)


############################################################
##### Meta-analysis for MD for Ypd, Ymd, Ptlp, P50 and HSM #
############################################################
#1. Estimate MD ####
plas<-escalc(measure="MD", n2i=plas$MeasR_N_cor_imp, n1i=plas$MeasT_N_cor_imp, m2i=plas$MeasR_mean_cor, m1i=plas$MeasT_mean_cor, sd2i=plas$MeasR_SD_cor_imp, sd1i=plas$MeasT_SD_cor_imp, 
             var.names=c("yi","vi"), add.measure=FALSE,
             append=TRUE, data=plas)

#2. Meta-analysis for Ypd, Ymd, Ptlp and P50 ####

# Function to perform bootstrap analysis for all combinations of traits and factors
run_analysis_with_bootstrap <- function(data, trait, factor_level, tot_rep = 10000, min_sample_size = 3) {
  
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
  
  # Attempt to fit the random-effects meta-analysis model using rma.mv
  model <- tryCatch(
    {
      rma.mv(yi, vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2), 
             control = list(optimizer = "optim", optmethod = "Nelder-Mead"), 
             data = subset_data)
    },
    error = function(e) NULL
  )
  
  # If the first attempt fails, retry with a different optimizer
  if (is.null(model)) {
    cat("  Warning: Initial model fitting failed for Trait:", trait, "and Factor:", factor_level, "\n")
    model <- tryCatch(
      {
        rma.mv(yi, vi, random = list(~ 1 | PaperID, ~ 1 | SpeciesFix2), 
               control = list(optimizer = "nlminb"), 
               data = subset_data)
      },
      error = function(e) NULL
    )
  }
  
  # If model fitting still fails, return NA values
  if (is.null(model)) {
    cat("  Error: Model fitting failed for Trait:", trait, "and Factor:", factor_level, "even after retrying\n")
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Extract specific parameters
  estimate <- coef(model)[1]
  
  # Bootstrap function
  boot_func <- function(plasboot, indices, prog) {
    prog$tick()
    sel <- plasboot[indices, ]
    res <- try(suppressWarnings(rma.mv(yi, vi, random = list(~ 1 | SpeciesFix, ~ 1 | PaperID), data = sel)), silent = TRUE)
    if (inherits(res, "try-error")) NA else c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
  # Initialize the progress bar
  pb <- progress_bar$new(total = tot_rep + 1)
  
  # Run the bootstrap
  boot_results <- tryCatch(
    boot(subset_data, boot_func, tot_rep, prog = pb),
    error = function(e) NULL
  )
  
  # If bootstrap fails, return NA values
  if (is.null(boot_results)) {
    cat("  Error: Bootstrap failed for Trait:", trait, "and Factor:", factor_level, "\n")
    return(data.frame(
      Estimate = NA, Lower_CI = NA, Upper_CI = NA, Sample_Size = nrow(subset_data)
    ))
  }
  
  # Calculate confidence intervals
  ci_results <- tryCatch(
    boot.ci(boot_results, type = "bca", index = 1:2),
    error = function(e) NULL
  )
  
  # If confidence interval calculation fails, return NA values
  if (is.null(ci_results)) {
    cat("  Error: Confidence interval calculation failed for Trait:", trait, "and Factor:", factor_level, "\n")
    return(data.frame(
      Estimate = estimate,
      Lower_CI = NA,
      Upper_CI = NA,
      Sample_Size = model$k
    ))
  }
  
  # Return the results
  data.frame(
    Estimate = estimate,
    Lower_CI = ci_results$bca[4],
    Upper_CI = ci_results$bca[5],
    Sample_Size = model$k
  )
}

# Function to run analysis for each combination of Trait and Factor level
run_analysis_with_subsets <- function(data, tot_rep = 10000) {
  
  # Define subsets
  subsets <- list(
    "No Spatial" = subset(data, VarContextAgg != "spatial"),
    "Only Spatial" = subset(data, VarContextAgg == "spatial"),
    "All Studies" = data
  )
  
  all_results <- list()
  
  # Filter for specific MeasName and VarFactorAgg
  filtered_data <- subset(data, MeasName %in% c("Ypd", "Ymd", "Ptlp", "P50") & VarFactorAgg == "W")
  
  for (subset_name in names(subsets)) {
    cat("Analyzing subset:", subset_name, "\n")
    subset_data <- subset(subsets[[subset_name]], MeasName %in% c("Ypd", "Ymd", "Ptlp", "P50") & VarFactorAgg == "W")
    
    # Get unique combinations of MeasName and VarFactorAgg
    meas_levels <- unique(filtered_data$MeasName)
    factor_levels <- unique(filtered_data$VarFactorAgg)
    
    for (trait in meas_levels) {
      for (factor_level in factor_levels) {
        
        cat("  Analyzing Trait:", trait, "and Factor:", factor_level, "\n")
        
        # Run the analysis for each combination
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

# Run the function with the filtered data
mdwaterpotentials <- run_analysis_with_subsets(plas, tot_rep = 10000)

# Display the final dataset
print(mdwaterpotentials)

# Save the results as a .tsv file
write.table(mdwaterpotentials, 
            file = "mdwaterpotentials.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


#3. Estimate MD of global Hydraulic safety margin (HSM) as the difference in the bootstrap distribution of P50 and Ymd ####

plas<-escalc(measure="MD", n2i=plas$MeasR_N_cor_imp, n1i=plas$MeasT_N_cor_imp, m2i=plas$MeasR_mean_cor, m1i=plas$MeasT_mean_cor, sd2i=plas$MeasR_SD_cor_imp, sd1i=plas$MeasT_SD_cor_imp, 
             var.names=c("yi","vi"), add.measure=FALSE,
             append=TRUE, data=plas)


plas.w<- subset(plas, subset=(VarFactorAgg =="W"))
plas.plas <- subset(plas.w, subset = (VarContextAgg != "spatial"))
plas.spatial <- subset(plas.w, subset = (VarContextAgg == "spatial"))
plas.all <- plas.w  # Include all studies


plasboot<-plas.plas.P50<- subset(plas.plas, subset=(MeasName =="P50"))  
plasboot<-plas.plas.Ymd<- subset(plas.plas, subset=(MeasName =="Ymd"))
plasboot<-plas.spatial.P50<- subset(plas.spatial, subset=(MeasName =="P50"))  
plasboot<-plas.spatial.Ymd<- subset(plas.spatial, subset=(MeasName =="Ymd"))
plasboot<-plas.all.P50<- subset(plas.all, subset=(MeasName =="P50"))  
plasboot<-plas.all.Ymd<- subset(plas.all, subset=(MeasName =="Ymd"))

#Analyses for obtain md for each trait

boot.func <- function(plasboot, indices, prog) {
  
  #display progress with each run of the function
  prog$tick()
  
  sel <- plasboot[indices,]
  res <- try(suppressWarnings(rma.mv(yi, vi, random = list(~ 1 | SpeciesFix, ~ 1 | PaperID), data=sel)), silent=TRUE)
  
  if (inherits(res, "try-error")) {
    NA
  } else {
    c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
}

tot_rep <- 10000

#initialize progress bar objec
pb <- progress_bar$new(total = tot_rep + 1) 

boot.plas.plas.P50<- boot(plas.plas.P50, boot.func, tot_rep, prog = pb)
boot.plas.plas.Ymd<- boot(plas.plas.Ymd, boot.func, tot_rep, prog = pb)
boot.plas.spatial.P50<- boot(plas.spatial.P50, boot.func, tot_rep, prog = pb)
boot.plas.spatial.Ymd<- boot(plas.spatial.Ymd, boot.func, tot_rep, prog = pb)
boot.plas.all.P50<- boot(plas.all.P50, boot.func, tot_rep, prog = pb)
boot.plas.all.Ymd<- boot(plas.all.Ymd, boot.func, tot_rep, prog = pb)

#Plasticity HSM
p50<-boot.plas.plas.P50[["t"]]
Ymd<-boot.plas.plas.Ymd[["t"]]
write.csv(p50, file = "plas.plas.P50.csv")
write.csv(Ymd, file = "plas.plas.Ymd.csv")
plas.plas.P50<-read.table("plas.plas.P50.csv", header=T, sep=",", fill = TRUE)
plas.plas.Ymd<-read.table("plas.plas.Ymd.csv", header=T, sep=",", fill = TRUE)
diffplas=plas.plas.Ymd$V1-plas.plas.P50$V1
hist(diffplas)
median(diffplas, na.rm=TRUE)
quantile(diffplas,c(0.05), na.rm=TRUE)
quantile(diffplas,c(0.95), na.rm=TRUE)

#All studies HSM
p50<-boot.plas.all.P50[["t"]]
Ymd<-boot.plas.all.Ymd[["t"]]
write.csv(p50, file = "plas.all.P50.csv")
write.csv(Ymd, file = "plas.all.Ymd.csv")
plas.all.P50<-read.table("plas.all.P50.csv", header=T, sep=",", fill = TRUE)
plas.all.Ymd<-read.table("plas.all.Ymd.csv", header=T, sep=",", fill = TRUE)
diffall=plas.all.Ymd$V1-plas.all.P50$V1
hist(diffall)
median(diffall, na.rm=TRUE)
quantile(diffall,c(0.05), na.rm=TRUE)
quantile(diffall,c(0.95), na.rm=TRUE)

#Spatial HSM
p50<-boot.plas.spatial.P50[["t"]]
Ymd<-boot.plas.spatial.Ymd[["t"]]
write.csv(p50, file = "plas.spatial.P50.csv")
write.csv(Ymd, file = "plas.spatial.Ymd.csv")
plas.spatial.P50<-read.table("plas.spatial.P50.csv", header=T, sep=",", fill = TRUE)
plas.spatial.Ymd<-read.table("plas.spatial.Ymd.csv", header=T, sep=",", fill = TRUE)
diffspatial=plas.spatial.Ymd$V1-plas.spatial.P50$V1
hist(diffspatial)
median(diffspatial, na.rm=TRUE)
quantile(diffspatial,c(0.05), na.rm=TRUE)
quantile(diffspatial,c(0.95), na.rm=TRUE)


# Combine results into a table
results_table_HSM <- data.frame(
  Subset = c("No Spatial", "Only Spatial", "All Studies"),
  MeasName = "HSM",
  Estimate = c(median(diffplas, na.rm=TRUE),
               median(diffspatial, na.rm=TRUE),
               median(diffall, na.rm=TRUE)),
  Lower_CI = c(quantile(diffplas, c(0.05), na.rm=TRUE),
               quantile(diffspatial, c(0.05), na.rm=TRUE),
               quantile(diffall, c(0.05), na.rm=TRUE)),
  Upper_CI = c(quantile(diffplas, c(0.95), na.rm=TRUE),
               quantile(diffspatial, c(0.95), na.rm=TRUE),
               quantile(diffall, c(0.95), na.rm=TRUE)),
  Sample_Size = NA
)

# Display the results table
print(results_table_HSM)

write.table(results_table_HSM, 
            file = "results_table_HSM.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

#4. Combine results from 2 and 3 ####
results_table_HSM <- read_tsv('results_table_HSM.tsv')
mdwaterpotentials <- read_tsv('mdwaterpotentials.tsv')
# Ensure both tables have the same column names
names(results_table_HSM) <- c("Subset", "MeasName", "Estimate", "Lower_CI", "Upper_CI", "Sample_Size")

# For mdwaterpotentials, ensure the relevant columns exist and match the same order
mdwaterpotentials <- mdwaterpotentials[, c("Subset", "MeasName", "Estimate", "Lower_CI", "Upper_CI", "Sample_Size")]

# Combine both tables into one
waterpotentials_HSM <- rbind(mdwaterpotentials, results_table_HSM)

# Display the combined results
print(waterpotentials_HSM)

# Save the combined results to a file
write.table(waterpotentials_HSM, 
            file = "waterpotentials_HSM.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

##### Creating Fig. 3a ##########
#Load the file with the md created in the previous section
waterpotentials_HSM <- read_tsv('waterpotentials_HSM.tsv')
waterpotentials_HSM$Subset <- factor(waterpotentials_HSM$Subset,levels = rev(c("Only Spatial", "All Studies", "No Spatial")))
#Plot mds only for plasticity and all_studies subsets: waterpotentials_HSM<- subset(waterpotentials_HSM))
waterpotentials_HSM2<-subset(waterpotentials_HSM, subset=(Subset =="No Spatial" |Subset =="All Studies")) #C,N,W,T_VPD,L
waterpotentials_HSM2$MeasName <- factor(waterpotentials_HSM2$MeasName,levels = rev(c("HSM","P50", "Ptlp","Ymd", "Ypd")))
md2 <- ggplot(waterpotentials_HSM2, aes(x = MeasName, y = Estimate, fill = factor(Subset))) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0, position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(shape = 21, size = 9, color = "black", stroke = 1, position = position_dodge(width = 0.6))+
  scale_fill_manual(values = c("#008DF2","white"))+
  #geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "", y = expression("Adjustment " * Delta * " (MPa)")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(size = 28, color = "black", margin = margin(t = 10)),  # Espacio entre letras y ticks en el eje x
        axis.text.y = element_text(size = 30, color = "black", margin = margin(r = 5)),  # Espacio entre números y ticks en el eje y
        axis.title.y = element_text(size = 37, color = "black", margin = margin(r = 10)),  # Espacio entre el título del eje y y los números
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.15, "cm")) +
  scale_y_continuous(limits = c(-1, 0), breaks = seq(-1, 0, by = 0.2))+
  scale_x_discrete(labels = c(expression(Ψ[PD]),expression(Ψ[MD]),expression(Π[TLP]),expression(P[50]),
                              "HSM"))
md2
ggsave(filename = "Fig3a.png", plot = md2, dpi = 300, width = 8, height = 6)

##### Creating Fig. S5 ##########
#Load the file with the md created in the previous section
waterpotentials_HSM <- read_tsv('waterpotentials_HSM.tsv')
waterpotentials_HSM$Subset <- factor(waterpotentials_HSM$Subset,levels = rev(c("Only Spatial", "All Studies", "No Spatial")))
#Plot mds for all subsets
waterpotentials_HSM$MeasName <- factor(waterpotentials_HSM$MeasName,levels = rev(c("HSM","P50", "Ptlp","Ymd", "Ypd")))
md3 <- ggplot(waterpotentials_HSM, aes(x = MeasName, y = Estimate, fill = factor(Subset))) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0, position = position_dodge(width = 0.6), size = 0.5) +
  geom_point(shape = 21, size = 9, color = "black", stroke = 1, position = position_dodge(width = 0.6))+
  scale_fill_manual(values = c("#008DF2","white", "lightblue"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "", y = expression("Adjustment " * Delta * " (MPa)")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(size = 28, color = "black", margin = margin(t = 10)),  # Espacio entre letras y ticks en el eje x
        axis.text.y = element_text(size = 30, color = "black", margin = margin(r = 5)),  # Espacio entre números y ticks en el eje y
        axis.title.y = element_text(size = 38, color = "black", margin = margin(r = 10)),  # Espacio entre el título del eje y y los números
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.15, "cm")) +
  scale_y_continuous(limits = c(-1, 0.3), breaks = seq(-1, 0.3, by = 0.2))+
  scale_x_discrete(labels = c(expression(Ψ[PD]),expression(Ψ[MD]),expression(Π[TLP]),expression(P[50]),
                              "HSM"))
md3
ggsave(filename = "FigS5.png", plot = md3, dpi = 300, width = 11, height = 6)



###################################################################################
##### Testing for differences in HSM between reference and treatment environments #
###################################################################################

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
HSM_two <- bind_rows(final_MeasT, final_MeasR)

# Step 9: Add 'HSM' as the difference between Ymd and P50
HSM_two <- final %>%
  mutate(HSM = if_else(!is.na(Ymd) & !is.na(P50), Ymd - P50, NA_real_))

# Step 10: Save the results
write.table(HSM_two, 
            file = "HSM_two_environments.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


# Display the final dataframe
print(HSM_two)


# Modelo con lmer
test <- lmer(HSM ~ Environment + (1|ID) + (1|Species) + (1|VarContextAgg), data = HSM_two)
Anova(test)
summary(test)
emmeans(test, pairwise~Environment)


##### Creating Fig. 3b ######
HSM_two <- read_tsv('HSM_two_environments.tsv')
ggplot(HSM_two, aes(x = Environment, y = HSM))+ geom_boxplot(aes(fill= Environment))
plot <- ggplot(HSM_two, aes(x = factor(Environment), y = HSM, color = factor(Environment))) + 
  geom_boxplot(width = 2, lwd = 1, outlier.shape = NA) +  # Remove outlier points
  scale_color_manual(values = c("#00B4F7", "#FF0000")) +
  geom_point(position = "jitter", size = 4, shape = 20, stroke = 1) +
  stat_summary(fun = mean, geom = "crossbar", size = 1, color = "black", fill = "red", linetype = "dashed") +
  scale_x_discrete(labels = c("R" = "Reference", "T" = "Treatment")) + # Replace labels
  scale_y_continuous(
    limits = c(-2.25, 6.5),         # Set the limits of the y-axis
    breaks = seq(-4, 7, 2)         # Define tick intervals
  ) +
  labs(x = NULL, y = expression(paste('HSM ', (MPa)))) +  # Remove x-axis label
  theme_bw() +
  theme(
    plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black", size = 2),
    axis.text.x = element_text(size = 44, color = "black", margin = margin(t = 16, r = 50, b = 20, l = 55)), # Adjust x-axis tick text size
    axis.text.y = element_text(size = 52, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 10)),
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 58, color = "black", face = "plain", margin = margin(t = 0, r = 1, b = 0, l = 0)),
    axis.ticks.x = element_line(color = "black", size = 1),
    axis.ticks.y = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.25, "cm"), 
    legend.margin = margin(l = 0.6, unit = 'cm'),
    legend.position = "none"
  )
plot

ggsave(filename = "Fig3b.png", plot = plot, dpi = 300, width = 9.2, height = 10)

###############################################################################################
############### Correlations of log ratios by paper and environmental factor ##################
###############################################################################################

# Create the 'Combined' variable in the original dataframe
plas$Combined <- str_c(plas$PaperID, "_", plas$SpeciesFix2, "_", plas$VarFactorAgg, "_", plas$VarContextAgg, "_", plas$MeasName)

# Perform meta-analysis
res <- rma(yi, vi, mods = ~Combined - 1, 
           control = list(optimizer="optim", optmethod="Nelder-Mead"), data = plas)
res

# Create a dataframe with the coefficients from the meta-analysis
data <- data.frame(Moderator = names(coef(res)), Estimate = coef(res))

# Clean the 'Moderator' column
data$Moderator <- gsub("Charra-Vaskou_ 2012", "Charra-Vaskou_2012", data$Moderator)
data$Moderator <- gsub("CombinedEpron et al_", "CombinedEpron_2015_TreePhysiol_", data$Moderator)
data$Moderator <- gsub(" ", "_", data$Moderator)  # Replace spaces with underscores
data$Moderator <- gsub(" ", "", data$Moderator)  # Remove spaces
data$Moderator <- gsub("_et_al", "", data$Moderator)  # Remove "_et_al"
data$Moderator <- gsub("CombinedCarins_Murphy", "CombinedCarins-Murphy", data$Moderator)
data$Moderator <- gsub("CombinedLamy_Pinus_pinaster", "CombinedLamy_2014_NewPhytol_Pinus_pinaster", data$Moderator)
data$Moderator <- gsub("_T_Phys_", "_TPhys_", data$Moderator)
data$Moderator <- gsub("Plant_soil", "Plantsoil", data$Moderator)
data$Moderator <- gsub("2014_plant_biology", "2014_plantbiology", data$Moderator)
data$Moderator <- gsub("Petruzzellis_2019_Tree_Physiology", "Petruzzellis_2019_TreePhysiology", data$Moderator)
data$Moderator <- gsub("Renninger_2007_", "Renninger_2007_TreePhysiol_", data$Moderator)
data$Moderator <- gsub("Stiller_2009_Environ._Exp._Bot._", "Stiller_2009_EnvironExpBot_", data$Moderator)
data$Moderator <- gsub("Tree_P_", "TreeP_", data$Moderator)

# Adjust "T_VPD" values manually to treat them as a single entity in VarFactorAgg
data$Moderator <- gsub("_T_VPD_", "_TVPD_", data$Moderator)

# Split the 'Moderator' column into parts separated by '_'
split_moderator <- strsplit(data$Moderator, "_")

# Extract columns from 'split_moderator'

# Extract 'PaperID' (combine values 1, 2, and 3, removing "Combined")
data$PaperID <- sapply(split_moderator, function(x) {
  paste(x[1:3], collapse = "_")
})
data$PaperID <- gsub("^Combined", "", data$PaperID)  # Remove "Combined"

# Extract 'Species' (combine values 4 and 5)
data$Species <- sapply(split_moderator, function(x) {
  paste(x[4:5], collapse = "_")
})

# Extract 'VarFactorAgg' (value 6, previously adjusted for "T_VPD" as "TVPD")
data$VarFactorAgg <- sapply(split_moderator, function(x) {
  x[6]
})

# Extract 'VarContextAgg' (value 7)
data$VarContextAgg <- sapply(split_moderator, function(x) {
  x[7]
})

# Extract 'MeasName' (last value)
data$MeasName <- sapply(split_moderator, function(x) {
  tail(x, 1)
})

# Verify the first rows of the dataframe with the new columns
unique(data$VarFactorAgg)
unique(data$VarContextAgg)
unique(data$MeasName)

# Reshape the dataframe to wide format
data <- data %>%
  pivot_wider(
    id_cols = c(VarFactorAgg, VarContextAgg, Species, PaperID), # Grouping variables
    names_from = MeasName,  # Column names will be levels of MeasName
    values_from = Estimate,  # Values for the new columns come from Estimate
    names_glue = "{MeasName}_yi"  # Add the "_yi" suffix to the column names
  )

# Check the new dataframe
head(data)

# Write the new dataframe to a file
write.table(data, 
            file = "fileforcorrelations.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Optionally, read the file back
data <- read_tsv('fileforcorrelations.tsv')
data.l<- subset(data, subset=(VarFactorAgg =="L"))
data.n<- subset(data, subset=(VarFactorAgg =="N"))
data.w<- subset(data, subset=(VarFactorAgg =="W"))
data.c<- subset(data, subset=(VarFactorAgg =="C"))
data.w.plas<- subset(data.w, VarContextAgg != "spatial")

#To test regressions between log-ratios. Change response variable and covariable
fit1 <- lmer(P50_yi ~ Ypd_yi +(1|VarContextAgg)+(1|PaperID), data = data) #change variables to obtain the different regressions
fit1 <- lmer(P50_yi ~ Ypd_yi +(1|VarContextAgg)+(1|PaperID), data = data.w) #change variables to obtain the different regressions
summary(fit1)
Anova(fit1)

r2 <- r.squaredGLMM(fit1)
r2


##### Creating Fig. 4 ####

data <- read_tsv('fileforcorrelations.tsv')
data.l<- subset(data, subset=(VarFactorAgg =="L"))
data.n<- subset(data, subset=(VarFactorAgg =="N"))
data.w<- subset(data, subset=(VarFactorAgg =="W"))
data.c<- subset(data, subset=(VarFactorAgg =="C"))
# Definition of colors
colors <- c("#9C9C9C", "#CCCC00", "#00C800", "#FF0200", "#008DF2")
names(colors) <- c("C", "L", "N", "TVPD", "W")
# Function to darken colors
darken_color <- function(color, amount = 0.3) {
  col <- col2rgb(color)
  col <- col * (1 - amount)
  col <- rgb(t(col), maxColorValue = 255)
  return(col)
}
dark_colors <- sapply(colors, darken_color, amount = 0.3)
names(dark_colors) <- names(colors)


# Filter data to remove missing values
data_filtered <- data %>%
  filter(!is.na(Ypd_yi) & !is.na(P50_yi))

# Fit the global mixed-effects model with lmerTest
fit_global <- lmer(P50_yi ~ Ypd_yi + (1|VarContextAgg) + (1|PaperID), data = data_filtered)

# Check if the model is singular
if (isSingular(fit_global)) {
  warning("The mixed model fit is singular. Consider simplifying the random effects structure.")
}

# Extract the p-value for the fixed effect
global_summary <- summary(fit_global)
global_p_value <- coef(global_summary)["Ypd_yi", "Pr(>|t|)"]  # Extract p-value

# Initialize the plot
plot <- ggplot(data_filtered, aes(x = Ypd_yi, y = P50_yi, fill = factor(VarFactorAgg))) +
  geom_point(shape = 21, size = 13, stroke = 1, color = "black") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_abline(slope = 1, linetype = "dashed", size = 2) +
  xlab(expression(paste('log-ratio ', Psi[PD]))) +
  ylab(expression(paste('log-ratio ', P[50]))) +
  theme_minimal()

# Add regression lines for significant subsets
for (group in unique(data_filtered$VarFactorAgg)) {
  subset_data <- filter(data_filtered, VarFactorAgg == group)
  if (nrow(subset_data) > 1) {
    fit_subset <- lm(P50_yi ~ Ypd_yi, data = subset_data)
    subset_summary <- summary(fit_subset)
    if (!is.null(subset_summary$coefficients[2, 4])) {
      subset_p_value <- subset_summary$coefficients[2, 4]
      if (subset_p_value < 0.05) {
        plot <- plot + geom_smooth(data = subset_data, method = "lm", size = 2, se = TRUE, 
                                   color = dark_colors[group])
      }
    }
  }
}

# Add the global regression line if significant
if (!is.null(global_p_value) && global_p_value < 0.05) {
  plot <- plot + geom_smooth(data = data_filtered, method = "lm", se = TRUE, linetype = "solid", 
                             size = 2, aes(group = 1), color = "black")
}

# Display the plot
print(plot)

reg <- plot +
  theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black", size = 2),
        axis.text.x = element_text(size = 54, color = "black", margin = margin(t = 5, r = 50, b = 20, l = 55)),
        axis.text.y = element_text(size = 54, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 0)),
        axis.title.x = element_text(size = 58, color = "black", face = "plain", vjust = 2, margin = margin(t = 10, r = 0, b = 2, l = 0)),
        axis.title.y = element_text(size = 58, color = "black", face = "plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x = element_line(color = "black", size = 1),
        axis.ticks.y = element_line(color = "black", size = 1),
        axis.ticks.length = unit(0.25, "cm"),
        legend.margin = margin(l = 0.6, unit = 'cm')) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(-1, 4), ylim = c(-1, 4)) + #For Ypd correlations
  scale_x_continuous(breaks = c(-1, 0, 1, 2 ,3, 4)) + #For Ypd correlations
  scale_y_continuous(breaks = c(-1, 0, 1, 2 ,3, 4)) #For Ypd correlations
  #coord_cartesian(xlim = c(-1, 1.6), ylim = c(-1, 1.6)) + #For correlations other than Ypd
  #scale_x_continuous(breaks = c(-1.5,-1, -0.5, 0, 0.5,1, 1.5)) + #For correlations other than Ypd
  #scale_y_continuous(breaks = c(-1.5,-1, -0.5, 0, 0.5,1, 1.5)) #For correlations other than Ypd
  
reg
# Save the final plot
#ggsave(filename = "Fig4-Ypd-Ptlp.png", plot = reg, dpi = 300, width = 10.5, height = 10.3) #For other trait-trait correlations

ggsave(filename = "Fig4-Ypd-P50.png", plot = reg, dpi = 300, width = 10.5, height = 11) #For Ypd-Other traits




###############################################################
#### Regressions between species means within environments #####
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

# Write the new dataframe to a file
write.table(final, 
            file = "regwithinenviroment.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)


##### Creating Fig. 5 ####
final <- read_tsv('regwithinenviroment.tsv')
reference <- subset(final, Environment == "R")
treatment <- subset(final, Environment == "T")
fit1 <- sma(log(abs(Ymd))~log(Ks+1), data=reference) #do it for each pair of traits
fit1 <- sma(log(abs(Ymd))~log(Ks+1), data=treatment) #do it for each pair of traits
summary(fit1)

fit_ref <- sma(log(Ks+1) ~ log(abs(Ymd)), data = reference)
summary_ref <- summary(fit_ref)
pval_ref <- fit_ref[["pval"]] 
pval_ref
fit_treat <- sma(log(Ks+1) ~ log(abs(Ymd)), data = treatment)
summary_treat <- summary(fit_treat)
pval_treat <- fit_treat[["pval"]]
pval_treat

#
pic <- ggplot(final, aes(x = log(abs(Ymd)), y = log(Ks+1), color = Environment)) +
  geom_point(size=10,shape=21,stroke = 3) + 
  labs(x = "Ymd", y = "Ks") +
  scale_color_manual(values = c("R" = "#0080FF", "T" = "#C00000")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
pic

# Reference
if (pval_ref < 0.05) {
  pic <- pic + geom_smooth(data = reference, aes(x = log(abs(Ymd)), y = log(Ks+1)),
                           method = "lm", se = TRUE, color = "#0080FF", fill = "#0080FF", 
                           linetype = "solid", size = 2.5, alpha = 0.1)
} else if (pval_ref < 0.10) {
  pic <- pic + geom_smooth(data = reference, aes(x = log(abs(Ymd)), y = log(Ks+1)),
                           method = "lm", se = TRUE, color = "#0080FF", fill = "#0080FF",
                           linetype = "dashed", size = 2.5, alpha = 0.1)
}

# Treatment
if (pval_treat < 0.05) {
  pic <- pic + geom_smooth(data = treatment, aes(x = log(abs(Ymd)), y = log(Ks+1)),
                           method = "lm", se = TRUE, color = "#C00000", fill = "#C00000",
                           linetype = "solid", size = 2.5, alpha = 0.3)
} else if (pval_treat < 0.10) {
  pic <- pic + geom_smooth(data = treatment, aes(x = log(abs(Ymd)), y = log(Ks+1)),
                           method = "lm", se = TRUE, color = "#C00000", fill = "#C00000",
                           linetype = "dashed", size = 2.5, alpha = 0.3)
}

print(pic)

plot=pic + theme(panel.grid = element_blank(),plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) + 
  theme(panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(), panel.background = element_rect(fill = "white", colour = "black", size=2),
        axis.text.x = element_text(size = 52, color = "black", margin = margin(t = 5, r = 50, b = 20, l = 55)), 
        axis.text.y = element_text(size = 52, color = "black", margin = margin(t = 0, r = 10, b = 5, l = 0)),
        axis.title.x = element_text(size = 58, color = "black", face="plain", vjust=2, margin = margin(t = 10, r = 0, b = 2, l = 0)),
        axis.title.y = element_text(size = 56, color = "black", face="plain", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.ticks.x=element_line(color = "black", size=1),
        axis.ticks.y=element_line(color = "black",  size=1),
        axis.ticks.length = unit(0.25, "cm"), 
        legend.margin=margin(l = 0.6, unit='cm')) + 
  theme(legend.position = "none")+
  coord_cartesian(xlim = c(-1,1.6), ylim = c(-0.2, 2.5)) + 
  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3,4,5,6)) + scale_y_continuous(breaks=c(-1,0,1,2,3))+
  xlab(expression(paste(log, " |",Psi[MD],"|"," ", (MPa))))+
  #ylab(expression(paste(log, " |",Pi[TLP],"|"," ", (MPa))))
  ylab(expression(paste(log," ", K[S], " (", kg, " ", m^-1, " ", s^-1, " ", MPa^-1, ")")))
plot

ggsave(filename = "Fig.5-logYmd-logKs.png", plot = plot, dpi = 300, width = 11, height = 11.5)
