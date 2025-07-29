#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse args, with fallback defaults
wd <- if (length(args) >= 1) args[1] else "/data/adalla/final/"
si <- if (length(args) >= 2) args[2] else 80    
coverage <- if (length(args) >= 3) args[3] else 40


# Load libraries
library(tidyverse)
library(survival)
library(survminer)
library(readxl)

setwd(wd)

# Load the binary matrix (species x samples)
filename = paste0("presence_absence_matrix_", si, "_", coverage, ".csv")
species_matrix <- read.csv(filename, row.names = 1)

# Load metadata
metadata <- read_excel("DLC380_RNAseq_seq_source_withBatches.xlsx")
# metadata[metadata$Library_ID == "A38160", ]

# Now if we checked that every samplename of metadata is present in species_matrix with 
# this line of code, we could obtain FALSE. That0s because in the filtered_infections_file
# only positive infections are stored. So the samples without infection would be missing in the 
# presence_absence_matrix and therefore in the species_matrix object.
# So in this step we have to re-add them
#all(colnames(species_matrix) %in% metadata[[1]])

sample_ids_matrix <- colnames(species_matrix)
sample_ids_metadata <- metadata$Library_ID

# Identify missing samples
missing_in_samples <- setdiff(sample_ids_metadata, sample_ids_matrix)
#missing_in_samples

# Add them as zero columns
for (s in missing_in_samples) {
  species_matrix[[s]] <- 0
}

# There are 2 samples missing in metadata: "A38160" "A38176"
missing_in_metadata <- setdiff(sample_ids_matrix, sample_ids_metadata)
#missing_in_metadata
#metadata[metadata$Library_ID == "A38160", ]
#metadata[metadata$Library_ID == "A38176", ]

# transpose the input binary matrix of filtered infection predictions 
species_matrix_t <- t(species_matrix)

# Remove the 2 NA rows from metadata 
metadata <- metadata[!is.na(metadata$Library_ID), ]

# Assign to rownames of metadata the samples names 
rownames(metadata) <- metadata$Library_ID

# Subset both metadata and species_matrix_t to shared samples
shared_ids <- intersect(rownames(species_matrix_t), rownames(metadata))
species_matrix_t <- species_matrix_t[shared_ids, ]
metadata <- metadata[shared_ids, ]

# Merge metadata with species data, so that in a single row we have all the informations for that 
# samples and all the presence/absence value for each species
data <- cbind(metadata, species_matrix_t)
#head(data)

# factor of the types of DLBCL present in the dataset
lymphoma_type <- factor(metadata$COO)


##########################################################################################################

#### ARE THERE IS MICROBIAL PRESENCE ASSOCIATED WITH ANY TYPE OF LYMPHOMA?

# Filter out 'Unclassified' or NA lymphoma types
filtered_indices <- lymphoma_type %in% c("ABC", "GCB")
species_matrix_t_filtered <- species_matrix_t[filtered_indices, ]
lymphoma_type_filtered <- droplevels(lymphoma_type[filtered_indices])

# Build contingency table
any_microbe_present <- factor(as.integer(rowSums(species_matrix_t_filtered) > 0), 
                             levels = c(0, 1))

# Create contingency table
tbl <- table(any_microbe_present, lymphoma_type_filtered)
rownames(tbl) <- c("Absent", "Present")  # Label rows clearly
#attr(tbl, "dimnames")[[1]] <- NULL       # Remove row dimension name

# Fisher's exact test
any_infection <- fisher.test(tbl)

# Output results
cat("P-value for any microbe vs lymphoma type:", any_infection$p.value, "\n")

# Format p-value string
pval_str <- ifelse(any_infection$p.value < 0.001, 
                   format(any_infection$p.value, scientific = TRUE, digits = 2),
                   sprintf("%.3f", any_infection$p.value))
param_str <- paste0("s.i.: ", si, "%\tcoverage: ", coverage, "%")

# Save plot
filename <- "fisher_any_infection.png"
file_path <- paste0(wd, "plots/", "si", si, "-cov", coverage, "/", filename)
print(paste("Fisher's exact test saved at", file_path))

# Open PNG device
png(filename = file_path, width = 800, height = 600)

print(t(tbl))

# Plot mosaic
mosaicplot(t(tbl),
           main = "",
           xlab = "Lymphoma Subtype",
           ylab = "Microbe Presence",
           color = c("skyblue", "salmon"),
           cex.axis = 1.25
           )

title(main = paste0(param_str, "\nP-value: ", pval_str), cex.main = 1.5)

library(survival)
library(survminer)

# Prepare event indicator (assuming metadata$OS.Status available)
metadata$OS.Event <- ifelse(metadata$OS.Status == "Dead", 1, 0)

# Filter lymphoma types to ABC/GCB only as before
filtered_indices <- lymphoma_type %in% c("ABC", "GCB")
species_matrix_t_filtered <- species_matrix_t[filtered_indices, ]
lymphoma_type_filtered <- droplevels(lymphoma_type[filtered_indices])

# Create the any_microbe_present variable
any_microbe_present <- as.integer(rowSums(species_matrix_t_filtered) > 0)

# Subset metadata for filtered samples
metadata_filtered <- metadata[filtered_indices, ]

# Prepare data for survival
valid <- !is.na(metadata_filtered$OS.Time) & !is.na(metadata_filtered$OS.Event) & !is.na(any_microbe_present)

df_surv <- data.frame(
  time      = metadata_filtered$OS.Time[valid],
  event     = metadata_filtered$OS.Event[valid],
  infection = factor(any_microbe_present[valid], levels = c(0,1), labels = c("Absent", "Present"))
)

# Survival object
surv_obj <- Surv(df_surv$time, df_surv$event)

# Log-rank test
lr <- survdiff(surv_obj ~ infection, data = df_surv)
p_km <- 1 - pchisq(lr$chisq, df = 1)

# Cox proportional hazards model (no covariates here, add if you want)
cox_fit <- coxph(surv_obj ~ infection, data = df_surv)
summary_cox <- summary(cox_fit)

# Extract hazard ratio and confidence intervals
HR <- summary_cox$conf.int["infectionPresent", "exp(coef)"]
HR_95L <- summary_cox$conf.int["infectionPresent", "lower .95"]
HR_95U <- summary_cox$conf.int["infectionPresent", "upper .95"]
p_cox <- summary_cox$coefficients["infectionPresent", "Pr(>|z|)"]

# Print results
cat("Log-rank p-value:", p_km, "\n")
cat(sprintf("Cox HR=%.2f (95%% CI: %.2f - %.2f), p=%.3g\n", HR, HR_95L, HR_95U, p_cox))

# Format p-value string
pval_str <- ifelse(p_cox < 0.001, format(p_cox, scientific = TRUE, digits = 2), sprintf("%.3f", p_cox))

# Plot Kaplan-Meier curves
fit <- survfit(surv_obj ~ infection, data = df_surv)

filename <- paste0("survival_analysis_any_microbe.png")
file_path <- paste0(wd, "plots/", "si", si, "-cov", coverage, "/", filename)
print(paste("Survival analysis saved at", file_path))

png(filename = file_path, width = 800, height = 600)

title_txt <- sprintf("Any Microbe Presence\nHR=%.2f (95%% CI: %.2f–%.2f), p=%s (Cox regression)", HR, HR_95L, HR_95U, pval_str)

g <- ggsurvplot(fit,
                data = df_surv,
                risk.table = TRUE,
                conf.int = TRUE,
                title = title_txt,
                legend.title = "Microbe Presence",
                legend.labs = c("Absent", "Present"),
                palette = c("blue", "red"),
                pval = TRUE,
                xlab = "Times (years)")

print(g)
dev.off()

##########################################################################################################

#### ARE THERE SPECIFIC MICROBES ASSOCIATED WITH A SPECIFIC TYPE OF LYMPHOMA???

### CHI2 TEST

# Run a chi2 test on each species. 
# The chi2 test determine if there is a significant association between categorical variables. 
# It compares the observed frequencies of occurrences in different categories with 
# the frequencies expected if there were no associations between the variables.
# To run a chi2 test we need a contingency table, which we can obtain in R with the function table()
# Here we use apply to run a chi2 test on each species independently.
# An example of contingency table is 
#       lymphoma_type
#     x   ABC GCB Unclass
#     0  97 173      35
#     1   0   2       0
#
# We detected this species in 2 sample, both of them of lymphoma type GCB.

# Filter to species present in at least N samples (e.g., 5)
#present_counts <- colSums(species_matrix_t > 0)
#species_matrix_t_filtered <- species_matrix_t[, present_counts >= 5]
#
## Then redo the test
#results <- apply(species_matrix_t_filtered, 2, function(x) {
#  test <- chisq.test(table(x, lymphoma_type))
#  return(test$p.value)
#})
#
## sig_results <- data.frame( species = names(p_adj), 
##                            p_value = results, 
##                            p_adj = p_adj )
#
#sig_results <- data.frame( species = names(results), 
#                           p_value = results)
#
#
#sig_results <- sig_results[sig_results$p_value < 0.05 & !is.na(sig_results$p_value), ]
#sig_results
#
# Sign of the relationship:
#for (species in sig_results$species) {
#  cat("\n=== ", species, " ===\n")
#  tbl <- table(species_matrix_t[, species], lymphoma_type)
#  print(tbl)
#}


##### FISHER'S EXACT TEST

# Since we're doing fihser's exact test, it is better to have a 2x2 contingency table
# Filter out 'Unclassified' samples
filtered_indices <- lymphoma_type %in% c("ABC", "GCB")

species_matrix_t_filtered <- species_matrix_t[filtered_indices, ]
lymphoma_type_filtered <- droplevels(lymphoma_type[filtered_indices])

# Now run Fisher's test only on ABC vs GCB
results_fisher <- sapply(colnames(species_matrix_t_filtered), function(species) {
  tbl <- table(species_matrix_t_filtered[, species], lymphoma_type_filtered)
  
  # Should be 2x2 now
  if (all(dim(tbl) == c(2, 2))) {
    test <- fisher.test(tbl)
    return(test$p.value)
  } else {
    return(NA)
  }
})

# Adjust p-values for multiple testing
p_adj_fisher <- p.adjust(results_fisher, method = "fdr")

# Compile results in a data frame
sig_results_fisher <- data.frame(
  species = names(p_adj_fisher),
  p_value = results_fisher,
  p_adj = p_adj_fisher
)

# Keep only significant species
sig_results_fisher <- sig_results_fisher[!is.na(sig_results_fisher$p_value) & sig_results_fisher$p_value < 0.1, ]

# Print significant species with their contingency tables (using filtered lymphoma_type!)
#for (species in sig_results_fisher$species) {
#  cat("\n=== ", species, " ===\n")
#  print(table(species_matrix_t_filtered[, species], lymphoma_type_filtered))
#}


library(ggplot2)

for (species in sig_results_fisher$species) {
  # convert it into a factor so it is guaranteed that when i assign rownames, "Present" will
  # always be assigned to 1 and "Absence" to 0.
  binary_species <- factor(species_matrix_t_filtered[, species], levels = c(0, 1))
  tbl <- table(binary_species, lymphoma_type_filtered)
  rownames(tbl) <- c("Absent", "Present")
  
  print(tbl)
  
  # Get raw p-value
  pval <- sig_results_fisher$p_value[sig_results_fisher$species == species]
  
  # Format p-value nicely (scientific if small)
  pval_str <- ifelse(pval < 0.001, format(pval, scientific = TRUE, digits = 2),
                     sprintf("%.3f", pval))
  
  filename <- paste0("fisher_", species, ".png")
  file_path <- paste0(wd, "plots/", "si", si, "-cov", coverage, "/", filename)
  print(paste("Fisher's exact test saved at ", file_path))
  
  # Open PNG device
  png(filename = file_path, width = 800, height = 600)
  
  # Create plot with p-value in title
  mosaicplot(t(tbl),
             main = "",
             xlab = "",
             color = c("skyblue", "salmon"),
             cex.axis = 1.2)
  
  title(xlab = "Lymphoma Subtype", cex.lab = 1.4)
  title(ylab = "Species presence", cex.lab = 1.4)
  title(main = paste0(species, "\nP-value: ", pval_str),
        cex.main = 1.5)

  dev.off()
}

print(" ")

#######################################################################################################################

# ARE THERE SPECIES ASSOCIATED WITH CHANGES IN THE SURVIVAL RATE?

#  To test the relationship between a cathegorical and a continuous variable i can choose to do the ttest
#  (parametrical) or the ranksum test (non-parametrical). I can do the ttest only if the 
#  distribution of the continuous variable is normal

# graphical representation
# hist(metadata$OS.Time, main="Histogram", xlab="Value")
# qqnorm(metadata$OS.Time)


# statistical test: H0 = the distribution of the continous variable is normal
#     if p-value < 0.05 -> reject the null hypothesis -> distribution is not normal
# shapiro.test(metadata$OS.Time)
    # p-value = 1.222e-09


# statistical test: H0 = the distribution of the continous variable is normal
#     if p-value < 0.05 -> reject the null hypothesis -> distribution is not normal
# library(nortest)
# sf.test(metadata$OS.Time)
    # p-value = 2.724e-08


# So our variable is not normally distributed -> non-parametric test wilcoxon ranksum test

### prova
#species = "Malassezia_restricta"
#all(rownames(species_matrix_t) == rownames(metadata))  # Should be TRUE
#
## If FALSE, reorder metadata to species_matrix_t sample order:
#metadata <- metadata[rownames(species_matrix_t), ]
#
## Extract infection status for Malassezia_restricta:
#infection_status <- species_matrix_t[, "Malassezia_restricta"]
#
# Check if infection_status has two groups (0 and 1):
#unique(infection_status)
#
## Extract survival times:
#survival_time <- metadata$OS.Time
#
## Check if both groups have enough data:
#table(infection_status)
## Make sure both groups have > 0 samples
#
## Run Wilcoxon rank-sum test:
#wilcox_result <- wilcox.test(survival_time ~ infection_status, conf.int = TRUE)
#
# See results:
#print(wilcox_result)
#
# p-value:
#wilcox_result$p.value

### versione finale
#results_wilcox <- sapply(colnames(species_matrix_t), function(species) {
#  infection_status <- species_matrix_t[, species]
#  
#  # Skip species with only one infection status group
#  if (length(unique(infection_status)) < 2) return(NA)
#  
#  # Extract survival times
#  survival_times <- metadata$OS.Time
#  
#  # Remove samples with NA survival time or infection status
#  valid_idx <- !is.na(survival_times) & !is.na(infection_status)
#  
#  if (sum(valid_idx & infection_status == 0) == 0 || sum(valid_idx & infection_status == 1) == 0) {
#    # If one group is empty after removing NAs, skip
#    return(NA)
#  }
#  
#  # Run Wilcoxon rank-sum test (Mann-Whitney U test)
#  test <- wilcox.test(survival_times ~ infection_status, conf.int = TRUE)
#  
#  return(test$p.value)
#})
#
#p_adj_wilcox <- p.adjust(results_wilcox, method = "fdr")
#
#sig_results_wilcox <- data.frame(
#  species = names(p_adj_wilcox),
#  p_value = results_wilcox,
#  p_adj = p_adj_wilcox
#)
#
## Keep only significant species
#sig_results_wilcox <- sig_results_wilcox[!is.na(sig_results_wilcox$p_value) & sig_results_wilcox$p_value < 0.05, ]
#
## Print significant species with their contingency tables
#for (species in sig_results_wilcox$species) {
#  cat("\n=== ", species, " ===\n")
#  print(table(species_matrix_t[, species], lymphoma_type))
#}


#############################################################################################################################

# KEPLAN-MEIER SURVIVAL PLOTS

library(survival)
library(survminer)

# Prepare event indicator
metadata$OS.Event <- ifelse(metadata$OS.Status == "Dead", 1, 0)

# 0) Filter the infection table so that only pretty common species are taken into consideration:
min_pres <- 16  # i.e. 5% of the whole dataset
common_species <- names(which(colSums(species_matrix_t == 1, na.rm=TRUE) >= min_pres))


# 1) Gather KM p‐values AND Cox model HR + p‐values
results <- lapply(common_species, function(species) {
  inf   <- species_matrix_t[, species]
  valid <- !is.na(metadata$OS.Time) & !is.na(metadata$OS.Event) & !is.na(inf)
  if (length(unique(inf[valid])) < 2) return(NULL)
  
  df_sub <- data.frame(
    time   = metadata$OS.Time[valid],
    event  = metadata$OS.Event[valid],
    infection = factor(inf[valid], levels=c(0,1), labels=c("Absent","Present"))
  )
 
  surv_obj <- Surv(df_sub$time, df_sub$event)
  
  
  # Log‐rank test to determine if the KM curves stratified by the infection are significantly different
  lr   <- survdiff(surv_obj ~ infection, data=df_sub)
  p_km <- 1 - pchisq(lr$chisq, df=1)
  
  # Cox PH
  cfit <- coxph(surv_obj ~ infection, data=df_sub)
  s <- summary(cfit)
  
  data.frame(
    species = species,
    p_km    = p_km,
    HR      = s$conf.int["infectionPresent", "exp(coef)"],
    HR_95L  = s$conf.int["infectionPresent", "lower .95"],
    HR_95U  = s$conf.int["infectionPresent","upper .95"],
    p_cox   = s$coefficients["infectionPresent","Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})


# bind and adjust
df <- do.call(rbind, results)
df$p_adj_km  <- p.adjust(df$p_km, method = "fdr")
df$p_adj_cox <- p.adjust(df$p_cox, method = "fdr")

# significant by KM or Cox
sig_km  <- subset(df, p_adj_km  < 0.1)
sig_cox <- subset(df, p_cox < 0.1)

# print(sig_km)
# print(sig_cox)

# 2) Plot KM + annotate with HR
for (sp in sig_cox$species) {
  inf   <- species_matrix_t[, sp]
  valid <- !is.na(metadata$OS.Time) & !is.na(metadata$OS.Event) & !is.na(inf)
  df_sub <- data.frame(
    time   = metadata$OS.Time[valid],
    event  = metadata$OS.Event[valid],
    infection = factor(inf[valid], levels = c(0,1), labels = c("Absent","Present"))
  )
  
  # KM curve
  fit <- survfit(Surv(time, event) ~ infection, data = df_sub)
  
  # save the plot in the output plot directory
  filename <- paste0("survival_analysis_", sp, ".png")
  file_path <- paste0(wd, "plots/", "si", si, "-cov", coverage, "/", filename)
  print(paste("Survival analysis saved at ", file_path))
  
  # Open PNG device
  png(filename = file_path, width = 800, height = 600)
  
  # pull HR/p for title
  row <- df[df$species == sp,]
  title_txt <- sprintf("%s\nHR=%.2f (95%% CI: %.2f–%.2f), p=%.3f (Cox regression)",
                       sp, row$HR, row$HR_95L, row$HR_95U, row$p_cox)
  
  # Plot
  g <- ggsurvplot(fit,
             data        = df_sub,
             risk.table  = TRUE,
             conf.int    = TRUE,
             title       = title_txt,
             legend.title= sp,
             legend.labs = c("Absent","Present"),
             palette     = c("blue","red"),
             pval        = TRUE,
             xlab = "Times (years)")
  print(g)
  
  dev.off()
}


#############################################################################################################################

# KEPLAN-MEIER SURVIVAL PLOTS

library(survival)
library(survminer)
filename = paste0("presence_absence_matrix_", si, "_", coverage, ".csv")
species_matrix <- read.csv(filename, row.names = 1)

# Load metadata
metadata <- read_excel("DLC380_RNAseq_seq_source_withBatches.xlsx")
# metadata[metadata$Library_ID == "A38160", ]

# Now if we checked that every samplename of metadata is present in species_matrix with 
# this line of code, we could obtain FALSE. That0s because in the filtered_infections_file
# only positive infections are stored. So the samples without infection would be missing in the 
# presence_absence_matrix and therefore in the species_matrix object.
# So in this step we have to re-add them
#all(colnames(species_matrix) %in% metadata[[1]])

sample_ids_matrix <- colnames(species_matrix)
sample_ids_metadata <- metadata$Library_ID

# Identify missing samples
missing_in_samples <- setdiff(sample_ids_metadata, sample_ids_matrix)
#missing_in_samples

# Add them as zero columns
for (s in missing_in_samples) {
  species_matrix[[s]] <- 0
}

# There are 2 samples missing in metadata: "A38160" "A38176"
missing_in_metadata <- setdiff(sample_ids_matrix, sample_ids_metadata)
#missing_in_metadata
#metadata[metadata$Library_ID == "A38160", ]
#metadata[metadata$Library_ID == "A38176", ]

# transpose the input binary matrix of filtered infection predictions 
species_matrix_t <- t(species_matrix)

# Remove the 2 NA rows from metadata 
metadata <- metadata[!is.na(metadata$Library_ID), ]

# Assign to rownames of metadata the samples names 
rownames(metadata) <- metadata$Library_ID

# Subset both metadata and species_matrix_t to shared samples
shared_ids <- intersect(rownames(species_matrix_t), rownames(metadata))
species_matrix_t <- species_matrix_t[shared_ids, ]
metadata <- metadata[shared_ids, ]

# Merge metadata with species data, so that in a single row we have all the informations for that 
# samples and all the presence/absence value for each species
data <- cbind(metadata, species_matrix_t)
#head(data)

# factor of the types of DLBCL present in the dataset
lymphoma_type <- factor(metadata$COO)


# Recalculate OS.Event to ensure it's numeric
metadata$OS.Event <- ifelse(metadata$OS.Status == "Dead", 1, 0)
metadata$OS.Event <- as.numeric(metadata$OS.Event)

# Merge metadata and infection data again to ensure all data is updated
data <- cbind(metadata, species_matrix_t)

# Remove rows with missing time or event values
data_clean <- data[!is.na(data$OS.Time) & !is.na(data$OS.Event), ]

# Confirm OS.Event is numeric
str(data_clean$OS.Event)

# Create the survival object
surv_obj <- Surv(time = data_clean$OS.Time, event = data_clean$OS.Event)

# Ensure infection column is a factor
#data_clean$Escherichia_coli <- factor(data_clean$Escherichia_coli, levels = c(0, 1), labels = c("Absent", "Present"))

cox_model <- coxph(surv_obj ~ COO + Escherichia_coli, data = data_clean)

summary(cox_model)

fit_km <- survfit(surv_obj ~ COO + Escherichia_coli, data = data_clean)

ggsurvplot(fit_km,
           data = data_clean,
           pval = TRUE,
           #conf.int = TRUE,
           #risk.table = TRUE,
           legend.title = "E. coli infection",
           xlab = "Time (years)")

