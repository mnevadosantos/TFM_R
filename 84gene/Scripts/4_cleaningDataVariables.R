#########################################################################################
## This script explores and cleans the variables of the dataset resulting for the       #
## previous script, 3_clinGenData.R.                                                    #
#########################################################################################

# First, loading the dataset:
datos_global <- read.csv("../Data/datos_global.csv", sep = ";", header = T, stringsAsFactors = T)
datos_global$X <- NULL
# The columnos corresponding to the variants are:
grep("^X+", names(datos_global)) # They are the columns 15 to 613

# Now exploring some variables:
nrow(datos_global[is.na(datos_global$vital_status),])
nrow(datos_global[datos_global$vital_status != "Alive",]) # Number of death patients
as.integer(datos_global$death_days_to[datos_global$vital_status != "Alive"]) # Survivance (days) of death patients
min(as.integer(datos_global$death_days_to[datos_global$vital_status != "Alive"]))
max(as.integer(datos_global$death_days_to[datos_global$vital_status != "Alive"]))

# Exploring causes of death:
summary(datos_global$cause_of_death)
datos_global$cause_of_death_source[datos_global$vital_status != "Alive"] # Uninteresting variable
datos_global$cause_of_death_source <- NULL # Deleting

# Exploring tumor status:
levels(datos_global$tumor_status)
summary(datos_global$tumor_status)

# Exploring adjuvancy with radiotherapy
summary(datos_global$radiation_treatment_adjuvant) # Too many NAs

# Exploring chemotherapic adjuvancy
summary(datos_global$pharmaceutical_tx_adjuvant) # Too many NAs

# Selecting clinical variables by NAs proportion (<= 10%:
clinical_var_index <- c(4:13, 614:684) # Clinical variables
sel_clinical_var <- c()
sel_clinical_index <- c()
for (index in clinical_var_index) {
  nas <- nrow(datos_global[is.na(datos_global[index]), ])
  if (nas < 50){
    sel_clinical_var <- c(sel_clinical_var, colnames(datos_global[index]))
    sel_clinical_index <- c(sel_clinical_index, index)
  }
}
# Studying the results:
length(sel_clinical_var)
sapply(sel_clinical_index, function(x){
  summary(datos_global[x])
})

# Selecting the columns initially considered useful for the study, and eliminating the rest:
sel_clinical_index <- c(1, sel_clinical_index, 14:612)
datos_global <- datos_global[sel_clinical_index]

# Dropping useless variables. First,defining which variables will be dropped:
drop_variables <- c("last_contact_days_to", "death_days_to", "gender", "prospective_collection", "retrospective_collection", 
                    "history_neoadjuvant_treatment", "initial_pathologic_dx_year", "lymph_nodes_examined", "days_to_initial_pathologic_diagnosis", 
                    "extranodal_involvement", "icd_10", "icd_o_3_histology", "icd_o_3_site", "informed_consent_verified", "pathologic_M", 
                    "pathologic_stage", "system_version", "tumor_tissue_site")
datos_global[colnames(datos_global) %in% drop_variables] <- NULL


# Creating a new clinical variable: the WHO's grade group:
datos_global$WHO_group[datos_global$gleason_score <= 6] <- 1
datos_global$WHO_group[datos_global$gleason_score == 8] <- 4
datos_global$WHO_group[datos_global$gleason_score > 8] <- 5
datos_global$WHO_group[datos_global$gleason_pattern_primary == 3 & 
                   datos_global$gleason_pattern_secondary == 4] <- 2
datos_global$WHO_group[datos_global$gleason_pattern_primary == 4 & 
                   datos_global$gleason_pattern_secondary == 3] <- 3

# Writing data set to a file:
write.csv2(datos_global, "datos_global_clean.csv")

