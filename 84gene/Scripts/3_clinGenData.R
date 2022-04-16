#####################################################################################
# From the results of the two previous scripts, this one creates a dataframe with   #
# the clinical data together with the variants presented by each patient.           #
# Each row represents a patient. The presence or absence of each of the 159         #
# variants is indicated by a Boolean variable. The clinical variables will also     #
# appear in columns.                                                                #
#####################################################################################

d <- dir('../Data/Variantes_estudio/')
setwd('../Data/Variantes_estudio/')

# Generating a vector with each patient's name:
nombre <- c()
for (i in 1:length(d)){
  nombre <- c(nombre, unlist(strsplit(d[i], '_'))[1])
}
variantes_pacientes <- data.frame(unique(nombre))

#####################################################################################
# It is evident that there are 54 patients with more than one vcf. To solve this,   #
# it is necessary to use the bash script called fileDate.sh, which gets the date    #
# from each csv file and stores it in fileDate.tsv.                                 #
#####################################################################################

setwd("../")
fileFechas <- data.frame(read.table('fileDate.tsv', header = F, sep = '\t'))
# After importing fileDate.tsv, naming columns, and deleting the fifth one:
names(fileFechas) <- c('id2', 'filename', 'fecha', 'size', 'trash')
fileFechas$trash <- NULL
fileFechas <- fileFechas[grep('oxoG', fileFechas$filename),] # In the imported dataframe, keep only the files of interest:
fileFechas <- fileFechas[order(fileFechas$filename, -abs(fileFechas$size)),] # Sorting by filename and descending size.
fileFechas <- fileFechas[!duplicated(fileFechas$filename),] # When duplicates are deleted, only those with the oldest date are kept.
# Creating a column with patient ids:
fileFechas$patient <- sapply(fileFechas$filename, function(i){
  unlist(strsplit(i, '_'))[1]
})
fileFechas <- fileFechas[order(fileFechas$patient, abs(fileFechas$fecha)),] # Now sorting by patient id and ascending date.
fileFechas <- fileFechas[!duplicated(fileFechas$patient),] # Keeping the earliest date.
# Changing the vcf extensions to csv.
fileFechas$filename <- sub('.vcf$', '.csv', fileFechas$filename)

# Problem of duplicate patients has been solved.
# The dataframe contains the name of the csv corresponding to the original vcf files 
# of larger size and older date in the ##fileDate field.
############################################################################################################
setwd("../Data/Variantes_estudio/")
variantes_pacientes <- data.frame()
for (i in fileFechas$filename){
  fich <- read.table(i, header = T, sep = ',')
  fich$nom_variant <- paste0(fich$seqnames, sep ='_', fich$start)
  fich$nom_variant <- paste0(fich$nom_variant, '_', fich$Allele)
  for (j in fvariantes$nom_variant){
    variantes_pacientes[i,j] <- j %in% fich$nom_variant
  }
}
variantes_pacientes$Muestra <- row.names(variantes_pacientes)
variantes_pacientes$Paciente <- sapply(variantes_pacientes$Muestra, function(i){
  unlist(strsplit(i, '_'))[1]
})
row.names(variantes_pacientes) <- NULL

######################################################################################################
## Processing of clinical data to obtain a single data frame to be merged with that of the variants. #
## Data from the TCGA_CLINICAL_PRAD_DATA.ods file of the Cabezuelo study.                            #
#####################################################################################################

setwd('../')
# First, reading all the clinical data:
datos_patient <- read.table('TCGA_CLINICAL_PRAD_DATA_PATIENT.csv', header = T, sep = ',')
datos_followup <- read.table('TCGA_CLINICAL_PRAD_DATA_FOLLOWUP.csv', header = T, sep = ',')
datos_drug <- read.table('TCGA_CLINICAL_PRAD_DATA_DRUG.csv', header = T, sep = ',')
datos_nte <- read.table('TCGA_CLINICAL_PRAD_DATA_NTE.csv', header = T, sep = ',')
datos_radiation <- read.table('TCGA_CLINICAL_PRAD_DATA_RADIATION.csv', header = T, sep = ',')
datos_v40 <- read.table('TCGA_CLINICAL_PRAD_DATA_V40.csv', header = T, sep = ',')

names(variantes_pacientes)[601] <- "bcr_patient_barcode" # Changing this name in order to make the "merge".

datos_global <- merge(x = variantes_pacientes, y = datos_patient, all.x = T)
datos_global <- merge(x = datos_global, y = datos_followup, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_drug, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_nte, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_radiation, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_v40, all.x = T)
# Now, correcting some things:
datos_global[datos_global == '[Not Available]'] <- NA
datos_global[datos_global == '[Unknown]'] <- NA
names(datos_global)[15:613] # These columns correspond to the variants.
write.csv2(datos_global, 'datos_global.csv')
