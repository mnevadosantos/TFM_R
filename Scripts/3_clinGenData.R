#####################################################################################
# From the results of the two previous scripts, this one creates a dataframe with   #
# the clinical data together with the variants presented by each patient.           #
# Each row represents a patient. The presence or absence of each of the 159         #
# variants is indicated by a Boolean variable. The clinical variables will also     #
# appear in columns.                                                                #
#####################################################################################

d <- dir('~/Documentos/Olavide/TFM/Variantes_estudio')
setwd('~/Documentos/Olavide/TFM/Variantes_estudio')

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

setwd("~/Documentos/Olavide/TFM/Descargar")
fileFechas <- data.frame(read.table('fileDate2.tsv', header = F, sep = '\t'))
names(fileFechas) <- c('id2', 'filename', 'fecha', 'size', 'trash') # Si importo directamente del tsv que crea el fileDate.sh, añadir el name 'trash', y luego hacerla NULL
fileFechas$trash <- NULL
fileFechas <- fileFechas[grep('oxoG', fileFechas$filename),] # Nos quedamos sólo con los ficheros de interés
fileFechas <- fileFechas[order(fileFechas$filename, -abs(fileFechas$size)),] # ordenamos por nombre y por tamaño descendente.
fileFechas <- fileFechas[!duplicated(fileFechas$filename),] # Eliminamos los duplicados, y nos quedamos con los de fecha más antigua
# Ahora vamos a crear un campo con el identificador del paciente:
fileFechas$patient <- sapply(fileFechas$filename, function(i){
  unlist(strsplit(i, '_'))[1]
})
fileFechas <- fileFechas[order(fileFechas$patient, abs(fileFechas$fecha)),] # Ordenamos ahora por nombre de paciente y por fecha ascendente
fileFechas <- fileFechas[!duplicated(fileFechas$patient),] # Nos quedamos con la fecha más antigua
# Y ahora, convertimos los nombres de los ficheros vcf en csv, que es como los tengo con las variantes.
fileFechas$filename <- sub('.vcf$', '.csv', fileFechas$filename)
# Y con este dataframe ya puedo seguir trabajando, sin el problema de los pacientes duplicados.
# Contiene el nombre de los csv corresondientes a los vcf de mayor tamaño y con fecha más antigua (en el campo ##fileDate=)
############################################################################################################
setwd("~/Documentos/Olavide/TFM/Variantes_estudio")
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
####P PASO 4 ####
######################################################################################################
## PROCESAMOS AHORA LOS CSV CON DATOS CLÍNICOS Y DE SEGUIMIENTO, PARA CONSEGUIR UN SÓLO DATAFRAME QUE
## MEZCLAREMOS CON EL DE LOS DATOS DE VARIANTES.
## Los datos clínicos proceden de la base de datos TCGA_CLINICAL_PRAD_DATA.ods del estudio de Cabezuelo
#####################################################################################################

setwd('~/Documentos/Olavide/TFM/Clinica/')
datos_patient <- read.table('TCGA_CLINICAL_PRAD_DATA_PATIENT.csv', header = T, sep = ',')
datos_followup <- read.table('TCGA_CLINICAL_PRAD_DATA_FOLLOWUP.csv', header = T, sep = ',')
datos_drug <- read.table('TCGA_CLINICAL_PRAD_DATA_DRUG.csv', header = T, sep = ',')
datos_nte <- read.table('TCGA_CLINICAL_PRAD_DATA_NTE.csv', header = T, sep = ',')
datos_radiation <- read.table('TCGA_CLINICAL_PRAD_DATA_RADIATION.csv', header = T, sep = ',')
datos_v40 <- read.table('TCGA_CLINICAL_PRAD_DATA_V40.csv', header = T, sep = ',')

names(variantes_pacientes)[161] <- "bcr_patient_barcode" # Para poder hacer el merge

datos_global <- merge(x = variantes_pacientes, y = datos_patient, all.x = T)
datos_global <- merge(x = datos_global, y = datos_followup, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_drug, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_nte, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_radiation, all.x = T)
# datos_global <- merge(x = datos_global, y = datos_v40, all.x = T)
datos_global[datos_global == '[Not Available]'] <- NA
datos_global[datos_global == '[Unknown]'] <- NA
names(datos_global)[15:173] # Éstas son las variables correspondientes a las variantes.
write.csv2(datos_global, '~/Documentos/Olavide/TFM/datos_global.csv')