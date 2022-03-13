########################################################################################################
# This script extracts and sorts the variants of the genes of interest found in the vcf, which had been 
# annotated in the previous step, and had been converted into csv (one for each vcf file). The result 
# is a dataframe in which each row is a variant, and each column corresponds to its name, 
# characteristics and annotations.
########################################################################################################
# The following loop will iterate through the csv files, generating the variant name by joining the 
# chromosome number with the start position, and the letter of the variant nucleotide. In addition, at 
# each iteration it removes duplicates from this column. Finally, it adds the result to a dataframe 
# called f_all, which we will continue processing.
library(ggplot2)

d <- dir('~/Documentos/Olavide/TFM/Variantes_estudio')
setwd('~/Documentos/Olavide/TFM/Variantes_estudio')

f_ini_all <- data.frame()
f_all <- data.frame()
for (csv in d){
  f_ini <- read.table(csv, header = T, sep = ',', row.names = NULL)
  f_ini$X <- NULL
  f_ini$variante <- paste0(f_ini$start, sep='_', f_ini$Allele)
  f_ini$variante <- paste0(f_ini$seqnames, sep='_', f_ini$variante)
  f <- f_ini[!(duplicated(f_ini$variante)),]
  f_ini_all <- rbind(f_ini_all, f_ini)
  f_all <- rbind(f_all, f)
}

# Finally, f_all contains all the variants of each csv file. Now, calculating their  
# absolute and relative frequencies, storing them in the dataframe fvariantes:

fvariantes <- data.frame(t(table(f_all$variante)), row.names = NULL)
fvariantes <- fvariantes[,-1]
names(fvariantes) <- c('nom_variant', 'frec_abs')
fvariantes$frec_rel <- fvariantes$frec_abs/length(f_all$variante)
fvariantes_dep <- fvariantes[fvariantes$frec_rel>0.001,] # Variants with relative frecuencia greater than 1/1000

# Expanding the dataframe with annotations of each variant:

# First, creating a column with the transcript that corresponds to each variant,
# taking into account that the same variant can affect more than one transcript.
fvariantes$transcript <- sapply(fvariantes$nom_variant, function(i){
  paste0(unique(f_ini_all$SYMBOL[i == f_ini_all$variante]), collapse = ':|:')
})
# Now, doing something similar with the columns 'Consequence', 'IMPACT', 'BIOTYPE', and 'PolyPhen'
fvariantes$polyphen <- sapply(fvariantes$nom_variant, function(i){
  paste0(unique(f_ini_all$PolyPhen[i == f_ini_all$variante]), collapse = ':|:')
})

fvariantes$impact <- sapply(fvariantes$nom_variant, function(i){
  paste0(unique(f_ini_all$IMPACT[i == f_ini_all$variante]), collapse = ':|:')
})

fvariantes$biotype <- sapply(fvariantes$nom_variant, function(i){
  paste0(unique(f_ini_all$BIOTYPE[i == f_ini_all$variante]), collapse = ':|:')
})

# Substituting the string 'NA' for the logical value NA:

fvariantes[fvariantes == "NA"] <- NA

#### Some descriptive charts ####
fvariantes$nombre_clave <- paste0("V", '', 1:159)
ggplot(data = fvariantes, aes(x=reorder(nombre_clave, -frec_abs), y=frec_abs, fill=impact)) + 
  geom_bar(stat="identity", position = "dodge")
ggplot(data = fvariantes, aes(x=reorder(nombre_clave, -frec_abs), y=frec_abs, fill=impact)) + 
  geom_bar(stat="identity", position = "dodge") +
  facet_grid(impact ~ .)
