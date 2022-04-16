################################################################################
### This script processes the vcf files to select the variants that affect #####
### the study genes.                                                       #####
################################################################################
library(VariantAnnotation)
library(vcfR) #Librería para manejar ficheros vcf como objetos vcf
library(ensemblVEP)

wdini <- '~/Documentos/Olavide/TFM' # Initial working directory
setwd(wdini)
# Loading the downloaded TCGA manifest
manifiesto <- read.table("gdc_manifest.2021-10-15.csv", header = T, sep = "\t")
# Selecting, in the manifest, the names of file with the "oxoG" extension
manifiesto2 <- manifiesto[grep("oxoG", manifiesto$filename),]
# Selecting the largest file, among two with the same name
# First, ordering by file name and, in descending order, by size:
manifiesto3 <- manifiesto2[order(manifiesto2$filename, -abs(manifiesto2$size)),]
# Files with the same name are now contiguos. The largest is always first.
# Deleting duplicated, the function eliminates always the second line:
manifiesto3 <- manifiesto3[!duplicated(manifiesto3$filename),]
# Now, the loop changes directories, reads each vcf, parses it, selects variants
# with quality greater than 99, selects germ variants, and writes this information
# to a new vcf in compressed format, with the extension germ.vcf.gz
for(carpeta in manifiesto3$id){
  wd <- paste0(wdini, '/Descargar/', carpeta)
  setwd(wd)
  archiv <- dir(pattern = '.vcf$') # the input file
  archiv2 <- paste0(sub('.vcf$','', archiv), '.germ.vcf.gz') # Output file
  # Parsing vcf:
  vcf2 <- read.vcfR(archiv)
  vcf3 <- vcf2[as.numeric(vcf2@fix[, 'QUAL']) > 99,] # Selecting variants by quality
  vcf4 <- vcf3[grepl('Germline', vcf3@fix[, 'INFO']),] # Selecting germ variants.
  write.vcf(vcf4, archiv2) # Writing a new compressed vcf.
}
# Selecting the genes of the study in the following list:

lista_genes <- list('CYP3A4', 'GSTM1', 'CYP17A1', 'UGT2B17', 'MTHFR', 'CYP1A1', 'NAT1', 
                    'NAT2', 'GSTT1', 'OGG1', 'TXNRD1', 'PRDX3', 'GSTP1', 'GPX1', 'SOD2', 
                    'GSTA1', 'GSTA2', 'GSTA3', 'GSTA4', 'GSTA5', 'GSTT2', 'GSTM5', 'GSTM2', 
                    'GSTM3', 'GSTM4', 'GSTK1', 'GSTO2', 'GSTO1', 'SOD1', 'SOD3', 'DLX1', 
                    'TDRD1', 'OR51C1P', 'OR51E2', 'HOXC6', 'APOF', 'SERPINA11', 'AMACR', 
                    'KRT13', 'ONECUT2', 'C1QTNF3-AMACR', 'SLIT1', 'SERPINA5', 'ACTC1', 'TGM3', 
                    'SCGB1A1', 'ACSM1', 'CLCA2', 'HPN', 'GPX2', 'ATP8A2', 'KCNJ15', 'C2orf72', 
                    'MSLN', 'CYP4F22', 'CRTAC1', 'HSD17B13', 'SBSPON', 'DUOXA2', 'KRT16', 'MATK', 
                    'CA14', 'QPRT', 'MFSD2A', 'C2orf88', 'CLU', 'DUOX2', 'STAC', 'AOX1', 'CXCR2', 
                    'TFF3', 'WFDC2', 'SRARP', 'ERG', 'PTGS1', 'TMEM132C', 'COL2A1', 'GCNT4', 
                    'GNAO1', 'KLHL14', 'SNAP25', 'NDRG4', 'FAM167A', 'CLCA4')
# The following loop annotates the variants with the VEP script and selects those'
# annotated variants that correspond to any of the genes in the previous list.
for (carpeta in manifiesto3$id) {
  wd <- paste0(wdini, '/Descargar2/', carpeta) # Reads the name of each directory in the manifest file
  setwd(wd) 
  archiv <- dir(pattern = '.vcf.gz$') # Name of the input file
  archiv2 <- paste0(sub('.germ.vcf.gz','', archiv), '.csv') # Name of the output file
  # The following configures the parameters of the VEP script:
  param <- VEPFlags(flags = list(regulatory=TRUE, offline=TRUE, cache=TRUE,
                                 assembly='GRCh38', polyphen = 'b', sift = 'b',
                                 symbol = TRUE, pubmed = TRUE))
  print(paste("Cargado fichero: ", archiv)) # This is just to note the progress of the loop
  gr <- ensemblVEP(archiv, param) # This command annotates the vcf using the vep script
  gr2 <- subset(gr, SYMBOL %in% lista_genes) # This selects the variants of the genes from the list
  gr_df <- as.data.frame(gr2, 'data.frame', row.names = NULL) # Convert to a dataframe
  write.csv(gr_df, file = archiv2) # Writes the output file
  print(paste("Procesado con exito: ", archiv2)) # Again to note progress
}
# So far, for each vcf file we have obtained a csv with the annotated variants of the genes of interest.

