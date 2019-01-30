################################################################################
#
#   File name: create_assoc.R
#
#   Authors: Stefano Pirro' ( s.pirro@qmul.ac.uk ) (*aka wynstep*)
#
#   Short Description:
#     Starting from a GWAS catalog association file, the script performs manipulations to make it ready to be used
#     with PRSice software
#
################################################################################

#===============================================================================
#    Load libraries
#===============================================================================
suppressMessages(library(optparse))
suppressMessages(library(biomaRt))
suppressMessages(library(reshape2))

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-g", "--gwas"), action="store", default=NA, type='character',
              help="File containing GWAS catalog file"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Default directory")
)

opt = parse_args(OptionParser(option_list=option_list))
gwasFile <- opt$gwas

#===============================================================================
#   Pre-processing of the input files
#===============================================================================

# Read gwas association data
gwasData <- read.table(gwasFile, sep = "\t", header = T, stringsAsFactors = F, row.names = NULL, quote="\"")

# Slimming gwasData by selecting just the columns we need for PRSice software
gwasSlimmedData <- gwasData[,c("SNPS","CHR_ID","CHR_POS","OR.or.BETA","P.VALUE")]

#===============================================================================
#   Mapping SNPs IDs with BiomaRt
#===============================================================================

## This process is necessary because the GWAS catalog file does not contain any information
## about the ref and alt allele. We retrieve the information by using biomaRt

# Loading Human Short Variants from Ensembl Variation 95
ensembl_snp <- useMart("ENSEMBL_MART_SNP")
ensembl_snp <- useDataset("hsapiens_snp", mart = ensembl_snp)

converted_ids <- getBM(attributes=c("refsnp_id", "chr_name", "allele"),
      filters = 'snp_filter',
      values = gwasSlimmedData$SNPS,
      mart = ensembl_snp)

# expand allele column by splitting ref (A1) and alt(s) (A2)
converted_ids<- cbind(converted_ids, colsplit(converted_ids$allele,"/",names=c("A1","A2")))
converted_ids$allele <- NULL
colnames(converted_ids) <- c("SNPS","CHR_ID","A1","A2")

#===============================================================================
#   Post-processing of the input files
#===============================================================================

gwasAssocData <- merge(converted_ids, gwasSlimmedData)
colnames(gwasAssocData) <- c("SNP","CHR","A1","A2","BP","OR","P")

# save table into a tab delimited file
outFile <- paste0(strsplit(gwasFile, "\\.(.+$)", perl=TRUE)[[1]],".prsice.assoc")
write.table(gwasAssocData, file = outFile, sep = "\t", quote = F, na = '0', col.names = T, row.names = F)
