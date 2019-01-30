#/usr/bin/env python
# Coder: -.- Dr. Stefano Pirro (aka wynstep)

# This script (PRSCalc) aims at performing a brief calculation of one PRS,
# starting from a VCF file calculated from a single patient.
# There are two different ways for calculating:
#   prsice: takes advantage of the PRSice software for calculating the Polygenic Risk Score
#   manual: use a custom calculation method for extracting a single PRS for the patient

# import libraries
import os, sys, optparse, datetime
sys.path.append('scripts/')

# importing fixed vars and methods
from vars import * # import file with vars
from functions import *

# Loading arguments necessary for the analysis
parser = optparse.OptionParser(description='Welcome to PRSCalc')
parser.add_option('-f','--vcf', dest="vcf", help='vcf file to analyse')
parser.add_option('-g','--gwas', dest="gwas", help='gwas association file to compare (PRSice format)')
parser.add_option('-m','--mode', dest="mode", help='mode of analysis <prsice, manual>')
(options, args) = parser.parse_args()

# checking required arguments
for ra in required_arguments:
    if options.__dict__[ra] is None:
        parser.error('Error! Required argument {0} has not been provided.'.format(ra))

### create directory with the results
current_datetime = datetime.datetime.now().strftime("%I%M%p_%Y%B%d")
outfolder = "{0}/{1}".format(results_dir,current_datetime)
os.makedirs(outfolder, exist_ok=True)

### Choosing the type of analysis according to the mode argument
if (options.mode == "prsice"):
    # create plink files from submitted VCF
    os.system("{0} --vcf {1} --double-id --vcf-filter --out {2}/sample".format(plink, options.vcf,outfolder))
    # Assign control/test to fam file
    EditFam("{0}/sample.fam".format(outfolder))
    # Launch PRSice
    os.system("Rscript {0} --prsice {1} --dir {2} --base {3} \
            --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --pvalue P \
            --target {4}/sample --thread 1 --out {4}/prsice_analysis".format(prsice, prsice_bin, prsice_dir, options.gwas, outfolder))
    # Running twice in case of exec halted Error
    os.system("Rscript {0} --prsice {1} --dir {2} --base {3} \
            --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat OR --pvalue P \
            --target {4}/sample --thread 1 --out {4}/prsice_analysis --extract {4}/prsice_analysis.valid".format(prsice, prsice_bin, prsice_dir, options.gwas, outfolder))

elif (options.mode == "manual"):
    # Loading VCF file inside a structured variable
    vcfData = LoadVcf(options.vcf)
    # Loading GWAS association data inside a structured variable
    gwasData = LoadGwas(options.gwas)
    # Filtering VCF file according to custom parameters
    filteredVcfData = FilterVcf(vcfData["full"], "FILTER--NULL--PASS")
    # Intersecting VCF file mutations with GWAS association data
    commonSNPs = list(set(vcfData["slim"]) & set([x[:-1] for x in gwasData["slim"]]))
    # Calculate Polygenic Risk Score (for common SNPs)
    PRS = CalculatePRS(commonSNPs, gwasData["slim"], "snps")
    # Calculate Max PRS (for submitted GWAS association)
    maxPRS = CalculatePRS(gwasData["slim"], gwasData["slim"], "gwas")
    # Saving results (PRS for single SNPs and MAX)
    prs_detail_fn = "{0}/prs_detail.tsv".format(outfolder)
    max_prs_detail_fn = "{0}/max_prs_detail.tsv".format(outfolder)
    WritePRS(PRS["detail"],prs_detail_fn)
    WritePRS(maxPRS["detail"],max_prs_detail_fn)
    # Report PRS and maxPRS to the user
    print("\n\n \
    \t\t###### The PRS for the submitted VCF file is {0} \n\n \
    \t\t###### Maximum PRS for this GWAS file is {1}".format(PRS["score"],maxPRS["score"]))
    # Assessing the distance with a normal distribution and produce a plot
    os.system("Rscript scripts/calculate_significance.R --prs {0} --max {1} --dir {2}".format(PRS["score"], maxPRS["score"], outfolder))
