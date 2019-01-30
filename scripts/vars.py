#/usr/bin/env python
# Coder: -.- Dr. Stefano Pirro (aka wynstep)

# Library of custom vars necessary for PRSCalc

#arguments required for performing the analysis
required_arguments = ["vcf", "mode", "gwas"]

# folder with results
results_dir = "results"

# plink exec
plink = "./tools/plink/plink"
# prsice exec
prsice = "./tools/prsice/PRSice.R"
prsice_bin = "./tools/prsice/PRSice_linux"
prsice_dir = "./tools/prsice/"
