#/usr/bin/env python
# Coder: -.- Dr. Stefano Pirro (aka wynstep)

# Library of custom functions necessary for PRSCalc
import gzip

# this function loads the vcf file into a structured variable
# The VCF file must respect the VCFv4.0+ file format
def LoadVcf(vcf_fn):
    # check the extension of the file
    if vcf_fn.endswith('.gz'):
        # if it's a gunzipped file, we take advantage of the gzip module to open it while extracting
        vcf_io = gzip.open(vcf_fn, "rt")
    else:
        # otherwise, we just parse it
        vcf_io = open(vcf_fn, "r")

    # iterating into the VCF file
    vcf = []
    for line in vcf_io:
        # check snp header
        if line.startswith('#CHR'):
            header = line.rstrip().replace("#","").split("\t")
        # check lines with snp information
        elif line.startswith('chr'): # Here we filter out the header of the VCF file
            line = line.rstrip().split('\t')
            dict_snp = dict(zip(header,line))
            # adding dict snp to vcf list
            vcf.append(dict_snp)

    ### Slim vcf file for comparison with GWAS
    slimmed_vcf = []
    for vcf_el in vcf:
        slimmed_vcf.append((vcf_el["CHROM"].replace("chr",""), \
                            vcf_el["POS"], vcf_el["ID"], vcf_el["REF"], vcf_el["ALT"]))

    return {"full": vcf, "slim": slimmed_vcf}

def LoadGwas(gwas_fn):
    # initialising the GWAS association list
    gwas = []
    with open(gwas_fn, "r") as gwas_io:
        # collecting header
        header = gwas_io.readline().rstrip().split("\t")
        for line in gwas_io:
            line = line.rstrip().split('\t')
            dict_gwas = dict(zip(header,line))
            # adding dict snp to gwas list
            gwas.append(dict_gwas)

    ### Slim gwas file for comparison with vcf
    slimmed_gwas = []
    for gwas_el in gwas:
        slimmed_gwas.append((gwas_el["CHR"], \
                            gwas_el["BP"], gwas_el["SNP"], gwas_el["A1"], gwas_el["A2"], gwas_el["OR"]))

    return {"full": gwas, "slim": slimmed_gwas}

def FilterVcf(vcf_data, filtering_string):
    # iterating into vcf elements (SNPs) and discard the items not satisfying the required_arguments
    # filtering string --> vcf_header--vcf_subheader(if any)--value,
    filtering_container = [value.split('--') for value in filtering_string.split(',')]
    # iterating into vcf data
    for fc in filtering_container:
        tmp_filtered_snps = []
        for vcf_el in vcf_data:
            if fc[0] in vcf_el.keys():
                if (fc[1] == "NULL" and fc[2] == vcf_el[fc[0]]):
                    tmp_filtered_snps.append(vcf_el)
        # append tmp filtered snps to bigger container
        try:
            filtered_snps = {x:filtered_snps[x] for x in filtered_snps if x in tmp_filtered_snps}
        except UnboundLocalError:
            filtered_snps = tmp_filtered_snps

    return filtered_snps

def CalculatePRS(slimFile, slimGWAS, calcType):
    # iterating into slim file and get the OR value (from the gwas)
    prs_score = 0
    prs_score_detail = ["CHR\tPOS\tID\tREF\tALT\tPRS_SCORE\n"]
    for index, slim_el in enumerate(slimFile):
        if calcType == "snps":
            if slim_el in [x[:-1] for x in slimGWAS]:
                prs_score_snp = len(slim_el[4].split(","))*float(slimGWAS[index][5])
                prs_score += prs_score_snp
                prs_score_detail.append("{0}\t{1}\n".format("\t".join(slim_el),prs_score_snp))
        else:
            if slim_el in slimGWAS:
                prs_score_snp = len(slim_el[4].split(","))*float(slimGWAS[index][5])
                prs_score += prs_score_snp
                prs_score_detail.append("{0}\t{1}\n".format("\t".join(slim_el),prs_score_snp))

    return {"score":prs_score, "detail":prs_score_detail}

def WritePRS(prs_score_detail, outFile):
    io_stream = open(outFile, "w")
    for prs_det in prs_score_detail:
        io_stream.write(prs_det)
    io_stream.close()

def EditFam(fam_file):
    new_fam = []
    with open(fam_file, "r") as fam_io:
        for index,line in enumerate(fam_io):
            if index % 2 == 0:
                new_fam.append(line.replace("-9", "0"))
            else:
                new_fam.append(line.replace("-9", "1"))
    fam_io.close()

    # overwrite new fam
    io_stream = open(fam_file, "w")
    for el in new_fam:
        io_stream.write(el)
    io_stream.close
