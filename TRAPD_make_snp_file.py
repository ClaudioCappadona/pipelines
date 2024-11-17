#!/usr/bin/python
import optparse
import operator
import re
import sys
import gzip
import bisect
import pandas as pd
import re
import numpy as np
import pickle

# Parse options
parser = optparse.OptionParser()
parser.add_option("-v", "--vcffile", action="store", dest="vcffilename")
parser.add_option("-o", "--outfile", action="store",
                  dest="outfilename", default="snpfile.txt")
parser.add_option("--genecolname", action="store", dest="genecolname")
parser.add_option("--includesamples", action="store", dest="samplesfilename")

# Filters
parser.add_option("--includeinfo", action="append", dest="includeinfo")
parser.add_option("--excludeinfo", action="append", dest="excludeinfo")
parser.add_option("--includevep", action="append", dest="includevep")
parser.add_option("--excludevep", action="append", dest="excludevep")
parser.add_option("--pass", action="store_true", dest="passfilter")
parser.add_option("--vep", action="store_true", dest="vep")
parser.add_option("--snponly", action="store_true", dest="snponly")
parser.add_option("--indelonly", action="store_true", dest="indelonly")
parser.add_option("--bedfile", action="store", dest="bedfilename")

parser.add_option("--snpformat", action="store",
                  dest="snpformat", default="EPACTS")
parser.add_option("--genenull", action="store",
                  dest="genenull", default=".,NA")

options, args = parser.parse_args()

# Try to catch potential errors for parsing parameters
if not options.vcffilename:   # if filename is not given
    parser.error('A vcf file is needed')
    sys.exit()

if (options.vcffilename.endswith(".gz") is False) and (options.vcffilename.endswith(".bgz") is False):   # if vcf filename is not given
    parser.error('Is your vcf file gzipped?')
    sys.exit()

if not options.genecolname:
    parser.error('An INFO field with the gene names to use must be provided')
    sys.exit()

if (options.includevep is not None) or (options.excludevep is not None):
    if not options.vep:
        parser.error('--vep option must be supplied if using VEP annotations')
        sys.exit()

if options.snpformat != "VCFID" and options.snpformat != "CHRPOSREFALT" and options.snpformat != "EPACTS":   # if filename is not given
    parser.error('SNP format must be "VCFID" or "CHRPOSREFALT" or "EPACTS"')
    sys.exit()

if options.snponly and options.indelonly:
    parser.error('Please select only --snponly or --indelonly')
    sys.exit()

# Check to make sure all the filters seem well formed


def checkfilter(infofilter):
    if ("[" not in infofilter) or (infofilter.startswith("]")) or (infofilter.endswith("]")) or str(infofilter.split("[")[1].split("]")[0]) not in ["<", ">", "<=", ">=", "=", "!=", "in", "%"]:
        return 0
    else:
        return 1

# Filter to make sure that values are either all numeric or all str


def consistent(option_value, field_value):
    try:
        float(option_value) or int(option_value)
        try:
            float(field_value) or int(field_value)
            field_out = float(field_value)
            option_out = float(option_value)
            c = 1
        except ValueError:
            option_out = option_value
            field_out = field_value
            c = 0
    except ValueError:
        field_out = str(field_value)
        option_out = str(option_value)
        c = 1
    return [option_out, field_out, c]


# Read in vcf header and extract all INFO fields
info_fields = []
chrformat = "number"
vcffile = gzip.open(options.vcffilename, 'rt', encoding='utf-8')
for line_vcf1 in vcffile:
    if line_vcf1[0] == "#":
        if "##INFO=<ID=" in line_vcf1:
            temp_field = line_vcf1.split("##INFO=<ID=")[1].split(",")[0]
            info_fields.append(temp_field)
        elif "##contig" in line_vcf1:
            if "ID=chr" in line_vcf1:
                chrformat = "chr"
    else:
        break
vcffile.close()

# Read in vcf header to get VEP CSQ fields
if options.vep:
    #vcffile = gzip.open(options.vcffilename, "rb")
    vcffile = gzip.open(options.vcffilename, 'rt', encoding='utf-8')
    csq_found = 0
    for line_vcf1 in vcffile:
        if line_vcf1[0] == "#":
            if ("INFO=<ID=CSQ" in line_vcf1) or ("ID=CSQ" in line_vcf1) or ("ID=vep" in line_vcf1):
                csq_anno = line_vcf1.rstrip('\n').replace(
                    '"', '').strip('>').split("Format: ")[1].split("|")
                csq_found = 1
                break
    if csq_found == 0:
        sys.stdout.write("VEP CSQ annotations not found in vcf header\n")
        sys.exit()
    vcffile.close()

if options.vep:
    if options.genecolname not in csq_anno:
        sys.stdout.write("Gene column name not found in VEP annotations\n")
        sys.exit()


# Run through all filters to make sure they're okay
if options.includeinfo is not None:
    for i in range(0, len(options.includeinfo), 1):
        if checkfilter(options.includeinfo[i]) == 0:
            sys.stdout.write(str(options.includeinfo[i])+" is malformed\n")
            sys.exit()
        if options.includeinfo[i].split("[")[0] not in info_fields:
            sys.stdout.write(
                str(options.includeinfo[i])+" is not in VCF file\n")
            sys.exit()
if options.excludeinfo is not None:
    for i in range(0, len(options.excludeinfo), 1):
        if checkfilter(options.excludeinfo[i]) == 0:
            sys.stdout.write(str(options.excludeinfo[i])+" is malformed\n")
            sys.exit()
        if options.excludeinfo[i].split("[")[0] not in info_fields:
            sys.stdout.write(
                str(options.excludeinfo[i])+" is not in VCF file\n")
            sys.exit()
if options.includevep is not None:
    for i in range(0, len(options.includevep), 1):
        if checkfilter(options.includevep[i]) == 0:
            sys.stdout.write(str(options.includevep[i])+" is malformed\n")
            sys.exit()
        if options.includevep[i].split("[")[0] not in csq_anno:
            sys.stdout.write(
                str(options.includevep[i])+" is not in VCF file\n")
            sys.exit()
if options.excludevep is not None:
    for i in range(0, len(options.excludevep), 1):
        if checkfilter(options.excludevep[i]) == 0:
            sys.stdout.write(str(options.excludevep[i])+" is malformed\n")
            sys.exit()
        if options.excludevep[i].split("[")[0] not in csq_anno:
            sys.stdout.write(
                str(options.excludevep[i])+" is not in VCF file\n")
            sys.exit()


# Test if something is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# Extract canonical vep annotation


def canonical_vep(vcfline):
    canonical_index = csq_anno.index("CANONICAL")
    out = ""
    if (vcfline.split("|")[canonical_index]).upper() == "YES":
        out = 1
    return out


def find_vep_gene(genecolname, annot, csq_anno):
    csq_index = csq_anno.index(genecolname)
    genename = annot.split("|")[csq_index]
    return genename


def find_info_gene(genecolname, vcfline):
    if genecolname in vcfline:
        genename = (";"+vcfline).split(";"+genecolname+"=")[1].split(";")[0]
    else:
        genename = ""
    return genename


def test_include_info(filter, vcfline):
    option_field = filter.split("[")[0]
    option_value = filter.split("]")[1]
    if (";"+option_field+"=") in (";"+vcfline):
        field_value = (";"+vcfline).split((";"+option_field+"=")
                                          )[1].split(";")[0].split(",")[0]
        consist_out = consistent(option_value, field_value)
        if consist_out[2] == 1:
            if filter.split("[")[1].split("]")[0] == "in":
                listvalues = option_value.lstrip("(").rstrip(")").split(',')
                counter = 0
                for i in range(0, len(listvalues), 1):
                    if operator.eq(field_value, listvalues[i]):
                        counter += 1
                if counter > 0:
                    return 1
                else:
                    return 0
            else:
                if get_operator_fn(filter.split("[")[1].split("]")[0])(consist_out[1], consist_out[0]):
                    return 1
                else:
                    return 0
        else:
            return 0
    else:
        return 0


def test_exclude_info(filter, vcfline):
    option_field = filter.split("[")[0]
    option_value = filter.split("]")[1]
    if (";"+option_field+"=") in (";"+vcfline):
        field_value = (";"+vcfline).split((";"+option_field+"=")
                                          )[1].split(";")[0].split(",")[0]
        consist_out = consistent(option_value, field_value)
        if consist_out[2] == 1:
            if filter.split("[")[1].split("]")[0] == "in":
                listvalues = option_value.lstrip("(").rstrip(")").split(',')
                counter = 0
                for i in range(0, len(listvalues), 1):
                    if operator.eq(field_value, listvalues[i]):
                        counter += 1
                if counter > 0:
                    return 0
                else:
                    return 1
            else:
                if get_operator_fn(filter.split("[")[1].split("]")[0])(consist_out[1], consist_out[0]):
                    return 0
                else:
                    return 1
        else:
            return 0
    else:
        return 0


def test_include_vep(filter, annot, csq_anno):
    option_field = filter.split("[")[0]
    csq_index = csq_anno.index(option_field)
    option_value = filter.split("]")[1]
    field_value = annot.split("|")[csq_index]
    consist_out = consistent(option_value, field_value)
    if consist_out[2] == 1:
        if filter.split("[")[1].split("]")[0] == "in":
            listvalues = option_value.lstrip("(").rstrip(")").split(',')
            counter = 0
            for i in range(0, len(listvalues), 1):
                if operator.eq(field_value, listvalues[i]):
                    counter += 1
            if counter > 0:
                return 1
            else:
                return 0
        else:
            if get_operator_fn(filter.split("[")[1].split("]")[0])(consist_out[1], consist_out[0]):
                return 1
            else:
                return 0
    else:
        return 0


def test_exclude_vep(filter, annot, csq_anno):
    option_field = filter.split("[")[0]
    csq_index = csq_anno.index(option_field)
    option_value = filter.split("]")[1]
    field_value = annot.split("|")[csq_index]
    consist_out = consistent(option_value, field_value)
    if consist_out[2] == 1:
        if filter.split("[")[1].split("]")[0] == "in":
            listvalues = option_value.lstrip("(").rstrip(")").split(',')
            counter = 0
            for i in range(0, len(listvalues), 1):
                if operator.eq(field_value, listvalues[i]):
                    counter += 1
            if counter > 0:
                return 0
            else:
                return 1
        else:
            if get_operator_fn(filter.split("[")[1].split("]")[0])(consist_out[1], consist_out[0]):
                return 0
            else:
                return 1
    else:
        return 0


# Function to match operator strings
def get_operator_fn(op):
    return {
        '<': operator.lt,
        '<=': operator.le,
        '>': operator.gt,
        '>=': operator.gt,
        '=': operator.eq,
        '!=': operator.ne,
        '%': operator.contains,
    }[op]


# Create empty snptable
snptable = {}

#read in bedfile
if options.bedfilename is not None:
    if str(options.bedfilename).endswith(".gz") is True:
        bed = gzip.open(options.bedfilename, "rb")
    else:
        bed = open(options.bedfilename, "r")
    bed_lower = {}
    bed_upper = {}
    for line_b1 in bed:
        line_b = line_b1.rstrip().split('\t')
        chr = str(line_b[0]).lower().replace("chr", "")
        if chr not in bed_lower:
            bed_lower[chr] = [chr, []]
            bed_upper[chr] = [chr, []]
        bed_lower[chr][1].append(int(line_b[1])+1)
        bed_upper[chr][1].append(int(line_b[2]))
    bed.close()

#vcffile = gzip.open(options.vcffilename, "rb")
vcffile = gzip.open(options.vcffilename, 'rt', encoding='utf-8')

for line_vcf1 in vcffile:
    line_vcf = line_vcf1.rstrip().split('\t')
    keep = 1

    if keep == 1 and options.samplesfilename:
        if line_vcf[0]== "#CHROM":
            header = line_vcf
            header = [r.strip() for r in header] #remove newline
            header = [i for i in header if i] # remove any empty element
            samples = open(options.samplesfilename, "rt", encoding='utf-8')
            samples=samples.read().split("\n")
            samples = [r.strip() for r in samples] #remove newline
            samples = [i for i in samples if i] # remove any empty element

            rev_intersec_samples=[ element for element in samples if element not in header[9:]]
            match_samples=list(set(header[9:]).intersection(samples))

            index=np.nonzero(np.in1d(header, samples))[0]
            #print(np.array(header)[index.astype(int)])
            if(len(samples) == len(match_samples)) and (len(rev_intersec_samples) == 0):
                sys.stdout.write(
                    "all samples in"+ str(options.samplesfilename.split("/")[-1])+ "are present in vcf file...\n")
                sys.stdout.write(
                    "total samples:"+ str(len(match_samples))+"\n")
            elif((len(samples) != len(match_samples)) and (len(rev_intersec_samples) != 0)):
                sys.stdout.write(
                    "warning: these samples from"+ str(options.samplesfilename.split("/")[-1])+ "are not present in vcf file:\n")
                sys.stdout.write(
                    str("\n".join(rev_intersec_samples)))
            elif(len(match_samples)==0):
                sys.stdout.write(
                    "error: no samples in"+ str(options.samplesfilename)+ "are present in vcf file\n")
                sys.exit()
                    
    
    if line_vcf[0][0] != "#":
        if keep == 1 and options.samplesfilename:
            indexed_line=np.array(line_vcf)[index.astype(int)] #get index of
            genotypes=[r.split(":")[0] for r in indexed_line] #get format
            genotypes=[r.split("/") for r in genotypes]
            genotypes=sum(genotypes, [])
            genotypes = [x for x in genotypes if x != "."] 
            if(sum([eval(i) for i in genotypes])==0):
                keep = 0 #ensure no row with total vars=0 is retained after samples filter
            else:
                cols=np.arange(start=0, stop=8, step=1)
                whole_index=np.concatenate((cols, index.astype(int)), axis=None)
                line_vcf=np.array(line_vcf)[whole_index.astype(int)] #new line with new samples subset
                
        if keep == 1 and options.passfilter:
            if line_vcf[6] != "PASS":
                keep = 0
        if keep == 1 and options.snponly:
            if len(line_vcf[3]) > 1 or len(line_vcf[4]) > 1:
                keep = 0
        if keep == 1 and options.indelonly:
            if len(line_vcf[3]) == 1 and len(line_vcf[4]) == 1:
                keep = 0
        if "," in line_vcf[4]:
            keep = 0

        # Subset on bedfile
        if options.bedfilename is not None:
            chr = str(line_vcf[0]).lower().replace("chr", "")
            temp_index = bisect.bisect(bed_lower[chr][1], int(line_vcf[1]))-1
            if temp_index < 0:
                keep = 0
            elif int(line_vcf[1]) > bed_upper[chr][1][temp_index]:
                keep = 0

 # Go through INFO field filters
        if keep == 1 and options.includeinfo is not None:
            iter = 0
            while keep == 1 and iter < len(options.includeinfo):
                filter = options.includeinfo[iter]
                keep = test_include_info(filter, line_vcf[7])
                iter = iter+1

        if keep == 1 and options.excludeinfo is not None:
            iter = 0
            while keep == 1 and iter < len(options.excludeinfo):
                filter = options.excludeinfo[iter]
                keep = test_exclude_info(filter, line_vcf[7])
                iter = iter+1

   # Go through INFO/VEP field filters
        if keep == 1 and options.vep:
            vcfline = line_vcf[7].replace("vep=", "CSQ=")
            if "CSQ=" in vcfline:
                annots = (
                    ";"+vcfline).split(";CSQ=")[1].split(";")[0].split(",")
                keep_a = [1] * len(annots)

                for i in range(0, len(annots), 1):
                    if len(csq_anno) == len(annots[i].split("|")) and canonical_vep(annots[i]) == 1:
                        if options.includevep is not None:
                            iter = 0
                            while keep_a[i] == 1 and iter < len(options.includevep):
                                filter = options.includevep[iter]
                                keep_a[i] = test_include_vep(
                                    filter, annots[i], csq_anno)
                                iter = iter+1
                        if options.excludevep is not None:
                            iter = 0
                            while keep_a[i] == 1 and iter < len(options.includevep):
                                filter = options.includevep[iter]
                                keep_a[i] = test_include_vep(
                                    filter, annots[i], csq_anno)
                                iter = iter+1
                    else:
                        keep_a[i] = 0
                if not 1 in keep_a:
                    keep = 0
# If variant meets all filters for at least one transcript, then extract gene name for all ok transcripts
        if keep == 1:
            vcfline = line_vcf[7].replace("vep=", "CSQ=")
            if options.vep and "CSQ=" in vcfline:
                gene = []
                for i in range(0, len(annots), 1):
                    if keep_a[i] == 1:
                        gene.append(find_vep_gene(
                            options.genecolname, annots[i], csq_anno))
            else:
                gene = []
                gene.append(find_info_gene(options.genecolname, line_vcf[7]))
            gene = list(set(gene))
            if len(gene) > 0:
                if options.snpformat == "VCFID":
                    snpid = str(line_vcf[2])
                    
                elif options.snpformat == "EPACTS":
                    snpid = chrformat+str(line_vcf[0].lstrip(
                        "chr"))+":"+str(line_vcf[1])+"_"+str(line_vcf[3])+"/"+str(line_vcf[4])
                else:
                    snpid = str(line_vcf[0].lstrip(
                       "chr" ))+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
                for i in range(0, len(gene), 1):
                    if gene[i] not in options.genenull.split(","):
                        if gene[i] not in snptable:
                            snptable[gene[i]] = [gene[i], [snpid]]
                        else:
                            snptable[gene[i]][1].append(snpid)

vcffile.close()


###################
# Write Output
outfile = open(options.outfilename, "w")
if options.snpformat != "EPACTS":
    outfile.write("#GENE\tSNPS\n")
#for genes in snptable:   
for x in snptable:
    if len(x) > 0:

        # Read through hash table and print out variants
        if options.snpformat != "EPACTS":
            snp_out = ','.join(snptable[x][1])
        elif options.snpformat == "EPACTS":
            df=pd.DataFrame()
            for ele in snptable[x][1]:
                ser=ele.replace("chr", "")
                ser=ser.replace("X", "23")
                ser=pd.Series(re.split('_|:',ser))
                ser=ser.to_frame().T
                ser[0]=np.array(ser[0].astype(int))
                ser[1]=np.array(ser[1].astype(int))
                df=df.append(ser,ignore_index=True)
            df.columns=["chr", "pos", "vars"]
            df = df.sort_values(['chr','pos'], ascending=True)
            c = [ snptable[x][1][i] for i in df.index]
            #ls=[]
            #for i in df.index:
            #    ls.append(chrformat+str(df.loc[i][0]).replace("23", "X")+":"+str(df.loc[i][1])+"_"+str(df.loc[i][2]))
            #c = [ snptable[x][1][i] for i in df.index]
            snp_out = '\t'.join(c)
            outfile.write(str(x)+"\t"+snp_out+"\n")
outfile.close()

# python make_snp_file.py -o test.out.txt -v gnomad.test.vcf.gz --vep --genecolname SYMBOL --snpformat CHRPOSREFALT --pass --includeinfo "AC[<]5"
