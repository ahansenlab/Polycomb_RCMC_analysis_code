#Script to make a bunch of bigwigs using bam files with names specified in a text file

from sys import exit
import subprocess as sp
import argparse
import uuid
import re

parser = argparse.ArgumentParser(description = "make BigWigs from bam files")
parser.add_argument("filename", help = "a file containing the names of all bams from which to make BigWigs")
parser.add_argument("readtype", help = "process paired or single end read data (only paired is supported in RNA mode currently)")
parser.add_argument("genome", help = "genome to get genome sizes for - mouse or human")
parser.add_argument("--threads", "-t", help = "number of threads to use - default is 1 - only relevant for RNA mode", default = "1")
parser.add_argument("--outdir", help = "a directory to store output files", default = "./")
parser.add_argument("-r", "--rna", help = "set this flag to make strand-specific bigwigs from RNA-seq data and similar", action = "store_true")
args = parser.parse_args()

filename = args.filename
readtype = args.readtype
genome = args.genome
outdir = args.outdir
rnaflag = args.rna
threads = args.threads

#Check requirements are fulfilled:
condapacks = sp.run("conda list".split(), capture_output=True)
condapacksstr = str(condapacks.stdout)
if "macs" not in condapacksstr:
    print("Please make sure macs3 is installed in your current conda environment (check conda list)")
    exit()

path_to_chromsizes = "/mnt/md0/DataRepository/chromsizes/"
path_to_wigtobigwig = "/mnt/md0/Tools/scripts/wigToBigWig"

if genome == "mm10" or genome == "mm39":
    gensize = "mm"
elif genome == "hg19" or genome == "hg38":
    gensize = "hs"
else:
    print("Genome option not recognised. Please use mm10/39 or hg19/38 or change the script to accommodate your new organism/genome.")
    exit()



bams_to_process_file = open(filename, "r")
bamfiles_to_process = bams_to_process_file.read().splitlines()
nfiles = len(bamfiles_to_process)
print("Making bigwigs for", nfiles, "file(s).")

# names = [bams.rsplit("/",1)[1].rsplit(".",1)[0] for bams in bamfiles_to_process]
# names = [bams.rsplit(".",1)[0] for bams in bamfiles_to_process]
names = [None] * len(bamfiles_to_process)
for count, bams in enumerate(bamfiles_to_process):
    if '/' in bams:
        # path = bams.rsplit("/",1)[0]
        name = bams.rsplit("/",1)[1].rsplit(".",1)[0]
    else:
        # path = "./"
        name = bams.rsplit(".",1)[0]
    names[count] = name


if args.outdir is not None and not outdir.endswith("/"):
    outdir = outdir + "/"
elif args.outdir is not None and " " in outdir:
    outdir = re.escape(outdir)


#Process ID (used to make unique sorttemp, so these are not overlapping for multiple processes in the same outdir)
uniqueid = str(uuid.uuid4())

if not rnaflag:
    if readtype == "single" or readtype == "Single":
        for i, bam in enumerate(bamfiles_to_process):
            if " " in names[i]:
                name = re.escape(names[i])
            else:
                name = names[i]
            #MACS3
            sp.run(f"macs3 pileup -i {bam} -f BAM -o {name}.bg --outdir {outdir} --verbose 0", shell=True, executable="/bin/bash")
            #Wigtobigwig
            sp.run(f"{path_to_wigtobigwig} -clip {outdir}{name}.bg {path_to_chromsizes}{genome}/{genome}.chrom.sizes {outdir}{name}.bw", shell=True, executable="/bin/bash")
            sp.run(f"rm {outdir}{name}.bg", shell=True, executable="/bin/bash")
            print(name, "finished!")
    elif readtype == "paired" or readtype == "Paired":
        #Make temp directory
        sp.run(f"mkdir {outdir}{uniqueid}/", shell=True, executable="/bin/bash")
        for i, bam in enumerate(bamfiles_to_process):
            if " " in names[i]:
                name = re.escape(names[i])
            else:
                name = names[i]
            print(f"Making bigwig for {name}")
            #MACS3
            sp.run(f"macs3 callpeak -t {bam} -f BAMPE -g {gensize} --bdg -n {name} --outdir {outdir} --verbose 0 --tempdir {outdir}{uniqueid}/", shell=True, executable="/bin/bash")
            #Wigtobigwig
            sp.run(f"{path_to_wigtobigwig} -clip {outdir}{name}_treat_pileup.bdg {path_to_chromsizes}{genome}/{genome}.chrom.sizes {outdir}{name}.bw", shell=True, executable="/bin/bash")
            sp.run(f"rm {outdir}{name}_control_lambda.bdg {outdir}{name}_peaks.narrowPeak {outdir}{name}_peaks.xls {outdir}{name}_summits.bed {outdir}{name}_treat_pileup.bdg", shell=True, executable="/bin/bash")
            print(bam, "finished!")
        sp.run(f"rm {outdir}{uniqueid}/ -r", shell=True, executable="/bin/bash")
    else:
        print("Please specify 'single' or 'paired' reads")
        exit()
elif rnaflag:
    if readtype == "single" or readtype == "Single":
        print("I'm sorry, I can't process single-end RNA-seq data yet.")
        exit()
    elif readtype == "paired" or readtype == "Paired":
        #Make temp directory (not needed here?)
        # sp.run(f"mkdir {outdir}{uniqueid}/", shell=True, executable="/bin/bash")
        for i, bam in enumerate(bamfiles_to_process):
            if " " in names[i]:
                name = re.escape(names[i])
            else:
                name = names[i]
            print(f"Making bigwig for {name}")
            #Separate strands
            #Forward
            sp.run(f"samtools view -b -f 83 {bam} > {outdir}{name}_83.bam", shell=True, executable="/bin/bash")
            sp.run(f"samtools view -b -f 163 {bam} > {outdir}{name}_163.bam", shell=True, executable="/bin/bash")
            #Reverse
            sp.run(f"samtools view -b -f 147 {bam} > {outdir}{name}_147.bam", shell=True, executable="/bin/bash")
            sp.run(f"samtools view -b -f 99 {bam} > {outdir}{name}_99.bam", shell=True, executable="/bin/bash")
            #Merge into two files
            sp.run(f"sambamba merge -t {threads} {outdir}{name}_FORWARD.bam {outdir}{name}_83.bam {outdir}{name}_163.bam", shell=True, executable="/bin/bash")
            sp.run(f"sambamba merge -t {threads} {outdir}{name}_REVERSE.bam {outdir}{name}_147.bam {outdir}{name}_99.bam", shell=True, executable="/bin/bash")
            #Index
            sp.run(f"sambamba index -t {threads} {outdir}{name}_FORWARD.bam", shell=True, executable="/bin/bash")
            sp.run(f"sambamba index -t {threads} {outdir}{name}_REVERSE.bam", shell=True, executable="/bin/bash")
            #Remove the intermediates
            sp.run(f"rm {outdir}{name}_83.bam {outdir}{name}_163.bam {outdir}{name}_147.bam {outdir}{name}_99.bam", shell=True, executable="/bin/bash")
            #Make bigwigs
            sp.run(f"genomeCoverageBed -bga -split -ibam {outdir}{name}_FORWARD.bam | {path_to_wigtobigwig} stdin {path_to_chromsizes}{genome}/{genome}.chrom.sizes {outdir}{name}_FORWARD.bw", shell=True, executable="/bin/bash")
            sp.run(f"genomeCoverageBed -bga -split -ibam {outdir}{name}_REVERSE.bam | {path_to_wigtobigwig} stdin {path_to_chromsizes}{genome}/{genome}.chrom.sizes {outdir}{name}_REVERSE.bw", shell=True, executable="/bin/bash")
            print(bam, "finished!")
        # sp.run(f"rm {outdir}{uniqueid}/ -r", shell=True, executable="/bin/bash")
    else:
        print("Please specify 'single' or 'paired' reads")
        exit()
