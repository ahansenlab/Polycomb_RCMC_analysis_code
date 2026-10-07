#Aligning and processing paired end fastq files for spike-in RNA-seq using bowtie2

from sys import exit
import subprocess as sp
import argparse
import multiprocessing
import uuid
import pandas as pd

parser = argparse.ArgumentParser(description = "run bowtie2 on paired end fastq files with spikein to produce aligned bam files")
parser.add_argument("--filename", "-f", help = "a tab-separated file containing on each line the path to the fastq with the first ends of pairs, the path to the fastq with the second ends of pairs, and the desired output name for the aligned file")
parser.add_argument("--genome", "-g", help = "genome build to align to - mouse or human - required")
parser.add_argument("--spikegenome", "-s", help = "spikein genome build to align to - mouse, human, or drosophila - required")
parser.add_argument("--threads", "-t", help = "number of threads to use for bowtie2 - default is 1", default = "1")
parser.add_argument("--outname", "-o", help = "name for the table to store counts", default = "alignmentcounts")
parser.add_argument("--outdir", help = "a directory to store output files - default is current directory", default = "./")
args = parser.parse_args()

genome = args.genome
spikegenome = args.spikegenome
threads = args.threads
outname = args.outname
outdir = args.outdir

#Read in file to determine what to process
files = pd.read_csv(args.filename, sep='\t', header=None, names = ['end1', 'end2', 'name'])
file1 = files.end1.values.tolist()
file2 = files.end2.values.tolist()
outnames = files.name.values.tolist()

#Check requirements are fulfilled:
condapacks = sp.run("conda list".split(), capture_output=True)
condapacksstr = str(condapacks.stdout)

if "bowtie2" not in condapacksstr or "samtools" not in condapacksstr or "sambamba" not in condapacksstr or "star" not in condapacksstr:
    print("Please make sure bowtie2, samtools, sambamba, and star are installed in your current conda environment (check conda list)")
    exit()

#Check that outdir ends with a /, add one if it doesn't
if args.outdir is not None and not outdir.endswith("/"):
    outdir = outdir + "/"

#Check a genome was specified
if args.genome is None:
    print("Genome not specified - check help for formatting")
    parser.print_usage()
    exit()

#Check if file1, file2, and outnames are present and the same lengths
if file1 is None or file2 is None or outnames is None:
    print("Input files not specified - check help for formatting")
    parser.print_usage()
    exit()
elif type(file1) == list and type(file1) == list and type(outnames) == list and len(file1) > 1:
    if len(file1) != len(file2) or len(file1) != len(outnames) or len(file2) != len(outnames):
        print("Mismatch in number of input files or output names, check arguments")
        exit()

#Check that a sensible number of threads has been requested - more protections here are possible - at the moment users are trusted to be sensible
cpucount = multiprocessing.cpu_count()
if args.threads is None:
    print("Defaulting to one thread")
    threads = 1
elif int(args.threads) >= cpucount:
    print("Too many threads requested, resetting to default")
    threads = 1
else:
    print(f"Running alignment with {threads} threads...")

genome_path = "/mnt/md0/DataRepository/genomes/"
STAR_genome_path = "/mnt/md0/DataRepository/STAR/"

#Check that the user has entered a valid genome to align to
if genome == "mm10" or genome == "mm39":
    print(f"Primary genome is mouse genome {genome}")
elif genome == "hg19" or genome == "hg38":
    print(f"Primary genome is human genome {genome}")
else:
    print("Genome option not recognised or not entered. Please use mm10/39 or hg19/38 or ask Miles to change the script to accommodate your new organism/genome")
    exit()

#Process ID (used to make unique sorttemp, so these are not overlapping for multiple processes in the same outdir)
uniqueid = str(uuid.uuid4())
sp.run(f"mkdir {outdir}{uniqueid}", shell=True, executable="/bin/bash")

#Create a place to store counts
allcountslist = []

print(f"Spike genome is {spikegenome}")

#Process the files

for i, fastq1 in enumerate(file1):
    fastq2 = file2[i]
    name = outnames[i]
    print(f"Aligning {name} to {genome} and {spikegenome}")
    #Explanation: align to rDNA from genome and spikein genome,
    # then remove multiply aligned reads with grep (XS: indicates the score of the next best aligning read, if it exists), then make a bam of all mapped reads (-F4 removes reads with a SAM flag of 4, which means unmapped)
    print("Removing rRNA reads...")
    sp.run(f"bowtie2 --very-fast -p {threads} --no-mixed --no-discordant -1 {fastq1} -2 {fastq2} -x {genome_path}{genome}.{spikegenome}.rDNA/{genome}.{spikegenome}.rDNA | samtools view -bSh -u -f 4 - | samtools sort -n -@ {threads} | samtools fastq -1 {outdir}{name}_1.fastq -2 {outdir}{name}_2.fastq -n", shell=True, executable="/bin/bash")
    #bamToFastq -i stdin -fq {outdir}{name}_1.fastq -fq2 {outdir}{name}_2.fastq", shell=True, executable="/bin/bash")
    print("Finished processing ribosomal DNA")

    print("Aligning ribosomal-depleted reads with STAR")
    sp.run(f"STAR --genomeDir {STAR_genome_path}{genome}.{spikegenome} --readFilesIn {outdir}{name}_1.fastq {outdir}{name}_2.fastq --runThreadN {threads} --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outFileNamePrefix {outdir}{name} --genomeLoad LoadAndKeep", shell=True, executable="/bin/bash")
    sp.run(f"samtools view -bhS -q 4 {outdir}{name}Aligned.out.sam > {outdir}{name}_STARmapped.bam", shell=True, executable="/bin/bash")
    sp.run(f"sambamba sort --tmpdir {outdir}{uniqueid}/ -t {threads} -m 30G {outdir}{name}_STARmapped.bam -o {outdir}{name}_STARmapped_sorted.bam", shell=True, executable="/bin/bash")
    print("Finished STAR alignment and filtering")
    sp.run(f"head -n 32 {outdir}{name}Log.final.out >> AlignmentStatistics.txt", shell=True, executable="/bin/bash")

    print("Now aligning with bowtie2")
    sp.run(f"bowtie2 --sensitive-local -p {threads} --no-mixed --no-discordant --mm -x {genome_path}{genome}.{spikegenome}/{genome}.{spikegenome} -1 {outdir}{name}Unmapped.out.mate1 -2 {outdir}{name}Unmapped.out.mate2 | grep -v XS: - | samtools view -bhS -F4 - > {outdir}{name}_Bowtie2mapped.bam", shell=True, executable="/bin/bash")

    sp.run(f"sambamba sort --tmpdir {outdir}{uniqueid}/ -t {threads} -m 30G {outdir}{name}_Bowtie2mapped.bam -o {outdir}{name}_Bowtie2mapped_sorted.bam", shell=True, executable="/bin/bash")
    print("Finished bowtie2 alignment")

    #Merge reads aligned with STAR and bowtie2, sort and remove duplicates
    sp.run(f"sambamba merge -t {threads} {outdir}{name}_mapped.bam {outdir}{name}_Bowtie2mapped_sorted.bam {outdir}{name}_STARmapped_sorted.bam", shell=True, executable="/bin/bash")
    sp.run(f"sambamba sort --tmpdir {outdir}{uniqueid}/ -t {threads} -m 30G {outdir}{name}_mapped.bam -o {outdir}{name}_mapped_sorted.bam", shell=True, executable="/bin/bash")
    sp.run(f"sambamba markdup -r --tmpdir {outdir}{uniqueid}/ -t {threads} --hash-table-size 1000000 --overflow-list-size 2000000 {outdir}{name}_mapped_sorted.bam {outdir}{name}_mapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")

    #Next, need to separate out reads from each genome:
    print(f"Extracting reads aligning uniquely to {genome}")
    #Use samtools view to open the file, grep to remove those with the spikein genome name in them, then change the chromosomes with the genome name to just the chromosome numbers as normally used, then put the file back into bam.
    sp.run(f"samtools view -@ {threads} -h {outdir}{name}_mapped_sorted_rmdup.bam | grep -v {spikegenome} | sed s/{genome}_chr/chr/g | samtools view -@ {threads} -bhS - > {outdir}{name}_{genome}.mapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
    #Do the same thing for the spikein
    print(f"Extracting reads aligning uniquely to {spikegenome}.")
    sp.run(f"samtools view -@ {threads} -h {outdir}{name}_mapped_sorted_rmdup.bam | grep -v {genome} | sed s/{spikegenome}_chr/chr/g | samtools view -@ {threads} -bhS - > {outdir}{name}_{spikegenome}.mapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
    #Index outputs
    sp.run(f"sambamba index -t {threads} {outdir}{name}_mapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
    sp.run(f"sambamba index -t {threads} {outdir}{name}_{genome}.mapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")
    sp.run(f"sambamba index -t {threads} {outdir}{name}_{spikegenome}.mapped_sorted_rmdup.bam", shell=True, executable="/bin/bash")

    #Now clean up
    sp.run(f"rm {outdir}{name}_1.fastq {outdir}{name}_2.fastq {outdir}{name}Aligned.out.sam {outdir}{name}Unmapped.out.mate1 {outdir}{name}Unmapped.out.mate2 {outdir}{name}Log.progress.out {outdir}{name}SJ.out.tab {outdir}{name}_Bowtie2mapped.bam {outdir}{name}_Bowtie2mapped_sorted.bam {outdir}{name}_STARmapped.bam {outdir}{name}_STARmapped_sorted.bam {outdir}{name}_mapped.bam {outdir}{name}_mapped_sorted.bam", shell=True, executable="/bin/bash")

    #Finally, count reads in each file:
    totalCount = sp.run(f"sambamba view -c -t {threads} {outdir}{name}_mapped_sorted_rmdup.bam", capture_output=True, shell=True, executable="/bin/bash")
    genomeCount = sp.run(f"sambamba view -c -t {threads} {outdir}{name}_{genome}.mapped_sorted_rmdup.bam", capture_output=True, shell=True, executable="/bin/bash")
    spikegenomeCount = sp.run(f"sambamba view -c -t {threads} {outdir}{name}_{spikegenome}.mapped_sorted_rmdup.bam", capture_output=True, shell=True, executable="/bin/bash")
    countsList = [totalCount.stdout.decode('ascii').strip(), genomeCount.stdout.decode('ascii').strip(), spikegenomeCount.stdout.decode('ascii').strip()]
    allcountslist.append(countsList)

countstable = pd.DataFrame(allcountslist, columns = ['allcounts', 'genomecounts', 'spikecounts'])
outtable = pd.concat([files, countstable], axis=1)
outtable.to_csv(outdir + outname + '.tsv', sep = '\t', index = False, header = True)

sp.run(f"rmdir {outdir}{uniqueid}", shell=True, executable="/bin/bash")