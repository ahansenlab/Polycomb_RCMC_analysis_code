#Analysing Polycomb loop strength

##################Imports###################
library(ggsignif)
library("tidyr")
library('plyr')
library('dplyr')
require('ggplot2')
require('reshape2')
library(ggExtra)
library('purrr')
library('stringr')
library(InteractionSet)

###########Read in data#############
#Loop strengths
#Rep 1
rep1.path <- "C:/Users/newto/Hansen Lab Dropbox/Members/Miles/HTS Data/Polycomb/loops/PcG_RCMC_rep1_PcGloopstrength.bedpe"
rep1 <- read.delim(rep1.path, header = TRUE, sep = "\t", dec = ".")

#Rep 2
rep2.path <- "C:/Users/newto/Hansen Lab Dropbox/Members/Miles/HTS Data/Polycomb/loops/PcG_RCMC_rep2_PcGloopstrength.bedpe"
rep2 <- read.delim(rep2.path, header = TRUE, sep = "\t", dec = ".")

##################Initial data organisation#############
#Merge the replicate data
merged.reps <- left_join(rep1, rep2)

#Remove loops that have an NA in any condition
merged.reps.nona <- na.omit(merged.reps)
merged.reps.nona <- merged.reps.nona %>%
  filter_all(all_vars(!is.infinite(.)))

################Averaging and getting fold changes##############
#These are fine, but they need to be averaged to hopefully control for some of the noisier data.
#Put in a list of cell line names that match columns of dataset
celllines = c()
conditions = c("UNT", "dTAG")
replicates = c("r1", "r2")
print(paste0(celllines, '_', conditions, '_', replicates))
condition.df <- crossing(celllines, conditions, replicates)
conditionnames <- unite(condition.df, "conditionnames", sep = "_")
samples <- crossing(celllines, conditions)
samplenames <- unite(samples, "samplenames", sep = "_")

#New df for storing stats with just loop info
strength.means <- merged.reps.nona[,c(1:6)]
#Convert tibble to vector
samplelist <- pull(samplenames, samplenames)
for (sample in samplelist) {
  #Get samples which match
  sampledata <- merged.reps.nona %>% select(contains(sample))
  outdata <- sampledata
  #Get geometric mean (since loop strengths are probably log norm distributed)
  outdata$mean <- rowMeans(log2(sampledata))
  outdata$mean <- 2^outdata$mean
  outdata$altmean <- outdata$mean - 1
  #Also get the variance
  outdata$sd <- apply(log2(sampledata),1,sd)
  outdata$sd <- 2^outdata$sd
  #Rename column to unique name based on sample
  outdata <- rename(outdata, !!paste0(sample, '_mean'):=mean)
  outdata <- rename(outdata, !!paste0(sample, '_sd'):=sd)
  outdata <- rename(outdata, !!paste0(sample, '_altmean'):=altmean)
  #Append to new df
  strength.means <- bind_cols(strength.means, outdata[,c(3,4,5)])
}

#New df for fold changes
strength.fcs <- merged.reps.nona[,c(1:6)]

#Next, calculate fold changes/frac changes
for (cellline in celllines) {
  #Get data which match from the mean df
  meandata.unt <- strength.means %>% select(contains(paste0(cellline, "_UNT_mean")))
  meandata.dtag <- strength.means %>% select(contains(paste0(cellline, "_dTAG_mean")))
  meandata <- bind_cols(meandata.unt, meandata.dtag)
  #Calculate FC
  meandata$FC <- meandata[,2]/meandata[,1]
  #Rename column
  meandata <- rename(meandata, !!paste0(cellline, '_meanFC'):=FC)
  strength.fcs <- bind_cols(strength.fcs, meandata %>% select(3))
}

#Merge the two together
strength.means.fcs <- join(strength.means, strength.fcs)

##############Adding strengths##############
#Looking at the data spreads to work out what to toss
#Combine all data together
strength.all <- join(strength.means.fcs, merged.reps.nona)

#Add a length column as well
strength.all$length <- ((strength.all$end2 + strength.all$start2)/2) - ((strength.all$end1 + strength.all$start1)/2)

##################Filtering based on sd################
#First, use an average, and filter out loops with a mean over 1.3 (using a stricter cutoff because I want to still restrict how many "often noisy" loops I get)
#Second, filter only based on UNT.
#Make a df only containing the UNT sd values
filter.table <- strength.all %>% rowwise()  %>% select(contains("UNT_sd"))
#Then, use this table to find the relevant lines in the list of all strengths to keep:
sd.cutoff.alt <- 1.3
strength.all.reduced.alt <- strength.all[which(rowMeans(filter.table) < sd.cutoff.alt),]
strength.all.noisy.alt <- strength.all[which(rowMeans(filter.table) >= sd.cutoff.alt),]
#Now this maintains 712 and tosses 335.
#Recombine for plotting
strength.all.labelled <- bind_rows("robust" = strength.all.reduced.alt, "noisy" = strength.all.noisy.alt, .id = 'set')
strength.all.labelled$set <- factor(strength.all.labelled$set, levels = c("robust", "noisy"))

#Looking at length:
ggplot(data = strength.all.labelled, aes(x = log10(length), colour = set))+
  theme_bw()+
  theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text.x = element_text(size=20, colour = 'black'), axis.text.y = element_text(size=20, colour = 'black'), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), plot.title = element_text(size=15,hjust = 0.5), legend.title=element_blank(),legend.position="right",legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  geom_histogram(alpha=0.15,size=1.75, pch=16, position = 'identity')+
  scale_y_continuous(limits=c(0,150),breaks = seq(0,140, by = 20)) +
  scale_x_continuous(limits=c(3,7),breaks = seq(3,7, by = 1))
#Uniformly very long loops that are filtered

################Filter based on loop length#################

#Short lengths are probably irrelevant. Remove loops where the end1 to start2 length is less than the length of either side
strength.all.reduced.noshort <- strength.all.reduced.alt %>% rowwise()  %>% filter(start2 - end1 > end1 - start1 & start2 - end1 > end2 - start2)
strength.all.reduced.short <- strength.all.reduced.alt %>% rowwise()  %>% filter(start2 - end1 < end1 - start1 | start2 - end1 < end2 - start2)

#Recombine for plotting
strength.all.labelled.short <- bind_rows(strength.all.reduced.noshort, strength.all.reduced.short, .id = 'short')

#Replot length hist
ggplot(data = strength.all.labelled.short, aes(x = log10(length), colour = short))+
  theme_bw()+
  theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text.x = element_text(size=20, colour = 'black'), axis.text.y = element_text(size=20, colour = 'black'), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), plot.title = element_text(size=15,hjust = 0.5), legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  geom_histogram(alpha=0.15,size=1.75, pch=16, position = 'identity')+
  # geom_vline(xintercept = 0, linetype = 2, color = 'black', size = 1) +
  # ylim(0, 5) + xlim(0, 5) +
  scale_y_continuous(limits=c(0,150),breaks = seq(0,140, by = 20)) +
  scale_x_continuous(limits=c(3,7),breaks = seq(3,7, by = 1))

###############Filter for UNT enrichment###################
#Also remove the ones that aren't really loops - i.e., they have minimal enrichment
#Look at UNT vs fold change to determine cutoff:
ggplot(data = strength.all.labelled.short, aes(x = .data[["dTAG.RING1AB_dTAG.SUZ12_UNT_mean"]], y = log2(.data[["dTAG.RING1AB_dTAG.SUZ12_meanFC"]])))+
  theme_bw()+
  theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), plot.title = element_text(size=15,hjust = 0.5), legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  geom_point(alpha=0.35,size=1.75, pch=16)+
  scale_y_continuous(limits=c(-3,3),breaks = seq(-3,3, by = 1))+
  scale_x_continuous(limits=c(-1,5), breaks = seq(-1,5, by = 1))+
  geom_hline(yintercept=0,linetype="dashed",lwd=2)
#Set to 1.1
enrichment.cutoff <- 1.1

filter.table <- strength.all.reduced.noshort %>% rowwise()  %>% select(contains("UNT_mean"))
strength.all.reduced.noshort.noweak <- strength.all.reduced.noshort[which(rowMeans(filter.table) > enrichment.cutoff),]
strength.all.reduced.noshort.weak <- strength.all.reduced.noshort[which(rowMeans(filter.table) <= enrichment.cutoff),]
#This culled ~150 loops to get to 517

################Anchor info###############
strength.all.noisy.alt #These are the ones with sd values too high
strength.all.reduced.short #These are the ones that are too short relative to their loop anchors
strength.all.reduced.noshort.weak #These are the weak loops
strength.all.reduced.noshort.noweak #These are the "good" ones

#Output these as bedpe files
write.table(strength.all.noisy.alt[,c(1:6)], file="C:/Users/newto/Hansen Lab Dropbox/Members/Miles/HTS Data/Polycomb/loops/PcG_RCMC_noisy.bedpe", quote=FALSE, sep='\t', row.names=FALSE)
write.table(strength.all.reduced.short[,c(1:6)], file="C:/Users/newto/Hansen Lab Dropbox/Members/Miles/HTS Data/Polycomb/loops/PcG_RCMC_short.bedpe", quote=FALSE, sep='\t', row.names=FALSE)
write.table(strength.all.reduced.noshort.weak[,c(1:6)], file="C:/Users/newto/Hansen Lab Dropbox/Members/Miles/HTS Data/Polycomb/loops/PcG_RCMC_weak.bedpe", quote=FALSE, sep='\t', row.names=FALSE)
write.table(strength.all.reduced.noshort.noweak[,c(1:6)], file="C:/Users/newto/Hansen Lab Dropbox/Members/Miles/HTS Data/Polycomb/loops/PcG_RCMC_good.bedpe", quote=FALSE, sep='\t', row.names=FALSE)

#Combine these together
strength.total <- bind_rows("noisy" = strength.all.noisy.alt, "short" = strength.all.reduced.short, "weak" = strength.all.reduced.noshort.weak, "good" = strength.all.reduced.noshort.noweak, .id = 'set')

write.table(strength.total[,c(1:7)], file="C:/Users/newto/Hansen Lab Dropbox/Members/Miles/HTS Data/Polycomb/loops/PcG_RCMC_all.bedpe", quote=FALSE, sep='\t', row.names=FALSE)





