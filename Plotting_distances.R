#Testing distance data from Sumin
library(readtext)
library(ggplot2)
library(ggsignif)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(ggbeeswarm)

outdir <- "/path/to/r_outputs"
if (!dir.exists(outdir)) {dir.create(outdir)}

#Read in data - comment out all but desired degron
depletion <- "PRC12"
#PRC1/2 pos
path <- "/path/to/Tbx5-Lhx5_PRC1_2_E1_Rep1_catdata.csv"
rep1 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-Lhx5_PRC1_2_E1_Rep2_catdata.csv"
rep2 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-Lhx5_PRC1_2_E1_Rep3_catdata.csv"
rep3 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
data <- bind_rows(rep1, rep2, rep3, .id = "rep")
pos <- bind_rows(rep1, rep2, rep3, .id = "rep")
#PRC1/2 neg
path <- "/path/to/Tbx5-NegCon_PRC1_2_A4_Rep1_catdata.csv"
rep1 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-NegCon_PRC1_2_A4_Rep2_catdata.csv"
rep2 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-NegCon_PRC1_2_A4_Rep3_catdata.csv"
rep3 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
data <- bind_rows(rep1, rep2, rep3, .id = "rep")
neg <- bind_rows(rep1, rep2, rep3, .id = "rep")

#PHC1/2/3 pos
depletion <- "PHC123"
path <- "/path/to/Tbx5-Lhx5_PHC1_2_3_D2_Rep1_catdata.csv"
rep1 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-Lhx5_PHC1_2_3_D2_Rep2_catdata.csv"
rep2 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-Lhx5_PHC1_2_3_D2_Rep3_catdata.csv"
rep3 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
data <- bind_rows(rep1, rep2, rep3, .id = "rep")
pos <- bind_rows(rep1, rep2, rep3, .id = "rep")
#PHC1/2/3 neg
path <- "/path/to/Tbx5-NegCon_PHC1_2_3_E10_Rep1_catdata.csv"
rep1 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-NegCon_PHC1_2_3_E10_Rep2_catdata.csv"
rep2 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-NegCon_PHC1_2_3_E10_Rep3_catdata.csv"
rep3 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
data <- bind_rows(rep1, rep2, rep3, .id = "rep")
neg <- bind_rows(rep1, rep2, rep3, .id = "rep")

#Read in data
depletion <- "PCGF24"
path <- "/path/to/Tbx5-Lhx5_PCGF2_4_B6_Rep1_catdata.csv"
rep1 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-Lhx5_PCGF2_4_B6_Rep2_catdata.csv"
rep2 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-Lhx5_PCGF2_4_B6_Rep3_catdata.csv"
rep3 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
data <- bind_rows(rep1, rep2, rep3, .id = "rep")
pos <- bind_rows(rep1, rep2, rep3, .id = "rep")

path <- "/path/to/Tbx5-NegCon_PCGF2_4_H10_Rep1_catdata.csv"
rep1 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-NegCon_PCGF2_4_H10_Rep2_catdata.csv"
rep2 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
path <- "/path/to/Tbx5-NegCon_PCGF2_4_H10_Rep3_catdata.csv"
rep3 <- read.delim(path, header = TRUE, sep = ",", dec = ".")
data <- bind_rows(rep1, rep2, rep3, .id = "rep")
neg <- bind_rows(rep1, rep2, rep3, .id = "rep")


line_data <- bind_rows(pos, neg, .id = "line")
line_data$line <- as.factor(line_data$line)

replist.new <- c('1', '2', '3')

stackmeans <- line_data %>% group_by(condition, stack, line) %>% summarise(distance = mean(distance..nm.), rep = unique(rep))
repmeans <- line_data %>% group_by(condition, rep, line) %>% summarise(distance = mean(distance..nm.))

stackmeds <- line_data %>% group_by(condition, stack, line) %>% summarise(distance = median(distance..nm.), rep = unique(rep))
repmeds <- line_data %>% group_by(condition, rep, line) %>% summarise(distance = median(distance..nm.))
medians <- line_data %>% group_by(condition, line) %>% summarise(distance = median(distance..nm.))
means <- line_data %>% group_by(condition, line) %>% summarise(distance = mean(distance..nm.))
counts.perrep <- line_data %>% group_by(condition, rep, line) %>% summarise(distance = length(distance..nm.))
counts <- line_data %>% group_by(condition, line) %>% summarise(distance = length(distance..nm.))

#Comment out as appropriate



#Plot CDF
#Tbx5Lhx5:
cellline <- "Tbx5Lhx5"
pdf(paste(outdir, "/", depletion, "_", cellline, "_CDF", ".pdf", sep = ""))

p = ggplot(data = line_data[which(line_data$line == 1),], aes(distance..nm., colour = condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_ecdf() +
  # geom_signif(comparisons = list(c(1, 2)),
  #             step_increase = 0.1,
  #             test = "ks.test") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(limits=c(0,800),breaks = seq(0,800, by = 100))
print(p)
dev.off()
#Tbx5negcon:
cellline <- "Tbx5negcon"
pdf(paste(outdir, "/", depletion, "_", cellline, "_CDF", ".pdf", sep = ""))

p = ggplot(data = line_data[which(line_data$line == 2),], aes(distance..nm., colour = condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_ecdf() +
  # geom_signif(comparisons = list(c(1, 2)),
  #             step_increase = 0.1,
  #             test = "ks.test") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(limits=c(0,800),breaks = seq(0,800, by = 100))
print(p)
dev.off()

#Plot histogram
pdf(paste(outdir, "/", depletion, "_histogram", ".pdf", sep = ""))
p <- ggplot(data = line_data, aes(distance..nm., colour = interaction(line, condition, sep = ":"), after_stat(density))) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_histogram(alpha = 0.5, position = 'identity', fill = NA, binwidth = 40, breaks = seq(0, 800, by = 40)) +
  coord_cartesian(ylim = c(0, 0.005)) +
  scale_x_continuous(limits=c(0,850),breaks = seq(0,800, by = 100))
print(p)
dev.off()





#Downsample the data
nrow(rep1[which(rep1$condition == "dTAG"),])

ggplot(data = rep2, aes(distance..nm., colour = condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_ecdf() +
  # geom_signif(comparisons = list(c(1, 2)),
  #             step_increase = 0.1,
  #             test = "ks.test") +
  coord_cartesian(ylim = c(0, 1))

rep2.sample <- rep2[sample(nrow(rep2), 300),]
ggplot(data = rep2.sample, aes(distance..nm., colour = condition)) +
  theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
  theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_ecdf() +
  # geom_signif(comparisons = list(c(1, 2)),
  #             step_increase = 0.1,
  #             test = "ks.test") +
  coord_cartesian(ylim = c(0, 1))








