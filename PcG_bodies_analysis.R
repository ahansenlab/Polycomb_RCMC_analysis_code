#Code for analysing Polycomb body data output from the segmentation code

library(readtext)
library(ggplot2)
library(ggsignif)
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
library(ggbeeswarm)

path <- "/path/to/output/files"

conditions <- c("UNT", "dTAG")


#Put cell line names here:
celllines <- c("CBX2",
               "CBX7",
               "PCGF24",
               "PCGF24SUZ12",
               "PCGF124",
               "PCGF246",
               "PHC123",
               "PHC123CBX27",
               "SUZ12"
)

for (cellline in celllines) {
  bodies <- data.frame()
  bodies.corrected <- data.frame()
  for (condition in conditions) {
    data <- read.delim(paste0(path, list.files(path = path, pattern = paste0("^", cellline, "_", condition))), header = TRUE, sep = ",", dec = ".")
    data$condition <- condition
    #Need an index for all nuclei in one condition
    data$condition_nuc_id <- interaction(data$nucleus_id, data$folder)
    for (i in 1:length(unique(data$condition_nuc_id))) {
      data$condition_nuc_index[data$condition_nuc_id == unique(data$condition_nuc_id)[[i]]] <- i
    }
    #Also produce a column that indicates the replicate
    data$replicate <- str_split_i(data$folder, "_", 1)
    bodies <- rbind(bodies, data)
  }
  
  #List of replicates:
  replist <- unique(bodies$replicate)
  replist.new <- c('rep1', 'rep2', 'rep3')
  #Rename the replicates using this
  bodies$replicate <- mapvalues(bodies$replicate, from = replist, to = replist.new)
  
  
  #Extract the nuclei data
  nuclei.dups <- data.frame(file = bodies$folder, nucleus_id = bodies$condition_nuc_index, nuc_mint = bodies$nuc_mean_intensity..au., nuc_vol = bodies$nuc_volume..px., condition = bodies$condition, replicate = bodies$replicate)
  #Remove all the duplicate entries
  nuclei <- unique(nuclei.dups)
  
  # calculate total intensity for each body by multiplying volume by mean intensity
  bodies$TotalInt <- bodies$mean_intensity..au. * bodies$body_volume..px.
  #QC here to check for distribution of Polycomb body sizes to decide how to filter:
  # ggplot(data = bodies, aes(x = log2(body_volume..px.)))+
  #   theme_bw()+
  #   theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), plot.title = element_text(size=15,hjust = 0.5), legend.title=element_blank(),legend.position=c(0.85,0.2),legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank())+
  #   geom_histogram(binwidth = 0.1)+
  #   scale_x_continuous(limits=c(0,10), breaks = seq(0,10, by = 1))
  #Based on this, filtering between 2^4 and 2^5 seems reasonable. 
  #Will use 27 px as the cutoff.
  
  # filter for bodies smaller than 27 px - this is 3x3x3 with isotropic pixels 
  bodies.corrected <- bodies[ which(bodies$body_volume..px. > 27), ]
  
  if (cellline == "SUZ12") {
    filternuclei <- c(92, 93, 94, 95)
    bodies.corrected <- bodies.corrected |> filter_out(condition_nuc_index %in% filternuclei & condition == "UNT")
    nuclei <- nuclei |> filter_out(nucleus_id %in% filternuclei & condition == "UNT")
  } else {filternuclei <- c()}
  
  nuclei$TotalInt <- nuclei$nuc_mint * nuclei$nuc_vol
  
  # count how PcG bodies per cell in each condition
  for (i in 1:max(nuclei$nucleus_id)) {
    unt.count <- length(which(bodies.corrected$condition == "UNT" & bodies.corrected$condition_nuc_index == i))
    dtag.count <- length(which(bodies.corrected$condition == "dTAG" & bodies.corrected$condition_nuc_index == i))
    nuclei$count[nuclei$condition == "UNT" & nuclei$nucleus_id == i] <- unt.count
    nuclei$count[nuclei$condition == "dTAG" & nuclei$nucleus_id == i] <- dtag.count
  }
  
  bodies.corrected.fraction <- data.frame(IntFraction = vector(mode = "numeric", length = length(nuclei$nucleus_id)), VolFraction = vector(mode = "numeric", length = length(nuclei$nucleus_id)), Enrichment = vector(mode = "numeric", length = length(nuclei$nucleus_id)), Condition = nuclei$condition, CellNumber = nuclei$nucleus_id)
  
  bodies.corrected.cellmeans <- data.frame(NucleoMeanInt = vector(mode = "numeric", length = length(nuclei$nucleus_id)), MeanBodyMeanInt = vector(mode = "numeric", length = length(nuclei$nucleus_id)), Enrichment = vector(mode = "numeric", length = length(nuclei$nucleus_id)), MeanBodyVol = vector(mode = "numeric", length = length(nuclei$nucleus_id)), MeanBodyTotalInt = vector(mode = "numeric", length = length(nuclei$nucleus_id)), Condition = nuclei$condition, CellNumber = nuclei$nucleus_id)
  
  cell.count <- 1
  
  for (condition in conditions) {
    all.condition.bodies.corrected <- bodies.corrected[ which(bodies.corrected$condition == condition), ] # filter for only ones with same condition
    all.condition.nucleus <- nuclei[ which(nuclei$condition == condition), ] # same for nuclei
    for (m in 1:max(all.condition.nucleus$nucleus_id)) {
      if (cellline == "SUZ12" & m %in% filternuclei & condition == "UNT") {
        next
      }
      cell.bodies.corrected <- all.condition.bodies.corrected[ which(all.condition.bodies.corrected$condition_nuc_index == m), ]  # take all bodies with matching cell number
      cell.nucleus <- all.condition.nucleus[ which(all.condition.nucleus$nucleus_id == m), ] # same for nucleus
      total.intensity.body.corrected <- sum(cell.bodies.corrected$TotalInt)
      cell.body.fraction.int.corrected <- total.intensity.body.corrected / cell.nucleus$TotalInt
      bodies.corrected.fraction$IntFraction[cell.count] <- cell.body.fraction.int.corrected
      
      total.vol.body.corrected <- sum(cell.bodies.corrected$body_volume..px.)
      cell.body.fraction.vol.corrected <- total.vol.body.corrected / cell.nucleus$nuc_vol
      bodies.corrected.fraction$VolFraction[cell.count] <- cell.body.fraction.vol.corrected
      
      mean.meanintensity.body.corrected <- mean(cell.bodies.corrected$mean_intensity..au.)
      mean.totalintensity.body.corrected <- mean(cell.bodies.corrected$TotalInt)
      mean.volume.body.corrected <- mean(cell.bodies.corrected$body_volume..px.)
      bodies.corrected.cellmeans$MeanBodyMeanInt[cell.count] <- mean.meanintensity.body.corrected
      bodies.corrected.cellmeans$MeanBodyVol[cell.count] <- mean.volume.body.corrected
      bodies.corrected.cellmeans$MeanBodyTotalInt[cell.count] <- mean.totalintensity.body.corrected
      
      nucleo.int.corrected <- cell.nucleus$TotalInt - total.intensity.body.corrected
      nucleo.vol.corrected <- cell.nucleus$nuc_vol - total.vol.body.corrected
      nucleo.mean.int.corrected <- nucleo.int.corrected / nucleo.vol.corrected
      body.mean.int.corrected <- total.intensity.body.corrected / total.vol.body.corrected
      bodies.corrected.fraction$Enrichment[cell.count] <- body.mean.int.corrected / nucleo.mean.int.corrected
      bodies.corrected.fraction$replicate[cell.count] <- cell.nucleus$replicate
      
      
      bodies.corrected.cellmeans$NucleoMeanInt[cell.count] <- nucleo.mean.int.corrected
      bodies.corrected.cellmeans$Enrichment[cell.count] <- mean.meanintensity.body.corrected / nucleo.mean.int.corrected
      bodies.corrected.cellmeans$replicate[cell.count] <- cell.nucleus$replicate
      
      cell.count <- cell.count + 1
    }
  }
  
  
  ################ PREPARING THE DFS #############
  
  bodies.corrected$condition <- factor(bodies.corrected$condition, levels = c("UNT", "dTAG"))
  nuclei$condition <- factor(nuclei$condition, levels = c("UNT", "dTAG"))
  bodies.corrected.fraction$Condition <- factor(bodies.corrected.fraction$Condition, levels = c("UNT", "dTAG"))
  bodies.corrected.cellmeans$Condition <- factor(bodies.corrected.cellmeans$Condition, levels = c("UNT", "dTAG"))

  ############################################# Enrichment calc #####################################
  plot.bodies.cellmeans <- melt(bodies.corrected.cellmeans, id.vars = c("Condition", "replicate"), measure.vars = c("NucleoMeanInt", "MeanBodyMeanInt"))
  plot.bodies.cellmeans$variable <- mapvalues(plot.bodies.cellmeans$variable, from = c("NucleoMeanInt", "MeanBodyMeanInt"), to = c("Nucleoplasm", "Polycomb Body"))
  plot.bodies.cellmeans$variable <- factor(plot.bodies.cellmeans$variable, levels = c("Polycomb Body", "Nucleoplasm")) # order so bodies first
  plot.bodies.cellmeans$Condition <- factor(plot.bodies.cellmeans$Condition, levels = c("UNT", "dTAG"))
  
  UNT.cellmeans <- bodies.corrected.cellmeans[which(bodies.corrected.cellmeans$Condition == "UNT"), ]
  UNT.enrichment <- mean(UNT.cellmeans$Enrichment)
  dTAG.cellmeans <- bodies.corrected.cellmeans[which(bodies.corrected.cellmeans$Condition == "dTAG"), ]
  dTAG.enrichment <- mean(dTAG.cellmeans$Enrichment)
  
  UNT.median.nucleo <- median(UNT.cellmeans$NucleoMeanInt) # measure median for non-body for each condition (so that median line is aligned with 1)
  dTAG.median.nucleo <- median(dTAG.cellmeans$NucleoMeanInt)
  norm.plot.bodies.cellmeans <- plot.bodies.cellmeans
  norm.plot.bodies.cellmeans$value <- ifelse(norm.plot.bodies.cellmeans$Condition == "UNT", norm.plot.bodies.cellmeans$value / UNT.median.nucleo, norm.plot.bodies.cellmeans$value) # normalise each condition's mean Polycomb body signal density with the corresponding median
  norm.plot.bodies.cellmeans$value <- ifelse(norm.plot.bodies.cellmeans$Condition == "dTAG", norm.plot.bodies.cellmeans$value / dTAG.median.nucleo, norm.plot.bodies.cellmeans$value)
  
  ############################################# Per replicate means #####################################
  
  #Get the per-replicate means
  repsmeans.bodies <- bodies.corrected %>% group_by(condition, replicate) %>% summarise(body_vol = mean(body_volume..px.), body_mint = mean(mean_intensity..au.), body_tint = mean(TotalInt))
  repsmeans_nuclei <- nuclei %>% group_by(condition, replicate) %>% summarise(nuc_vol = mean(nuc_vol), nuc_tint = mean(TotalInt), percell = mean(count), vol_Q1 = mean(Q1ct), vol_Q2 = mean(Q2ct), vol_Q3 = mean(Q3ct), vol_Q4 = mean(Q4ct), mint_Q1 = mean(Q1ctmint), mint_Q2 = mean(Q2ctmint), mint_Q3 = mean(Q3ctmint), mint_Q4 = mean(Q4ctmint), nuc_mint = mean(nuc_mint))
  repmeans_percell_measures <- bodies.corrected.cellmeans %>% group_by(Condition, replicate) %>% summarise(body_mint = mean(MeanBodyMeanInt, na.rm = TRUE), body_vol = mean(MeanBodyVol, na.rm = TRUE))
  cellmeans.means <- norm.plot.bodies.cellmeans %>% group_by(Condition, replicate, variable) %>% summarise(nucleo_normalised_signal = mean(value, na.rm = TRUE))
  fractions.means <- bodies.corrected.fraction %>% group_by(Condition, replicate) %>% summarise(IntFraction = mean(IntFraction, na.rm = TRUE), VolFraction = mean(VolFraction, na.rm = TRUE))
  
  ############################################## Overall means #######################################
  outdir <- "C:/Users/newto/OneDrive/Documents/Work/Data/Hansen Lab/Polycomb bodies/means/"
  if (!dir.exists(outdir)) {dir.create(outdir)}
  #Also get the means for each condition
  means_percell_measures <- bodies.corrected.cellmeans %>% group_by(Condition) %>% summarise(body_vol = mean(MeanBodyVol, na.rm = TRUE), count = length(MeanBodyVol))
  means_nuclei <- nuclei %>% group_by(condition) %>% summarise(percell = mean(count), nuc_mint = mean(nuc_mint))
  
  write.table(means_percell_measures, paste0(outdir, cellline , '_volmeans.txt'), quote=FALSE, col.names=TRUE)
  write.table(means_nuclei, paste0(outdir, cellline , '_countmeans.txt'), quote=FALSE, col.names=TRUE)
  
  
  
  ######################################### PLOTS #####################################################
  outdir <- "C:/Users/newto/OneDrive/Documents/Work/Data/Hansen Lab/Polycomb bodies/plots/"
  ############ Nuclear volumes
  pdf(paste(outdir, cellline,  "_nuclearvol", ".pdf", sep = ""))

  p <- ggplot(data = nuclei, aes(x = condition, y = nuc_vol * (0.0558809^3), fill = condition)) +
    theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
    theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    # geom_violin(lwd=0.5, colour="black") +
    geom_beeswarm(alpha = 0.5, cex = 1.5, aes(x = condition, y = nuc_vol * (0.0558809^3), colour = replicate)) +
    geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA, width = 0.05) +
    geom_beeswarm(data = repsmeans_nuclei, size = 3, cex = 3, colour = "black", stroke = 1, pch = 21, aes(y = nuc_vol * (0.0558809^3), fill = replicate)) +
    scale_fill_manual(values=c("#66ccee", "#ee6677", "#725def", "#dd217d", "#ffb00d"), breaks=c("UNT", "dTAG", replist.new)) +
    scale_colour_manual(values=c("#725def", "#dd217d", "#ffb00d"), breaks=c(replist.new)) +
    geom_signif(comparisons = list(c(1, 2)),
                step_increase = 0.1,
                test = "wilcox.test") +
    expand_limits(y = 0)
  print(p)
  dev.off()
  ############### Nucleus mean intensity
  pdf(paste(outdir, cellline,  "_meannuclearfluor", ".pdf", sep = ""))

  p <- ggplot(data = nuclei, aes(x = condition, y = nuc_mint)) +  # nuclei total intensity boxplot
    theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
    theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_beeswarm(stroke = NA, alpha = 0.5, cex = 1, aes(colour = replicate)) +
    geom_boxplot(fill = NA, lwd=0.5, fatten = 2, colour="black", outlier.size=NA, outlier.colour=NA, width = 0.05) +
    geom_beeswarm(data = repsmeans_nuclei, size = 3, cex = 3, colour = "black", stroke = 1, pch = 21, aes(y = nuc_mint, fill = replicate)) +
    scale_fill_manual(values=c("#66ccee", "#ee6677", "#725def", "#dd217d", "#ffb00d"), breaks=c("UNT", "dTAG", replist.new)) +
    scale_colour_manual(values=c("#725def", "#dd217d", "#ffb00d"), breaks=c(replist.new)) +
    geom_signif(comparisons = list(c(1, 2)),
                step_increase = 0.1,
                test = "wilcox.test") +
    coord_cartesian(ylim = c(0, 5000))
  print(p)
  dev.off()
  
  ############## bodies per cell
  pdf(paste(outdir, cellline,  "_bodycount", ".pdf", sep = ""))

  p <- ggplot(data = nuclei, aes(x = condition, y = count)) +
    theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
    theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_beeswarm(alpha = 0.5, cex = 1, stroke = NA, aes(colour = replicate)) +
    geom_boxplot(fill = NA, lwd=0.5, fatten = 2, colour="black", outlier.size=NA, outlier.colour=NA, width = 0.05) +
    geom_beeswarm(data = repsmeans_nuclei, size = 3, cex = 3, colour = "black", stroke = 1, pch = 21, aes(y = percell, fill = replicate)) +
    scale_fill_manual(values=c("#66ccee", "#ee6677", "#725def", "#dd217d", "#ffb00d"), breaks=c("UNT", "dTAG", replist.new)) +
    scale_colour_manual(values=c("#725def", "#dd217d", "#ffb00d"), breaks=c(replist.new)) +
    geom_signif(comparisons = list(c(1, 2)),
                step_increase = 0.1,
                test = "wilcox.test", y_position = 190) +
    coord_cartesian(ylim = c(0, 200))
  print(p)
  dev.off()

  ##############  cell mean body mean int
  pdf(paste(outdir, cellline,  "_cellmeanbodymint", ".pdf", sep = ""))

  p <- ggplot(data = bodies.corrected.cellmeans, aes(x = Condition, y = MeanBodyMeanInt)) +
    theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
    theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_beeswarm(stroke = NA, alpha = 0.5, cex = 1, aes(colour = replicate)) +
    geom_boxplot(fill = NA, lwd=0.5, fatten = 2, colour="black", outlier.size=NA, outlier.colour=NA, width = 0.05) +
    geom_beeswarm(data = repmeans_percell_measures, size = 3, cex = 3, colour = "black", stroke = 1, pch = 21, aes(y = body_mint, fill = replicate)) +
    scale_fill_manual(values=c("#66ccee", "#ee6677", "#725def", "#dd217d", "#ffb00d"), breaks=c("UNT", "dTAG", replist.new)) +
    scale_colour_manual(values=c("#725def", "#dd217d", "#ffb00d"), breaks=c(replist.new)) +
    geom_signif(comparisons = list(c(1, 2)),
                step_increase = 0.1,
                test = "wilcox.test", y_position = 10000) +
    coord_cartesian(ylim = c(0, 11000))
  print(p)
  dev.off()
  ##############  cell mean body vol
  pdf(paste(outdir, cellline,  "_cellmeanbodyvol", ".pdf", sep = ""))

  p <- ggplot(data = bodies.corrected.cellmeans, aes(x = Condition, y = MeanBodyVol * (0.05509^3))) +
    theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
    theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_beeswarm(stroke = NA, alpha = 0.5, cex = 0.75, aes(colour = replicate)) +
    geom_boxplot(fill = NA, lwd=0.5, fatten = 2, colour="black", outlier.size=NA, outlier.colour=NA, width = 0.05) +
    geom_beeswarm(data = repmeans_percell_measures, size = 3, cex = 3, colour = "black", stroke = 1, pch = 21, aes(y = body_vol * (0.05509^3), fill = replicate)) +
    scale_fill_manual(values=c("#66ccee", "#ee6677", "#725def", "#dd217d", "#ffb00d"), breaks=c("UNT", "dTAG", replist.new)) +
    scale_colour_manual(values=c("#725def", "#dd217d", "#ffb00d"), breaks=c(replist.new)) +
    geom_signif(comparisons = list(c(1, 2)),
                step_increase = 0.1,
                test = "wilcox.test", y_position = 0.055) +
    coord_cartesian(ylim = c(0, 0.06))
  print(p)
  dev.off()
  ############## nucleo mean int vs body mean int
  pdf(paste(outdir, cellline,  "_bodyvnucleo", ".pdf", sep = ""))

  p <- ggplot(data = norm.plot.bodies.cellmeans[which(norm.plot.bodies.cellmeans$variable == "Polycomb Body"),], aes(x = Condition, y = value)) +
    theme_bw()+ theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
    theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_beeswarm(stroke = NA, alpha = 0.5, cex = 1, aes(colour = replicate)) +
    geom_boxplot(fill = NA, lwd=0.5, fatten = 2, colour="black", outlier.size=NA, outlier.colour=NA, width = 0.05) +
    geom_beeswarm(data = cellmeans.means[which(cellmeans.means$variable == "Polycomb Body"),], size = 3, cex = 3, colour = "black", stroke = 1, pch = 21, aes(y = nucleo_normalised_signal, fill = replicate)) +
    scale_fill_manual(values=c("#66ccee", "#ee6677", "#725def", "#dd217d", "#ffb00d"), breaks=c("UNT", "dTAG", replist.new)) +
    scale_colour_manual(values=c("#725def", "#dd217d", "#ffb00d"), breaks=c(replist.new)) +
    coord_cartesian(ylim = c(0, 4.5))
  print(p)
  dev.off()
  
  #Extra plot to get significance values
  pdf(paste(outdir, cellline,  "_bodyvnucleowrap", ".pdf", sep = ""))
  # plot individual panels - stats
  p <- ggplot(data = norm.plot.bodies.cellmeans, aes(x = Condition, y = value, fill = Condition)) +
    theme_bw()+
    theme(aspect.ratio=1, axis.title.x=element_blank(), legend.text=element_text(size=10, colour="black"), axis.ticks = element_line(colour = "black", size = 0.75), axis.text.x = element_text(size=10,colour="black"), axis.text.y = element_text(size=10, colour="black"), plot.title=element_text(size=15, colour="black", hjust=0.5, face = 'bold'), axis.title.y=element_text(size=10, colour="black"))+
    theme(panel.border=element_rect(colour="black",size=1), panel.margin=unit(1, "lines"), strip.text = element_text(size=10, colour="black", face = "bold"), strip.background = element_blank(), legend.title=element_blank(),legend.position='right',legend.key = element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_boxplot(lwd=0.5, fatten = 0.8, colour="black", outlier.size=NA, outlier.colour=NA) + scale_fill_manual(values=c("blue", "green", "purple"), breaks=c("UNT", "dTAG", "dTAGTAM")) +
    ylab("Bodies per cell") +
    expand_limits(y = 0) +
    facet_wrap(vars(variable), scales = "free_x") +
    geom_signif(comparisons = list(c(1, 2)),
                step_increase = 0.1,
                test = "wilcox.test",
                test.args = list(exact = TRUE))
  print(p)
  dev.off()
  ##########################################
}
