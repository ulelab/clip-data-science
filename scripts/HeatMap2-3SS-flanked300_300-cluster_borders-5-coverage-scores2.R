library(gplots)
library("ggplot2")
library("smoother")
library("cowplot")

args<-commandArgs(TRUE)

###################
# regulated exons #
###################

# import and merge of 3SS with 5SS
data3SSinput <- read.table(args[1], header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
data5SSinput <- read.table(args[2], header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
#data3SSinput <- read.table("/Users/nhaberma/UCL/2017.04.28-CLIPo-iONMF-ENCODE/iCount2-For_Nejc/RNA-maps/PTB.enhanced.cassetteexons.hg38.bed-3SS-flanked300-clusters-map_positions-uniq.csv", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
#data5SSinput <- read.table("/Users/nhaberma/UCL/2017.04.28-CLIPo-iONMF-ENCODE/iCount2-For_Nejc/RNA-maps/PTB.enhanced.cassetteexons.hg38.bed-5SS-flanked300-clusters-map_positions-uniq.csv", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
data3SS.5SS.input <- merge(data3SSinput, data5SSinput, by="V1", all.x=TRUE)
data3SS.5SS.input <- data3SS.5SS.input[order(data3SS.5SS.input$V2.x),]

num.reg.exons <- nrow(data3SS.5SS.input)

# removal of transcripts with clusters outside of heat map region
data3SS.5SS <- data3SS.5SS.input[complete.cases(data3SS.5SS.input),]
data3SS.5SS$max3SS <- ""
data3SS.5SS$max5SS <- ""
for (i in 1:nrow(data3SS.5SS)) {  #we want to find the maximum value on ploted HeatMap so we can remove those which are out of our heatmap range
  data3SS.5SS$max3SS[i] <- max(data3SS.5SS[i,4:354])
  data3SS.5SS$max5SS[i] <- max(data3SS.5SS[i,857:1207])
}
data3SS.5SS <- data3SS.5SS[which(data3SS.5SS$max3SS >= 2 | data3SS.5SS$max5SS >= 2),]

# calculate cluster coverage
data3SS.5SS.input <- data3SS.5SS.input
data3SS.5SS.input[data3SS.5SS.input < 2] <- 0
data3SS.5SS.input[data3SS.5SS.input >= 2] <- 1
sum.data3SSinput <- colSums(data3SS.5SS.input[,4:354], na.rm = FALSE, dims = 1)
sum.data5SSinput <- colSums(data3SS.5SS.input[,857:1207], na.rm = FALSE, dims = 1)
sum.3SSA.norm <- sum.data3SSinput / num.reg.exons 
sum.5SSA.norm <- sum.data5SSinput / num.reg.exons 
sum.3SSA.norm.smooth <- smth(sum.3SSA.norm, window = 25, method = "gaussian") #SMOOTHING
sum.5SSA.norm.smooth <- smth(sum.5SSA.norm, window = 25, method = "gaussian") #SMOOTHING
ymax.regulated <- max(c(sum.5SSA.norm.smooth[!is.na(sum.5SSA.norm.smooth)]), max(sum.3SSA.norm.smooth[!is.na(sum.3SSA.norm.smooth)])) #find the y-max peak for any of the 

#calculate average coverage to add it as a score (only from intronic region)
norm.3SSA.coverage <- sum(sum.data3SSinput[1:300]) / num.reg.exons
norm.5SSA.coverage <- sum(sum.data5SSinput[50:301]) / num.reg.exons
norm.exon.coverage <- (sum(sum.data3SSinput[300:350]) + sum(sum.data5SSinput[1:51])) / num.reg.exons

#################
# control exons #
#################

# import and merge of 3SS with 5SS
control3SSinput <- read.table(args[3], header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
control5SSinput <- read.table(args[4], header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
#control3SSinput <- read.table("HepG2_PTBP1-dPSI_0.1-unchanged-sorted-merged.bed-3SS-flanked300-clusters-map_positions-merged.csv", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
#control5SSinput <- read.table("HepG2_PTBP1-dPSI_0.1-unchanged-sorted-merged.bed-5SS-flanked300-clusters-map_positions-merged.csv", header=FALSE, sep=",",stringsAsFactors=F,na.strings=c("NA", "-", "?"))
control3SS.5SSA <- merge(control3SSinput, control5SSinput, by="V1", all.x=TRUE)
control3SS.5SSA <- control3SS.5SSA[order(control3SS.5SSA$V2.x),]

num.reg.control.exons <- nrow(control3SS.5SSA)

# removal of transcripts with clusters outside of heat map region
control3SS.5SS <- control3SS.5SSA[complete.cases(control3SS.5SSA),]
control3SS.5SS$max3SS <- ""
control3SS.5SS$max5SS <- ""
for (i in 1:nrow(control3SS.5SS)) {  #we want to find the maximum value on ploted HeatMap so we can remove those which are out of our heatmap range
  control3SS.5SS$max3SS[i] <- max(control3SS.5SS[i,4:354])
  control3SS.5SS$max5SS[i] <- max(control3SS.5SS[i,857:1207])
}
control3SS.5SS <- control3SS.5SS[which(control3SS.5SS$max3SS >= 2 | control3SS.5SS$max5SS >= 2),]

# calculate cluster coverage
control.coverage3SS.5SS <- control3SS.5SS
control.coverage3SS.5SS[control.coverage3SS.5SS < 2] <- 0
control.coverage3SS.5SS[control.coverage3SS.5SS >= 2] <- 1
sum.control3SSinput <- colSums(control.coverage3SS.5SS[,4:354], na.rm = FALSE, dims = 1)
sum.control5SSinput <- colSums(control.coverage3SS.5SS[,857:1207], na.rm = FALSE, dims = 1)
sum.control.3SSA.norm <- sum.control3SSinput / num.reg.control.exons
sum.control.5SSA.norm <- sum.control5SSinput / num.reg.control.exons
sum.control.3SSA.norm.smooth <- smth(sum.control.3SSA.norm, window = 25, method = "gaussian") #SMOOTHING
sum.control.5SSA.norm.smooth <- smth(sum.control.5SSA.norm, window = 25, method = "gaussian") #SMOOTHING
ymax.control <- max(c(sum.control.5SSA.norm.smooth[!is.na(sum.control.5SSA.norm.smooth)]), max(sum.control.3SSA.norm.smooth[!is.na(sum.control.3SSA.norm.smooth)])) #find the y-max peak for any of the 

#calculate average coverage to add it as a score (only from intronic region)
norm.control.3SSA.coverage <- sum(sum.control3SSinput[1:300]) / num.reg.control.exons
norm.control.5SSA.coverage <- sum(sum.control5SSinput[50:301]) / num.reg.control.exons
norm.control.exon.coverage <- (sum(sum.control3SSinput[300:350]) + sum(sum.control5SSinput[1:51])) / num.reg.control.exons

ymax <- max(ymax.regulated, ymax.control)

# ggplot cluster density
library("ggplot2")
gg.3SS <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_area(aes(c(-300:50), sum.control.3SSA.norm.smooth), color="grey", alpha=0.3) + 
  geom_line(aes(c(-300:50), sum.3SSA.norm.smooth)) + 
  geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
  ggtitle("Cluster coverage") + 
  xlab("position relative to 3'SS") + 
  ylab("normalized density of XL clusters") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-300, 50)) +
  scale_y_continuous(limits = c(0, ymax))

gg.5SS <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_area(aes(c(-50:300), sum.control.5SSA.norm.smooth), color="grey", alpha=0.3) + 
  geom_line(aes(c(-50:300), sum.5SSA.norm.smooth)) + 
  geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
  ggtitle("Cluster coverage") + 
  xlab("position relative to 5'SS") + 
  ylab("normalized density of XL clusters") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_x_continuous(limits = c(-50, 300)) +
  scale_y_continuous(limits = c(0, ymax))

# separate 3SS from 5SS and convert it to matrix
rnames <- data3SS.5SS[,2]
data3SS <- data3SS.5SS[,4:354]
data3SS[,300] <- -1  #exon start
data5SS <- data3SS.5SS[,857:1206]
data5SS[,50] <- -1  #exon start

data3SS <- as.data.frame(data3SS)
data_matrix <- data.matrix(data3SS)

pdf(paste(args[5], sep=""), width = 8.50, height = 11)

h1 <- heatmap.2(
  data_matrix+1,
  dendrogram = "none",
  scale      = "none",
  trace      = "none",
  Rowv = FALSE,
  Colv = FALSE,
  key        = FALSE,
  labRow     = FALSE,
  labCol     = c(-300:50),
  #labCol = FALSE,
  col    = c("#FFFFFF", "#FFBDC0", "#F27F89", "#A12931", "#270004"),
  main="",
  key.xlab="pentamer denstiy",
  key.ylab="",
  key.title="",
  cexRow=0.2,
  cexCol=0.2
)

gg.3SS

# 5SS HeatMap
data5SS <- as.data.frame(data5SS)
data_matrix <- data.matrix(data5SS)

heatmap.2(
  data_matrix+1,
  dendrogram = "none",
  scale      = "none",
  trace      = "none",
  Rowv = FALSE,
  Colv = FALSE,
  key        = FALSE,
  labRow     = FALSE,
  labCol     = c(-50:300),
  #labCol = FALSE,
  col    = c("#FFFFFF", "#FFBDC0", "#F27F89", "#A12931", "#270004"),
  main="",
  key.xlab="pentamer denstiy",
  key.ylab="",
  key.title="",
  cexRow=0.2,
  cexCol=0.2
)

gg.5SS

# distance between regulated and controls
dist3SS.reg.vs.control <- norm.3SSA.coverage - norm.control.3SSA.coverage
dist5SS.reg.vs.control <- norm.5SSA.coverage - norm.control.5SSA.coverage
dist.exon.reg.vs.control <- norm.exon.coverage - norm.control.exon.coverage

# enrichment between regulated and controls
enr3SS.reg.vs.control <- norm.3SSA.coverage / norm.control.3SSA.coverage
enr5SS.reg.vs.control <- norm.5SSA.coverage / norm.control.5SSA.coverage
enr.exon.reg.vs.control <- norm.exon.coverage / norm.control.exon.coverage

# Table 1
reg.exons.cov <- c(norm.3SSA.coverage, norm.exon.coverage, norm.5SSA.coverage)
control.exons.cov <- c(norm.control.3SSA.coverage, norm.control.exon.coverage, norm.control.5SSA.coverage)
enrichment <- c(norm.3SSA.coverage/norm.control.3SSA.coverage, norm.exon.coverage/norm.control.exon.coverage, norm.5SSA.coverage/norm.control.5SSA.coverage)
distance <- c(dist3SS.reg.vs.control, dist.exon.reg.vs.control, dist5SS.reg.vs.control)
table1 = data.frame(reg.exons.cov, control.exons.cov, enrichment, distance) 
rownames(table1) <- c("3ss","exon","5ss")
colnames(table1) <- c("reg.exons.cov","control.exons.cov","enrichment","distance")

# Table 2
ratio.regulated <- c(norm.3SSA.coverage/norm.5SSA.coverage, norm.exon.coverage/norm.3SSA.coverage, norm.exon.coverage/norm.5SSA.coverage, norm.5SSA.coverage/norm.3SSA.coverage)
ratio.contols <- c(norm.control.3SSA.coverage/norm.control.5SSA.coverage, norm.control.exon.coverage/norm.control.3SSA.coverage, norm.control.exon.coverage/norm.control.5SSA.coverage, norm.control.5SSA.coverage/norm.control.3SSA.coverage)
table2 = data.frame(ratio.regulated, ratio.contols) 
rownames(table2) <- c("3ss vs. 5ss","exon vs. 3ss", "exon vs. 5ss","5ss vs. 3ss")
colnames(table2) <- c("regulated","control")

library(gridExtra)
gt1 <- tableGrob(round(table1,2))
gt2 <- tableGrob(round(table2,2))
grid.arrange(gt1, gt2)
