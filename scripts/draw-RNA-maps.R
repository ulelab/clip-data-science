library("ggplot2")
library("smoother")
library("cowplot")

args<-commandArgs(TRUE)
smoothing_window <- 25

# RNA maps
# enhanced 3' splice sites
######################################
#enhanced.3SS <- read.table("PTB.enhanced.cassetteexons.bed-3SS-flanked300-selected-crosslink-coverage.bed", sep='\t')
enhanced.3SS <- read.table(args[1], sep='\t')

enhanced.3SS.plus <- enhanced.3SS[which(enhanced.3SS$V6 == "+"),]
enhanced.3SS.minus <- enhanced.3SS[which(enhanced.3SS$V6 == "-"),]

enhanced.3SS.plus.sum <- aggregate(enhanced.3SS.plus$V8, list(enhanced.3SS.plus$V7), FUN=sum)
enhanced.3SS.minus.sum <- aggregate(enhanced.3SS.minus$V8, list(enhanced.3SS.minus$V7), FUN=sum)
enhanced.3SS.minus.sum$Group.1 <- c(601:1)

enhanced.3SS.sum <- merge(enhanced.3SS.plus.sum, enhanced.3SS.minus.sum, by="Group.1")
enhanced.3SS.sum$x <- enhanced.3SS.sum$x.x + enhanced.3SS.sum$x.y

enhanced.3SS.sum$norm <- enhanced.3SS.sum$x / nrow(unique(enhanced.3SS[c("V5")]))
enhanced.3SS.sum$smooth <- smth(enhanced.3SS.sum$norm, window = smoothing_window, method = "gaussian") #SMOOTHING
enhanced.3SS.sum$map <- c(-300:300)
remove(enhanced.3SS)

# enhanced 3' splice sites
######################################
#enhanced.5SS <- read.table("PTB.enhanced.cassetteexons.bed-5SS-flanked300-selected-crosslink-coverage.bed", sep='\t')
enhanced.5SS <- read.table(args[2], sep='\t')

enhanced.5SS.plus <- enhanced.5SS[which(enhanced.5SS$V6 == "+"),]
enhanced.5SS.minus <- enhanced.5SS[which(enhanced.5SS$V6 == "-"),]

enhanced.5SS.plus.sum <- aggregate(enhanced.5SS.plus$V8, list(enhanced.5SS.plus$V7), FUN=sum)
enhanced.5SS.minus.sum <- aggregate(enhanced.5SS.minus$V8, list(enhanced.5SS.minus$V7), FUN=sum)
enhanced.5SS.minus.sum$Group.1 <- c(601:1)

enhanced.5SS.sum <- merge(enhanced.5SS.plus.sum, enhanced.5SS.minus.sum, by="Group.1")
enhanced.5SS.sum$x <- enhanced.5SS.sum$x.x + enhanced.5SS.sum$x.y

enhanced.5SS.sum$norm <- enhanced.5SS.sum$x / nrow(unique(enhanced.5SS[c("V5")]))
enhanced.5SS.sum$smooth <- smth(enhanced.5SS.sum$norm, window = smoothing_window, method = "gaussian") #SMOOTHING
enhanced.5SS.sum$map <- c(-300:300)
remove(enhanced.5SS)

# silenced 3' splice sites
######################################
#silenced.3SS <- read.table("PTB.silenced.cassetteexons.bed-3SS-flanked300-selected-crosslink-coverage.bed", sep='\t')
silenced.3SS <- read.table(args[3], sep='\t')

silenced.3SS.plus <- silenced.3SS[which(silenced.3SS$V6 == "+"),]
silenced.3SS.minus <- silenced.3SS[which(silenced.3SS$V6 == "-"),]

silenced.3SS.plus.sum <- aggregate(silenced.3SS.plus$V8, list(silenced.3SS.plus$V7), FUN=sum)
silenced.3SS.minus.sum <- aggregate(silenced.3SS.minus$V8, list(silenced.3SS.minus$V7), FUN=sum)
silenced.3SS.minus.sum$Group.1 <- c(601:1)

silenced.3SS.sum <- merge(silenced.3SS.plus.sum, silenced.3SS.minus.sum, by="Group.1")
silenced.3SS.sum$x <- silenced.3SS.sum$x.x + silenced.3SS.sum$x.y

silenced.3SS.sum$norm <- silenced.3SS.sum$x / nrow(unique(silenced.3SS[c("V5")]))
silenced.3SS.sum$smooth <- smth(silenced.3SS.sum$norm, window = smoothing_window, method = "gaussian") #SMOOTHING
silenced.3SS.sum$map <- c(-300:300)
remove(silenced.3SS)

# silenced 3' splice sites
######################################
#silenced.5SS <- read.table("PTB.silenced.cassetteexons.bed-5SS-flanked300-selected-crosslink-coverage.bed", sep='\t')
silenced.5SS <- read.table(args[4], sep='\t')

silenced.5SS.plus <- silenced.5SS[which(silenced.5SS$V6 == "+"),]
silenced.5SS.minus <- silenced.5SS[which(silenced.5SS$V6 == "-"),]

silenced.5SS.plus.sum <- aggregate(silenced.5SS.plus$V8, list(silenced.5SS.plus$V7), FUN=sum)
silenced.5SS.minus.sum <- aggregate(silenced.5SS.minus$V8, list(silenced.5SS.minus$V7), FUN=sum)
silenced.5SS.minus.sum$Group.1 <- c(601:1)

silenced.5SS.sum <- merge(silenced.5SS.plus.sum, silenced.5SS.minus.sum, by="Group.1")
silenced.5SS.sum$x <- silenced.5SS.sum$x.x + silenced.5SS.sum$x.y

silenced.5SS.sum$norm <- silenced.5SS.sum$x / nrow(unique(silenced.5SS[c("V5")]))
silenced.5SS.sum$smooth <- smth(silenced.5SS.sum$norm, window = smoothing_window, method = "gaussian") #SMOOTHING
silenced.5SS.sum$map <- c(-300:300)
remove(silenced.5SS)

# control 3' splice sites
######################################
#control.3SS <- read.table("PTB.control.cassetteexons.bed-3SS-flanked300-selected-crosslink-coverage.bed", sep='\t')
control.3SS <- read.table(args[5], sep='\t')

control.3SS.plus <- control.3SS[which(control.3SS$V6 == "+"),]
control.3SS.minus <- control.3SS[which(control.3SS$V6 == "-"),]

control.3SS.plus.sum <- aggregate(control.3SS.plus$V8, list(control.3SS.plus$V7), FUN=sum)
control.3SS.minus.sum <- aggregate(control.3SS.minus$V8, list(control.3SS.minus$V7), FUN=sum)
control.3SS.minus.sum$Group.1 <- c(601:1)

control.3SS.sum <- merge(control.3SS.plus.sum, control.3SS.minus.sum, by="Group.1")
control.3SS.sum$x <- control.3SS.sum$x.x + control.3SS.sum$x.y

control.3SS.sum$norm <- control.3SS.sum$x / nrow(unique(control.3SS[c("V5")]))
control.3SS.sum$smooth <- smth(control.3SS.sum$norm, window = smoothing_window, method = "gaussian") #SMOOTHING
control.3SS.sum$map <- c(-300:300)
remove(control.3SS)

# control 3' splice sites
######################################
#control.5SS <- read.table("PTB.control.cassetteexons.bed-5SS-flanked300-selected-crosslink-coverage.bed", sep='\t')
control.5SS <- read.table(args[6], sep='\t')

control.5SS.plus <- control.5SS[which(control.5SS$V6 == "+"),]
control.5SS.minus <- control.5SS[which(control.5SS$V6 == "-"),]

control.5SS.plus.sum <- aggregate(control.5SS.plus$V8, list(control.5SS.plus$V7), FUN=sum)
control.5SS.minus.sum <- aggregate(control.5SS.minus$V8, list(control.5SS.minus$V7), FUN=sum)
control.5SS.minus.sum$Group.1 <- c(601:1)

control.5SS.sum <- merge(control.5SS.plus.sum, control.5SS.minus.sum, by="Group.1")
control.5SS.sum$x <- control.5SS.sum$x.x + control.5SS.sum$x.y

control.5SS.sum$norm <- control.5SS.sum$x / nrow(unique(control.5SS[c("V5")]))
control.5SS.sum$smooth <- smth(control.5SS.sum$norm, window = smoothing_window, method = "gaussian") #SMOOTHING
control.5SS.sum$map <- c(-300:300)
remove(control.5SS)

tans <- 1.0
g_size <- 0.4
adj <- 0.4
# ggCAGE <- ggplot() + theme_bw() + 
#   geom_line(aes(enhanced.3SS.sum$map, as.vector(enhanced.3SS.sum$norm), colour="rep1"),size=g_size, alpha=tans) + 
#   geom_line(aes(enhanced.5SS.sum$map, as.vector(enhanced.5SS.sum$norm), colour="rep2"),size=g_size, alpha=tans) + 
#   ggtitle("RNA map") + 
#   xlab("position relative to 5'SS") + 
#   ylab("normalised coverage of cDNA-starts") + 
#   theme(text=element_text(size=8),axis.text=element_text(size=8), axis.title=element_text(size=8,face="plain")) + 
#   scale_colour_manual(values=c("#08519c", "#47B6D4")) +
#   scale_x_continuous(limits = c(-75, 25)) 

# find max for y-axis limit
ylim <- max(c(enhanced.3SS.sum$smooth, silenced.3SS.sum$smooth, control.3SS.sum$smooth, enhanced.5SS.sum$smooth, silenced.5SS.sum$smooth, control.5SS.sum$smooth), na.rm = TRUE)

# calculate the enrichment: regulated(enhanced or silenced) / control
enhanced.3SS.sum$smooth.enrichment <- enhanced.3SS.sum$smooth / control.3SS.sum$smooth
enhanced.5SS.sum$smooth.enrichment <- enhanced.5SS.sum$smooth / control.5SS.sum$smooth
silenced.3SS.sum$smooth.enrichment <- silenced.3SS.sum$smooth / control.3SS.sum$smooth
silenced.5SS.sum$smooth.enrichment <- silenced.5SS.sum$smooth / control.5SS.sum$smooth

write.table(enhanced.3SS.sum, paste(args[7],".enhanced.3SS.tab", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(enhanced.3SS.sum, paste(args[7],".enhanced.5SS.tab", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(enhanced.3SS.sum, paste(args[7],".silenced.3SS.tab", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(enhanced.3SS.sum, paste(args[7],".silenced.5SS.tab", sep=""), sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#ggRNAmap.3ss.smooth <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) +
ggRNAmap.3ss.smooth <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(enhanced.3SS.sum$map, as.vector(enhanced.3SS.sum$smooth), colour="enhanced"),size=g_size, alpha=tans) + 
  geom_line(aes(silenced.3SS.sum$map, as.vector(silenced.3SS.sum$smooth), colour="silenced"),size=g_size, alpha=tans) + 
  geom_line(aes(control.3SS.sum$map, as.vector(control.3SS.sum$smooth), colour="control"),size=g_size, alpha=tans) + 
  geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
  ggtitle(args[7]) + 
  xlab("position relative to 3'SS") + 
  ylab("normalised coverage of CL-clusters") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_colour_manual(values=c("#788285","#C90000","#0952B0")) +
  scale_x_continuous(limits = c(-300, 50)) +
  scale_y_continuous(limits = c(0, ylim))

#ggRNAmap.5ss.smooth <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) +
ggRNAmap.5ss.smooth <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(enhanced.5SS.sum$map, as.vector(enhanced.5SS.sum$smooth), colour="enhanced"),size=g_size, alpha=tans) + 
  geom_line(aes(silenced.5SS.sum$map, as.vector(silenced.5SS.sum$smooth), colour="silenced"),size=g_size, alpha=tans) + 
  geom_line(aes(control.5SS.sum$map, as.vector(control.5SS.sum$smooth), colour="control"),size=g_size, alpha=tans) + 
  geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
  ggtitle(args[7]) + 
  xlab("position relative to 5'SS") + 
  ylab("normalised coverage of CL-clusters") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_colour_manual(values=c("#788285","#C90000","#0952B0")) +
  scale_x_continuous(limits = c(-50, 300)) +
  scale_y_continuous(limits = c(0, ylim))

########################
# enrichment RNA-maps

ylim2 <- max(c(enhanced.3SS.sum$smooth.enrichment, silenced.3SS.sum$smooth.enrichment, enhanced.5SS.sum$smooth.enrichment, silenced.5SS.sum$smooth.enrichment), na.rm = TRUE)

ggRNAmap.3ss.smooth.enrichment <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(enhanced.3SS.sum$map, as.vector(enhanced.3SS.sum$smooth.enrichment), colour="enhanced"),size=g_size, alpha=tans) + 
  geom_line(aes(silenced.3SS.sum$map, as.vector(silenced.3SS.sum$smooth.enrichment), colour="silenced"),size=g_size, alpha=tans) + 
  geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
  ggtitle(args[7]) + 
  xlab("position relative to 3'SS") + 
  ylab("normalised enrichment of CL-clusters") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_colour_manual(values=c("#C90000","#0952B0")) +
  scale_x_continuous(limits = c(-300, 50)) +
  scale_y_continuous(limits = c(0, ylim2))

ggRNAmap.5ss.smooth.enrichment <- ggplot() + theme_cowplot() + theme(plot.title = element_text(size = 6, face = "plain")) + background_grid(major = "xy", minor = "none") +
  geom_line(aes(enhanced.5SS.sum$map, as.vector(enhanced.5SS.sum$smooth.enrichment), colour="enhanced"),size=g_size, alpha=tans) + 
  geom_line(aes(silenced.5SS.sum$map, as.vector(silenced.5SS.sum$smooth.enrichment), colour="silenced"),size=g_size, alpha=tans) + 
  geom_vline(xintercept = 0, alpha = 0.5, linetype="longdash") +
  ggtitle(args[7]) + 
  xlab("position relative to 5'SS") + 
  ylab("normalised enrichment of CL-clusters") + 
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain")) + 
  scale_colour_manual(values=c("#C90000","#0952B0")) +
  scale_x_continuous(limits = c(-50, 300)) +
  scale_y_continuous(limits = c(0, ylim2))


### plot ###
pdf(paste(args[7],".pdf", sep=""))

# Multiple plot function
#
# ggplot objects can be passed in .., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(ggRNAmap.3ss.smooth, ggRNAmap.3ss.smooth.enrichment, FALSE, ggRNAmap.5ss.smooth, ggRNAmap.5ss.smooth.enrichment, FALSE, cols=2)


