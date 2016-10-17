library(grid)
library(ggplot2)
library(gridExtra)
library(reshape2)

# exonCountMat = read.table("~/analysis/allemand/counting/count_matrix", header=TRUE)
exonCountMat = read.table("~/analysis/allemand/v1/bam_clustering/final_consensus_7_long_Aligned/counts", header=FALSE, row.names=1)
colnames(exonCountMat)=paste('E', 1:ncol(exonCountMat), sep='')

# Heatmap from a matrix
# guides(fill=FALSE) to remove legend
z = as.data.frame(exonCountMat)
isoformMat = t(aggregate(z, by=z, length)[1:(ncol(z)+1)])
nIso = ncol(isoformMat)
colnames(isoformMat) = paste(rep("rawIso", nIso), 1:nIso, sep="_")

freq = isoformMat[nrow(isoformMat),]
isoforms=colnames(isoformMat)

isoFreq = data.frame(freq, factor(isoforms, levels=isoforms[order(freq)]))

heatDat = melt(isoformMat[-nrow(isoformMat),])
colors <- c("white", "lightblue")
p1=ggplot(heatDat, aes(Var1, Var2, fill=factor(value)))  + geom_tile(colour="grey50") + scale_fill_manual(values=colors) + guides(fill=FALSE) + ylim(levels(heatDat$Var2)[order(isoFreq$freq)]) # + theme(axis.text.x = element_text(angle = 45))  #+ theme(axis.ticks = element_blank(), axis.text.y = element_blank())

#iso_freq=(c(rep("iso1", 4), rep("iso2", 2), rep("iso3", 12), rep("iso4", 1), rep("iso5", 9), rep("iso6", 2)))
p2 = ggplot(isoFreq, aes(isoforms, freq)) + geom_bar(stat="identity") + coord_flip() + xlim(levels(heatDat$Var2)[order(isoFreq$freq)])
#p2 = qplot(iso_freq, geom="bar", xlab = "",) + coord_flip() + theme(axis.ticks = element_blank())

grid.arrange(p1,p2)

grid.newpage()
p3 = grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))

