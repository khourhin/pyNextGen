library(grid)
library(ggplot2)
library(gridExtra)
library(data.table)

count.files = list.files('/home/ekornobis/analysis/allemand/v1/bam_clustering_total_exons/counts',
                         full.names=T)
count.list = sapply(count.files, function(x) read.table(x, row.names=1), simplify=FALSE, USE.NAMES=TRUE)

# Concatenating all data in one dataframe
all.counts = rbindlist(count.list)

## Delete the reads which didn't align (not giving any exon counts)
all.counts.no.0 = all.counts[rowSums(all.counts)!=0,]

# Get the exon information (NOT USED HERE YET)
bedfile = read.table("~/analysis/allemand/v1/bam_clustering/final_consensus_7_long_Aligned/potential_exons.bed")
exons_length = bedfile[,3] - bedfile[,2] + 1

isoforms = unique(all.counts.no.0)

# Vector -> String, then count the frequencies of the isoforms and, after, sort them according to frequencies
isoforms.freq = table(apply(all.counts.no.0, 1, function(x) paste(x, collapse="")))

## Filtering out by minimum a read coverage
#isoforms.freq = isoforms.freq[isoforms.freq > 1]

isoforms.freq = isoforms.freq[order(isoforms.freq, decreasing=F)]
isoforms.freq = as.data.frame(cbind.data.frame(freqs=isoforms.freq,
                                               ids=paste("GPHN-", length(isoforms.freq):1, sep="")))
# as.data.frame is better apparently (not changing to numeric to "levels" for ex )

head(isoforms.freq)
# Prepare data for ggplot
heatDat = strsplit(rownames(isoforms.freq), "")
heatDat = as.data.frame(heatDat)
rownames(heatDat) = paste("E", 1:nrow(heatDat), sep="")
colnames(heatDat) = isoforms.freq[,2]
heatDat = melt(t(heatDat))

## heatDat$value = heatDat$value + rep(c(0,1),length(heatDat$value)/2)
head(heatDat)

# Presence/Absence tiles
colors <- c("white", "lightblue","lightgrey", "blue")
p1=ggplot(heatDat, aes(Var2, Var1, fill=factor(value)))  +
    geom_tile(colour="grey50") +
    scale_fill_manual(values=colors) +
    guides(fill=FALSE) +
    xlab("Exons") + ylab("Isoforms") + 
    theme(axis.text.y = element_text(size=1),
          axis.text.x = element_text(size=7),
          axis.ticks.y = element_blank(),
          plot.margin=unit(c(1,0,1,0), "cm"))

# Frequencies barplot
p2 = ggplot(isoforms.freq, aes(ids, freqs)) +
    geom_bar(stat="identity") +
    coord_flip() + xlim(levels(heatDat$Var1)[order(isoforms.freq$freqs)]) + # This xlim is wierd and gonna bite at some point
    xlab("") + ylab("Frequencies (log scale)") + 
    scale_y_log10() +
    scale_y_continuous(trans='log10',
                       breaks = c(1,2,5,10,25,50,100,200,500,1000,2000,5000,10000,20000)) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin=unit(c(1,0,1,0), "cm"))

grid.arrange(p1,p2)

grid.newpage()
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))


